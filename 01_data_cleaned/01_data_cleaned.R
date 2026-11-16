# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    renv::restore()
    #install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

# functions to filter birds
# Assigns an administrative unit + GlobalNorth to each record
atribuir_geografia <- function(df, shapes) {
  
  valid <- !is.na(df$decimalLatitude) & !is.na(df$decimalLongitude) &
    df$decimalLatitude  >= -90  & df$decimalLatitude  <= 90 &
    df$decimalLongitude >= -180 & df$decimalLongitude <= 180
  
  # Padrão: sem coordenada válida
  df$name_en     <- "unassigned_coord"
  df$GlobalNorth <- "unassigned"
  
  coords <- df[valid, ] %>%
    distinct(decimalLongitude, decimalLatitude)
  
  pts <- st_as_sf(coords,
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326, remove = FALSE) %>%
    st_transform(st_crs(shapes))
  
  hits   <- st_intersects(pts, shapes)
  n_hits <- lengths(hits)
  idx    <- vapply(hits, function(h) if (length(h) == 1) h else NA_integer_,
                   integer(1))
  
  coords$name_en <- ifelse(n_hits == 1, shapes$name_en[idx],
                           ifelse(n_hits == 0, "unassigned_no_polygon",
                                  "unassigned_overlap"))
  coords$GlobalNorth <- ifelse(n_hits == 1,
                               as.character(shapes$GlobalNorth[idx]),
                               "unassigned")
  
  # Junta de volta (apenas nos registros com coordenada válida)
  geo <- df[valid, c("decimalLongitude", "decimalLatitude")] %>%
    left_join(coords, by = c("decimalLongitude", "decimalLatitude"))
  
  df$name_en[valid]     <- geo$name_en
  df$GlobalNorth[valid] <- geo$GlobalNorth
  df
}

# Marks the loss stage (only for those still alive).
# An 'NA' in 'keep' counts as a removal, just like in filter().
marcar <- function(df, manter, etapa) {
  remover <- is.na(df$lost_at) & !coalesce(manter, FALSE)
  df$lost_at[remover] <- etapa
  df
}

processar_ordem <- function(ordem_dir, parquet_dir, parquet_clean,
                            dir_stats, shapes_adj, arrow_threads = NULL) {
  
  ordem_nome <- sub(
    paste0(".*", nivel, "=([^/]+).*"),
    "\\1",
    ordem_dir
  )
  arq_rds    <- file.path(dir_stats, paste0(ordem_nome, ".rds"))
  
  # Recap: the .rds file only exists if the Parquet file has already been written
  if (file.exists(arq_rds)) {
    return(list(order = ordem_nome, ok = TRUE, pulou = TRUE, msg = NA_character_))
  }
  
  res <- tryCatch({
    sf::sf_use_s2(FALSE)
    if (!is.null(arrow_threads)) {
      arrow::set_cpu_count(arrow_threads)
      arrow::set_io_thread_count(arrow_threads)
    }
    
    aves_all <- open_dataset(file.path(parquet_dir, ordem_dir)) %>%
      collect() %>%
      mutate(.rid = row_number(),
             lost_at = NA_character_,
             origin_of_data = "gbif",
             eventDate = na_if(eventDate, ""))
    
    n_original <- nrow(aves_all)
    
    # Geography BEFORE any exclusion (one time only)
    aves_all <- atribuir_geografia(aves_all, shapes_adj)
    
    # Steps
    aves_all <- marcar(aves_all,
                       !is.na(aves_all$decimalLatitude) &
                         !is.na(aves_all$decimalLongitude) &
                         aves_all$decimalLatitude  >= -90  & aves_all$decimalLatitude  <= 90 &
                         aves_all$decimalLongitude >= -180 & aves_all$decimalLongitude <= 180,
                       "01_NA_lat_lon")
    
    aves_all <- marcar(aves_all,
                       !aves_all$name_en %in% c("unassigned_no_polygon",
                                                "unassigned_overlap"),
                       "02_sem_localizacao")
    
    aves_all <- marcar(aves_all, !is.na(aves_all$year), "03_sem_ano")
    
    aves_all <- marcar(aves_all,
                       aves_all$taxonRank %in% c("SPECIES", "SUBSPECIES") &
                         aves_all$occurrenceStatus != "ABSENT",
                       "04_taxon_occurrenceStatus")
    
    aves_all$countryCode <- countrycode(
      aves_all$countryCode,
      origin = "iso2c", destination = "iso3c",
      custom_match = c("XK" = "XKX"))
    
    idx_vivos <- which(is.na(aves_all$lost_at))
    if (length(idx_vivos) > 0) {
      flags <- clean_coordinates(
        x = aves_all[idx_vivos, ],
        lon = "decimalLongitude", lat = "decimalLatitude",
        countries = "countryCode", species = "species",
        tests = c("equal", "gbif", "institutions", "seas", "zeros")
      )
      aves_all$lost_at[idx_vivos[!flags$.summary]] <- "05_coordinateCleaner"
    }
    
    aves_all <- marcar(aves_all,
                       aves_all$coordinateUncertaintyInMeters / 1000 <= 10 |
                         is.na(aves_all$coordinateUncertaintyInMeters),
                       "06_incerteza_coord")
    
    aves_all <- marcar(aves_all, aves_all$year < 2026, "07_year_2026")
    
    # rounding, sorting, and duplicates
    aves_clean <- aves_all[is.na(aves_all$lost_at), ] %>%
      mutate(
        decimalLatitude  = round(decimalLatitude, 2),
        decimalLongitude = round(decimalLongitude, 2),
        month = ifelse(month > 12 | month < 1, NA, month),
        day   = ifelse(day > 31 | day < 1, NA, day)
      ) %>%
      arrange(speciesKey, year, month, day)
    
    dup_ok <- if (nrow(aves_clean) > 0) {
      cc_dupl(aves_clean,
              lon = "decimalLongitude", lat = "decimalLatitude",
              species = "speciesKey", value = "flagged")
    } else logical(0)
    
    # .rid = original position in aves_all (which was never reordered)
    aves_all$lost_at[aves_clean$.rid[!dup_ok]] <- "08_duplicatas"
    aves_clean <- aves_clean[dup_ok, ]
    n_final <- nrow(aves_clean)
    
    # Losses table for this order + checks
    perdas_ordem <- aves_all %>%
      mutate(etapa = coalesce(lost_at, "99_retido")) %>%
      count(etapa, GlobalNorth, name_en, name = "n") %>%
      mutate(order = ordem_nome, .before = 1)
    
    stopifnot(
      sum(perdas_ordem$n) == n_original,
      sum(perdas_ordem$n[perdas_ordem$etapa == "99_retido"]) == n_final,
      !any(perdas_ordem$GlobalNorth == "unassigned" &
             !perdas_ordem$etapa %in% c("01_NA_lat_lon", "02_sem_localizacao"),
           na.rm = TRUE)
    )
    
    # Save Parquet (without auxiliary columns) and, finally, the checkpoint
    aves_out <- aves_clean %>% select(-.rid, -lost_at)
    write_parquet(aves_out,
                  sink = file.path(parquet_clean,
                                   paste0("clean_", ordem_nome, ".parquet")))
    
    estat <- data.frame(order = ordem_nome,
                        registros_originais = n_original,
                        registros_clean = n_final)
    saveRDS(list(estat = estat, perdas = perdas_ordem), arq_rds)
    
    list(order = ordem_nome, ok = TRUE, pulou = FALSE, msg = NA_character_)
    
  }, error = function(e) {
    msg <- conditionMessage(e)
    cat(sprintf("[%s] %s: %s\n", Sys.time(), ordem_nome, msg),
        file = file.path(dir_stats, "erros_log.txt"), append = TRUE)
    list(order = ordem_nome, ok = FALSE, pulou = FALSE, msg = msg)
  })
  
  # release the order's memory before returning
  rm(list = setdiff(ls(), "res"))
  invisible(gc())
  res
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "cowplot",
  "tidyverse",
  "CoordinateCleaner",
  "sf",
  "countrycode",
  "arrow"
)

# executing install & load for each package
# TODO adaptar para uso do ambiente virtual renv
sapply(package_vec, install.load.package)

# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# FILTER AND CLEAN ----
## GBIF ----
# load data
rm(list = setdiff(ls(), c("local_directory",
                          "atribuir_geografia",
                          "marcar",
                          "processar_ordem"))); gc()
load(file.path(
  local_directory,
  "00_raw_data",
  "tetrapodstraits_data.RData")
)

# key of each occurrence
anyDuplicated(data_tetrapodstraits$gbifID) 
nrow(data_tetrapodstraits) 
# 24.591.154 Global species occurrences
table(data_tetrapodstraits$class)

# how much data is na?
data_tetrapodstraits %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year)),
            missing_date = sum(is.na(eventDate)))

# exists a difference between "year" and "event date" 
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
data_tetrapods_clean <- data_tetrapodstraits %>%
  drop_na(decimalLatitude,
          decimalLongitude) %>%
  select(-eventDate)

data_tetrapods_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)), 
            missing_long = sum(is.na(decimalLongitude)), 
            missing_day = sum(is.na(day)), 
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year))) 

nrow(data_tetrapods_clean) # 15489662 Global species occurrences
table(data_tetrapods_clean$class) # occurrences per classes 

# Extract geolocation
# shapefile to crop adm unit
load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData"))

# Verificar validade das geometrias
sf_use_s2(FALSE)
st_is_valid(geographic_shape_data) %>% table()

# extracting geolocation per specie
data_occurences_gbif <- st_as_sf(data_tetrapods_clean,
                            coords = c("decimalLongitude", "decimalLatitude"),
                            crs = 4326) %>%    # WGS84
  st_transform(st_crs(geographic_shape_data))  # equal area

data_occurences_gbif <- data_occurences_gbif%>%
  mutate(ID = 1:nrow(data_occurences_gbif))
any(duplicated(data_occurences_gbif$ID))

geographic_shape_data_adj <- geographic_shape_data %>%
  select(name_en, GlobalNorth, geometry)

data_occurences_geo <- st_join(data_occurences_gbif,
                               geographic_shape_data_adj, 
                               join = st_intersects,
                               left = TRUE) %>% 
  filter(!is.na(name_en)) %>% # remove NAs
  group_by(ID) %>%
  filter(n() == 1) %>% # remove records in overlapping polygons
  ungroup()

any(duplicated(data_occurences_geo$ID))
# 907,533 occurrence records outside of polygons (ocean, disputed areas, etc.)
nrow(data_occurences_gbif)-nrow(data_occurences_geo) 
names(data_occurences_geo)
table(data_occurences_geo$GlobalNorth)

# Filters
data_tetrapods_filter <- data_occurences_geo %>%
  drop_na(year)

nrow(data_tetrapods_filter)
table(data_tetrapods_filter$GlobalNorth)

data_tetrapods_filter <- data_tetrapods_filter %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "gbif") %>%
  #mutate(eventDate = na_if(eventDate, "")) %>%
  # U. ricardinii appears without a class classification
  mutate(class = if_else(class == "", "Squamata", class))

nrow(data_tetrapods_filter) # 24.402.844
table(data_tetrapods_filter$GlobalNorth)

# coordsclean
data_tetrapods_filter <- data_tetrapods_filter %>%
  sf::st_transform(4326) %>%
  dplyr::mutate(
    decimalLongitude = sf::st_coordinates(geometry)[, 1],
    decimalLatitude  = sf::st_coordinates(geometry)[, 2]
  ) %>%
  sf::st_drop_geometry()

flags <- clean_coordinates(x = data_tetrapods_filter,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           #countries = "countryCode",
                           #species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))
# Flagged 678.415 of 14141511 records, EQ = 0.05.

#Exclude problematic records
data_tetrapods_filter_spatialpoints <- data_tetrapods_filter[
  flags$.summary,
] 

nrow(data_tetrapods_filter_spatialpoints) # 13389364 occ
table(data_tetrapods_filter_spatialpoints$GlobalNorth)

# filter coordinates uncertain > 100km
data_tetrapods_filtered <- data_tetrapods_filter_spatialpoints %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 10 | 
           is.na(coordinateUncertaintyInMeters)) 
nrow(data_tetrapods_filtered) # 13319194 
table(data_tetrapods_filtered$GlobalNorth)

## Human observation non-birds ----
rm(list = setdiff(ls(), c("local_directory",
                          "atribuir_geografia",
                          "marcar",
                          "processar_ordem",
                          "data_tetrapods_filtered",
                          "geographic_shape_data_adj"))); gc()

load(file.path(
  local_directory,
  "00_raw_data",
  "nonbirds_humanobservation_data.RData")
)

# key of each occurence
anyDuplicated(data_nonbirds$gbifID) 
nrow(data_nonbirds) # 38.702.092 Global species occurrences

# how much data is na?
data_nonbirds %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year)),
            missing_date = sum(is.na(eventDate)))

# exists a difference between "year" and "event date" of 7280 occurences
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
#View(data_tetrapods_clean %>%
#       filter(!is.na(eventDate) & is.na(year)))
data_nonbirds_filter <- data_nonbirds %>%
  drop_na(decimalLatitude,
          decimalLongitude) 

data_nonbirds_filter %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)), # without missing data =)
            missing_long = sum(is.na(decimalLongitude)), 
            missing_day = sum(is.na(day)), # 3.080.541 NA's
            missing_month = sum(is.na(month)), # 2.374.207 NA's
            missing_year = sum(is.na(year))) # without missing data =)

nrow(data_nonbirds_filter) # 38349941 Global species occurrences
table(data_nonbirds_filter$class) # occurences per classes 

table(data_nonbirds$class)

# spatial id
data_occurences_hononbirds <- st_as_sf(data_nonbirds_filter,
                                 coords = c("decimalLongitude", "decimalLatitude"),
                                 crs = 4326) %>%    # WGS84
  st_transform(st_crs(geographic_shape_data_adj))  # equal area

data_occurences_hononbirds <- data_occurences_hononbirds %>%
  mutate(ID = 1:nrow(data_occurences_hononbirds))
any(duplicated(data_occurences_hononbirds$ID))

# Without spatial id
data_occurences_geo <- st_join(data_occurences_hononbirds,
                               geographic_shape_data_adj, 
                               join = st_intersects,
                               left = TRUE) %>% 
  filter(!is.na(name_en)) %>% 
  group_by(ID) %>%
  filter(n() == 1) %>%
  ungroup()

any(duplicated(data_occurences_geo$ID))
nrow(data_occurences_geo)
table(data_occurences_geo$GlobalNorth)

# Without data
data_occurences_geo_date <- data_occurences_geo %>%
  drop_na(year) 

nrow(data_occurences_geo_date)
table(data_occurences_geo_date$GlobalNorth)

# filter
data_nonbirds_clean <- data_occurences_geo_date %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "HumanObservation") %>%
  mutate(eventDate = na_if(eventDate, "")) 

nrow(data_nonbirds_clean) # 38.317.234 Global species occurrences
table(data_nonbirds_clean$GlobalNorth)

# Coordsclean
data_nonbirds_clean <- data_nonbirds_clean %>%
  sf::st_transform(4326) %>%
  dplyr::mutate(
    decimalLongitude = sf::st_coordinates(geometry)[, 1],
    decimalLatitude  = sf::st_coordinates(geometry)[, 2]
  ) %>%
  sf::st_drop_geometry()

flags <- clean_coordinates(x = data_nonbirds_clean,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))

#Exclude problematic records
data_nonbirds_filter_spatialpoints <- data_nonbirds_clean[
  flags$.summary,
]
nrow(data_nonbirds_filter_spatialpoints) # 35.209.287 Global species occ
table(data_nonbirds_filter_spatialpoints$GlobalNorth)

# filter coordinates uncertain > 10km
data_nonbirds_filtered <- data_nonbirds_filter_spatialpoints %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 10 | 
           is.na(coordinateUncertaintyInMeters)) 

nrow(data_nonbirds_filtered) # 35058217 Global species occ
table(data_nonbirds_filtered$GlobalNorth)

## Human observation birds ----
#
# Data are large and were split into multiple .parquet files; 
# Passeriformes records were further partitioned by family.
#
# nested_cols <- c("speciesKey", "class", "order", 
#                 "family", "species", "gbifID",
#                 "day", "month", "year",
#                 "decimalLatitude", "decimalLongitude",
#                 "coordinateUncertaintyInMeters", "occurrenceStatus",
#                 "taxonRank", "countryCode", "eventDate")
#
#ebird_data <- arrow::open_delim_dataset(
#  file.path(temp_dir, "0018107-251025141854904.csv"), 
#  parse_options = csv_parse_options(
#    delimiter = "\t", 
#    newlines_in_values = TRUE),
#  read_options = list(
#    block_size = 1000000000L, 
#    use_threads = TRUE)) %>%
#  select(all_of(nested_cols))

# create multiple .parquet files
#ebird_data |> 
#  group_by(order) |> 
#  write_dataset(path = parquet_dir,
#                format = "parquet")

#tibble(
#  arquivos = list.files(parquet_dir, recursive = TRUE),
#  tamanho_MB = file.size(file.path(parquet_dir, arquivos)) / 1024^2
#) %>% View()

nivel <- "order"  # select whether to process 'order' or 'family'

if (nivel == "order") {
  
  # Birds orders non-passeriformes
  parquet_dir <- file.path(local_directory, 
                           "00_raw_data",
                           "parquet_por_ordem")
  
  ordens_disponiveis <- list.dirs(
    parquet_dir,
    recursive = FALSE,
    full.names = FALSE
  )
  ordens_disponiveis <- ordens_disponiveis[
    ordens_disponiveis != "order=Passeriformes"
  ]
} else if (nivel == "family") {
  
  # Family Passeriformes
  parquet_dir <- file.path(local_directory, 
                           "00_raw_data",
                           "parquet_por_ordem",
                           "Order=Passeriformes",
                           "family")
  
  ordens_disponiveis <- list.dirs(
    parquet_dir,
    recursive = FALSE,
    full.names = FALSE
  )
  
} else {
  stop("Nível taxonômico inválido. Escolha 'order' ou 'family'.")
}

# dir stats of .rds
dir_stats <- file.path(
  local_directory,
  "01_data_cleaned",
  "parquet_clean",
  "clean_statistics",
  if (nivel == "order") "por_ordem" else "por_familia"
)

# save clean .parquet files
parquet_clean <- file.path(local_directory,
                           "01_data_cleaned",
                           "parquet_clean")

# shapes_adj, parquet_dir, parquet_clean, local_directory; sf_use_s2(FALSE)
dir.create(dir_stats,     recursive = TRUE, showWarnings = FALSE)
dir.create(parquet_clean, recursive = TRUE, showWarnings = FALSE)

# Recap: the .rds file only exists if the Parquet file has already been written
resultados <- vector("list", length(ordens_disponiveis))
for (i in seq_along(ordens_disponiveis)) {
  cat(sprintf("[%d/%d] %s\n", i, length(ordens_disponiveis), ordens_disponiveis[i]))
  resultados[[i]] <- processar_ordem(ordens_disponiveis[i], parquet_dir,
                                     parquet_clean, dir_stats,
                                     geographic_shape_data_adj)
}

# Falhas
falhas <- Filter(function(x) !isTRUE(x$ok), resultados)
if (length(falhas) > 0) print(bind_rows(lapply(falhas, as.data.frame)))

# Consolidar a partir dos checkpoints
todos <- lapply(list.files(dir_stats,
                           pattern = "\\.rds$", 
                           full.names = TRUE), readRDS)
estatisticas <- bind_rows(lapply(todos, `[[`, "estat"))
perdas       <- bind_rows(lapply(todos, `[[`, "perdas"))

# Conferência: algum táxon ficou sem checkpoint?
print(setdiff(
  sub(paste0("^", nivel, "="), "", ordens_disponiveis),
  estatisticas[[nivel]]
))

print(estatisticas)

# Salvar estatísticas e perdas
write.table(
  estatisticas,
  file.path(local_directory,
            "01_data_cleaned",
            "parquet_clean",
            "clean_statistics",
            paste0("clean_statistics_", nivel, ".txt"))
)
write.table(
  perdas,
  file.path(local_directory,
            "01_data_cleaned",
            "parquet_clean",
            "clean_statistics",
            paste0("clean_losses_", nivel, "_step_region.txt"))
)

# Resumo 1: perdas por etapa × GlobalNorth
perdas %>%
  group_by(etapa, GlobalNorth) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = GlobalNorth, values_from = n, values_fill = 0)

# Resumo 2: detalhe das etapas 01 e 02
perdas %>%
  filter(etapa %in% c("01_NA_lat_lon", "02_sem_localizacao")) %>%
  group_by(etapa, name_en) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  filter(name_en %in% c("unassigned_coord", "unassigned_no_polygon",
                        "unassigned_overlap"))

perdas_etapa <- perdas %>%
  group_by(etapa, GlobalNorth) %>%
  summarise(n = sum(n), .groups = "drop")

# Total original de cada grupo (perdidos + retidos)
total_g <- perdas_etapa %>%
  group_by(GlobalNorth) %>%
  summarise(total = sum(n), .groups = "drop")

# Retidos após cada etapa = total - perdas acumuladas
retidos <- perdas_etapa %>%
  filter(etapa != "99_retido") %>%
  tidyr::complete(etapa, GlobalNorth, fill = list(n = 0)) %>%
  arrange(GlobalNorth, etapa) %>%
  group_by(GlobalNorth) %>%
  mutate(perdidos_acum = cumsum(n)) %>%
  ungroup() %>%
  left_join(total_g, by = "GlobalNorth") %>%
  mutate(retidos = total - perdidos_acum) %>%
  select(etapa, GlobalNorth, retidos)

# Linha de partida (antes de qualquer filtro)
original <- total_g %>%
  transmute(etapa = "00_original", GlobalNorth, retidos = total)

tab_retidos <- bind_rows(original, retidos) %>%
  tidyr::pivot_wider(names_from = GlobalNorth, values_from = retidos) %>%
  mutate(total = rowSums(across(-etapa))) %>%
  arrange(etapa)

print(tab_retidos, n = Inf)

## BioTIME 2.0v ----
rm(list = setdiff(ls(), c("local_directory",
                          "validar_dados",
                          "comparative_plot",
                          "data_tetrapods_filtered",
                          "data_nonbirds_filtered",
                          "geographic_shape_data_adj"))); gc()
load(file.path(
  local_directory,
  "00_raw_data",
  "biotime_data.RData")
) 

data_biotime_clean <- biotime_data_key_precleaned %>%
  select("speciesKey", "class", "order", "family", "valid_name", 
         "LATITUDE", "LONGITUDE", "YEAR", "MONTH", "DAY") %>%
  rename(decimalLatitude = LATITUDE, 
         decimalLongitude = LONGITUDE,
         species = valid_name,
         day = DAY,
         month = MONTH,
         year = YEAR) %>%
  mutate(origin_of_data = "biotime")
nrow(data_biotime_clean)

# extracting geolocation per specie
data_occurences_biotime <- st_as_sf(data_biotime_clean,
                                 coords = c("decimalLongitude", "decimalLatitude"),
                                 crs = 4326) %>%    # WGS84
  st_transform(st_crs(geographic_shape_data_adj))  # equal area

data_occurences_biotime <- data_occurences_biotime %>%
  mutate(ID = 1:nrow(data_occurences_biotime))
any(duplicated(data_occurences_biotime$ID))

data_occurences_geo <- st_join(data_occurences_biotime,
                               geographic_shape_data_adj, 
                               join = st_intersects,
                               left = TRUE) %>% 
  filter(!is.na(name_en)) %>% 
  group_by(ID) %>%
  filter(n() == 1) %>% 
  ungroup()

any(duplicated(data_occurences_geo$ID))
nrow(data_occurences_biotime)-nrow(data_occurences_geo) 
nrow(data_occurences_geo)
table(data_occurences_geo$GlobalNorth)

data_biotime_clean %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_day = sum(is.na(day)),     # 863.686
            missing_month = sum(is.na(month)), # 766.064
            missing_year = sum(is.na(year)))   # without missing data

# Filter
data_nonbirds_filter <- data_occurences_geo %>%
  sf::st_transform(4326) %>%
  dplyr::mutate(
    decimalLongitude = sf::st_coordinates(geometry)[, 1],
    decimalLatitude  = sf::st_coordinates(geometry)[, 2]
  ) %>%
  sf::st_drop_geometry()

flags <- clean_coordinates(x = data_nonbirds_filter,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))

#Exclude problematic records
data_biotime_filtered <- data_nonbirds_filter[
  flags$.summary,
]
nrow(data_biotime_filtered)
table(data_biotime_filtered$GlobalNorth)

## speciesLink ----
rm(list = setdiff(ls(), c("local_directory",
                          "validar_dados",
                          "comparative_plot",
                          "data_tetrapods_filtered",
                          "data_nonbirds_filtered",
                          "data_biotime_filtered",
                          "geographic_shape_data_adj"))); gc()
load(file.path(
  local_directory,
  "00_raw_data",
  "splinks_data.RData")
)

# names(splink_data_key_precleaned)
data_splink_clean <- splink_data_key_precleaned %>%
  select("speciesKey", "scientificname","monthcollected","daycollected",
         "latitude", "longitude", "yearcollected", "coordinateprecision") %>%
  rename(decimalLatitude = latitude, 
         decimalLongitude = longitude,
         species = scientificname,
         day = daycollected,
         month = monthcollected,
         year = yearcollected) %>%
  mutate(origin_of_data = "splink") %>%
  mutate(
    day = as.numeric(na_if(trimws(as.character(day)), "")),
    month = as.numeric(na_if(trimws(as.character(month)), "")),
    year = as.numeric(na_if(trimws(as.character(year)), ""))
  ) %>%
  mutate(
    year = ifelse(year > 2025 | year < 1758, NA_real_, year),
    month = ifelse(month > 12 | month < 1, NA_real_, month),
    day = ifelse(day > 31 | day < 1, NA_real_, day)
  ) %>%
  mutate(
    coordinateprecision = na_if(
      trimws(as.character(coordinateprecision)), ""
    )
  )

data_splink_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),     # 438.648
            missing_month = sum(is.na(month)), # 473.340
            missing_year = sum(is.na(year)))
nrow(data_splink_clean)

# invalidity coordinates
data_splink_cleaned <- data_splink_clean %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    decimalLatitude >= -90 & decimalLatitude <= 90,
    decimalLongitude >= -180 & decimalLongitude <= 180
  ) 
nrow(data_splink_cleaned)

# extracting geolocation per specie
data_occurences_splink <- st_as_sf(data_splink_cleaned,
                                    coords = c("decimalLongitude", "decimalLatitude"),
                                    crs = 4326) %>%    # WGS84
  st_transform(st_crs(geographic_shape_data_adj))  # equal area

data_occurences_splink <- data_occurences_splink %>%
  mutate(ID = 1:nrow(data_occurences_splink))
any(duplicated(data_occurences_splink$ID))

data_occurences_geo <- st_join(data_occurences_splink,
                               geographic_shape_data_adj, 
                               join = st_intersects,
                               left = TRUE) %>% 
  filter(!is.na(name_en)) %>% 
  group_by(ID) %>%
  filter(n() == 1) %>% 
  ungroup()

any(duplicated(data_occurences_geo$ID))
nrow(data_occurences_splink)-nrow(data_occurences_geo) 
nrow(data_occurences_geo)
table(data_occurences_geo$GlobalNorth)

# Without date
data_occurences_geo_date <- data_occurences_geo %>%
  drop_na(year)
nrow(data_occurences_geo_date)
table(data_occurences_geo_date$GlobalNorth)

# Filter
data_splink_cleaned <- data_occurences_geo_date %>%
  sf::st_transform(4326) %>%
  dplyr::mutate(
    decimalLongitude = sf::st_coordinates(geometry)[, 1],
    decimalLatitude  = sf::st_coordinates(geometry)[, 2]
  ) %>%
  sf::st_drop_geometry()

flags <- clean_coordinates(x = data_splink_cleaned,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros")) 
# Flagged 244567 of 646315 records, EQ = 0.38.

#Exclude problematic records
data_splink_filter_spatialpoints <- data_splink_cleaned[
  flags$.summary,
]

nrow(data_splink_filter_spatialpoints)
table(data_splink_filter_spatialpoints$GlobalNorth)

# filter coordinates uncertain > 100km
# filter years from 1900
data_splink_filtered <- data_splink_filter_spatialpoints %>%
  mutate(coordinateprecision = as.numeric(coordinateprecision)) %>%
  filter(coordinateprecision / 1000 <= 10 | is.na(coordinateprecision))

nrow(data_splink_filter_spatialpoints)
table(data_splink_filter_spatialpoints$GlobalNorth)

# If data its ok, save data
save(
  data_tetrapods_filtered,
  data_nonbirds_filtered,
  data_biotime_filtered,
  data_splink_filtered, 
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "datasets_filtered.RData")
)

# UNIFYING OCCURENCES DATA ---- 
# load data
rm(list = setdiff(ls(), c("local_directory"))); gc()

# Retain only keys present in `species_list_tetrapods_filter`, which contains
# all species recorded in TetrapodTraits. Some species may have been excluded
# from `data_tetrapods_filtered` because their expert observations lack dates,
# even though dated records exist from human observations. For example,
# Cristidorsa otai must be retained because it has dated human observation
# records.
load(file.path(
  local_directory,
  "00_raw_data",
  "tetrapodstraits_data.RData")
)
rm(list = setdiff(ls(), c("local_directory",
                          "species_list_tetrapods_filter"
))); gc()

parquet_dir <- file.path(local_directory,
                         "01_data_cleaned",
                         "parquet_clean")
data_birds <- open_dataset(
  sources = parquet_dir,
  unify_schemas = TRUE
) %>%
  collect()

nrow(data_birds)
table(data_birds$GlobalNorth)

load(file.path(
  local_directory,
  "01_data_cleaned",
  "datasets_filtered.RData")
)

# BioTIME lost keys
chaves_perdidas_biotime <- setdiff(
  data_biotime_filtered$speciesKey,
  species_list_tetrapods_filter$speciesKey
) 
# Retain ONLY speciesKeys that are NOT included in the lost species list.
data_biotime_filtered <- data_biotime_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month))

# speciesLink lost keys
chaves_perdidas_splink <- setdiff(
  data_splink_filtered$speciesKey,
  species_list_tetrapods_filter$speciesKey
) 

data_splink_filtered <- data_splink_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month))

# Non-Birds HumanObsevation lost keys
chaves_perdidas_nonbirdsHO <- setdiff(
  data_nonbirds_filtered$speciesKey,
  species_list_tetrapods_filter$speciesKey
) 
# No speciesKeys are missing, as all records were downloaded based on speciesKey. 

data_nonbirds_filtered <- data_nonbirds_filtered %>% 
  #filter(!(speciesKey %in% chaves_perdidas_nonbirdsHO)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month)) %>%
  mutate(origin_of_data = "HUMAN_OBSERVATION")

# Birds HumanObsevation lost keys
chaves_perdidas_birdsHO <- setdiff(
  data_birds$speciesKey,
  species_list_tetrapods_filter$speciesKey
) 
# No speciesKeys are missing, as all records were downloaded based on speciesKey.

data_birds <- data_birds %>% 
  #filter(!(speciesKey %in% chaves_perdidas_birdsHO)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month),
         GlobalNorth = as.double(GlobalNorth)) %>%
  mutate(origin_of_data = "HUMAN_OBSERVATION")

data_occurences_precleaned <- bind_rows(
  data_tetrapods_filtered,
  data_biotime_filtered,
  data_splink_filtered,
  data_nonbirds_filtered,
  data_birds)

nested_cols <- c("speciesKey", "class", "order", "family", "species",
                 "gbifID", "day", "month", "year",
                 "decimalLatitude", "decimalLongitude", "origin_of_data", 
                 "name_en", "GlobalNorth")

data_occurences_precleaned <- data_occurences_precleaned %>% 
  select(all_of(nested_cols))

write_parquet(data_occurences_precleaned,  # Cleaned and integrated data
              sink = file.path(
                local_directory,
                "01_data_cleaned",
                "data_occurences_filtered.parquet"))

# Create a dataframe combining all vectors.
# TODO: Use these keys to download occurrence records for additional species
# and retrieve more polygons if needed.
chaves_perdidas_combinadas <- bind_rows(
  data.frame(speciesKey = chaves_perdidas_biotime, Database = "Biotime"),
  data.frame(speciesKey = chaves_perdidas_splink, Database = "splink"),
  #data.frame(speciesKey = chaves_perdidas_nonbirdsHO, Database = "nonbirdsHO"),
  #data.frame(speciesKey = chaves_perdidas_birdsHO, Database = "birdsHO")
) %>%
  arrange(speciesKey, Database)

# A total of 852 keys are missing from TetrapodTraits but are present in
# BioTime and/or speciesLink. These records can be checked to determine
# whether they should be included in our dataset.
chaves_perdidas_combinadas %>% 
  distinct(speciesKey) %>% nrow() 

write.table(chaves_perdidas_combinadas, file.path(
  local_directory, "chaves_perdidas.txt"))

# REMOVE DUPLICATES ----
rm(list = setdiff(ls(), c("local_directory"))); gc()

# occurences data
data_occurences_precleaned <- open_dataset(
  sources = file.path(
    local_directory,
    "01_data_cleaned",
    "data_occurences_filtered.parquet")
)

#data_occurences_precleaned %>%
#  count(origin_of_data) %>%
#  compute() %>% 
#  collect()

data_occurences_precleaned <- data_occurences_precleaned %>%
  mutate(
    # duas casas decimais = 1.11 km 
    decimalLatitude = round(decimalLatitude, 2), 
    decimalLongitude = round(decimalLongitude, 2),
    month = ifelse(month > 12 | month < 1, NA, month),
    day = ifelse(day > 31 | day < 1, NA, day)          
  ) %>%
  arrange(speciesKey, year, month, day) %>%
  collect()

nrow(data_occurences_precleaned)
table(data_occurences_precleaned$GlobalNorth)

# Remove duplicates
duplicated_flags <- cc_dupl(
  data_occurences_precleaned,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "speciesKey",
  value = "flagged"
)

data_occurences_precleaned_duplicate <- data_occurences_precleaned[
  duplicated_flags == TRUE, ]

nrow(data_occurences_precleaned_duplicate)
table(data_occurences_precleaned_duplicate$GlobalNorth)

length(unique(data_occurences_precleaned_duplicate$speciesKey)) 
# 28.455 global scale spp

# FILTER POINTS BY RANGE POLYGONS ----
# Shapefile base
load(file = file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_polygons_key.RData")) # From TetrapodTraits

# The 110 km and 220 km scenarios represent expansions of each species' original range,
# aiming to capture occurrences that may fall outside the mapped polygons due to
# spatial limitations in the original range boundaries.
load(file = file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_polygons_sensitivity.RData")
)

cenario <- "110km" # "original", "110km", "220km"
# shapefile
if (cenario == "base") {
  data_tetrapods_sa <- tetrapods_polygons_key
} else if (cenario == "110km") {
  data_tetrapods_sa <- tetrapods_polygons_sensitivity %>%
    filter(scenario == "buffer_110km")
} else if (cenario == "220km") {
  data_tetrapods_sa <- tetrapods_polygons_sensitivity %>%
    filter(scenario == "buffer_220km")
} else {
  stop("Cenário inválido. Escolha: original, 110km ou 220km.")
}

# occurence filter
data_wallacean_knownledge_sa <- data_occurences_precleaned_duplicate %>%
  filter(speciesKey %in% data_tetrapods_sa$speciesKey)
# polygons filter
data_tetrapods_sa <- data_tetrapods_sa %>%
  filter(speciesKey %in% data_occurences_precleaned_duplicate$speciesKey)

# Number of species
length(unique(data_wallacean_knownledge_sa$speciesKey))
length(unique(data_tetrapods_sa$speciesKey))

# geometry
data_tetrapods_sa_transf <- st_cast(
  data_tetrapods_sa,
  "MULTIPOLYGON"
)
geometry_types <- st_geometry_type(data_tetrapods_sa_transf)
table(geometry_types)

# check if list_occurences_clean is empty before run
unique_species_keys <- unique(data_wallacean_knownledge_sa$speciesKey)
list_occurences_clean <- data.frame()
species_keys_with_errors <- c()

for (species_key in unique_species_keys) {
  # Use tryCatch to catch errors
  tryCatch({
    # Select the polygon corresponding to the species
    especie_poligono <- data_tetrapods_sa_transf %>%
      filter(speciesKey == species_key)
    
    subset_data <- data_wallacean_knownledge_sa %>%
      filter(speciesKey == species_key)
    
    # Apply the cc_iucn function
    range_flags <- cc_iucn(
      x = subset_data,
      range = especie_poligono,
      species = "speciesKey",
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      value = "flagged"
    )
    
    # Filter the data using the values returned by range_flags
    subset_data_filtered <- subset_data[range_flags == TRUE, ]
    
    # Add the filtered data to the final list
    list_occurences_clean <- rbind(list_occurences_clean, subset_data_filtered)
    
  }, error = function(e) {
    # If an error occurs, save the species_key to the vector
    species_keys_with_errors <<- c(species_keys_with_errors, species_key)
    # Optional: Print the error or the species_key for debugging
    message(paste("Erro no speciesKey:", species_key, "-", e$message))
  })
}
# Check for speciesKeys that encountered errors
print("SpeciesKey com erro:")
print(species_keys_with_errors)
#end_time <- Sys.time()
#print(end_time - start_time) # Time difference of 1h56min

#data_tetrapods_sa_transf <- data_tetrapods_sa %>%
#    filter(speciesKey %in% species_keys_with_errors)

nrow(list_occurences_clean) # 124068850
nrow(data_wallacean_knownledge_sa)-nrow(list_occurences_clean) # 14926904 removed
length(unique(list_occurences_clean$speciesKey)) # 27402 spp 
length(
  unique(
    data_wallacean_knownledge_sa$speciesKey))-length(
      unique(list_occurences_clean$speciesKey)) # 1053 spp removed

table(list_occurences_clean$GlobalNorth)
list_species <- list_occurences_clean %>%
 distinct(speciesKey, .keep_all=TRUE)
table(list_species$class)

# cleaned and integrated data
write_parquet(
  list_occurences_clean,  
  sink = file.path(
    local_directory,
    "01_data_cleaned",
    paste0(
      "list_occurences_clean",
      if (cenario != "base") paste0("_", cenario),
      ".parquet"
    )
  )
)

# NESTED DATAFRAME WITH BIOLOGICAL TRAITS ----
rm(list = setdiff(ls(), c("local_directory"))); gc()

# occurences data
data_occurences_units <- open_dataset(
  sources = file.path(
    local_directory,
    "01_data_cleaned",
    "list_occurences_clean.parquet")
) %>% collect()

# To access sensitivity analysis
#load(file = file.path(
#  local_directory,
#  "01_data_cleaned",
#  "data_occurences_geolocation_220.RData")
#)

load(file.path(
  local_directory,
  "00_raw_data",
  "trait_data.RData")
) # match taxonomic and trait information with speciesKey

# Algumas especies presentes no TetrapodTraits sao consideradas sinonimos,
# entao precisamos manter apenas as espécies mais antigas, e lidar com os atributos delas
#trait_data %>%
#  arrange(speciesKey, YearOfDescription) %>%
# filter(duplicated(speciesKey) | duplicated(speciesKey, fromLast = TRUE)) %>%
#  View()
tetrapods_key_species <- trait_data %>%
  # Especies com o mesmo speciesKey fico com a mais antiga
  arrange(speciesKey, YearOfDescription) %>%
  select("speciesKey","scientificName","Class", "Order","Family", 
         "YearOfDescription","BodyLength_mm","ImputedLength","BodyMass_g",
         "ImputedMass","Diu","Noc","Nocturnality", "Fos", "Ter", "Aqu", "Arb", "Aer",
         "ImputedHabitat","Verticality","ImputedActTime","ImputedMajorHabitat","RangeSize",
         "HumanDensity", "AssessedStatus", "Latitude", "Elevation", "ETA50K", "AnnuMeanTemp") %>%
  distinct(speciesKey, .keep_all = TRUE) # 31594 spp

anyDuplicated(tetrapods_key_species$speciesKey) # ok

data_tetrapods_nested <- data_occurences_units %>%
  group_by(speciesKey) %>%  # Agrupa pelo speciesKey
  nest() %>%
  rename(event_table = data)
anyDuplicated(data_tetrapods_nested$speciesKey) # ok

data_wallacean_nested <- left_join(
  data_tetrapods_nested,
  tetrapods_key_species,
  by = "speciesKey"
) %>% 
  relocate(event_table, .after = Family) %>%
  arrange(Class, Order, scientificName)

# check data
table(data_wallacean_nested$Class)

anyDuplicated(data_wallacean_nested$speciesKey)

# unnest data
data_wallacean_unnested <- data_wallacean_nested %>%
  unnest(cols = c(event_table))
nrow(data_wallacean_unnested) 

save(data_wallacean_nested,
     data_wallacean_unnested,
     file = file.path(
       local_directory,
       "01_data_cleaned",
       "dataset_occurences.RData")
)

# SPECIES LIST PER ADMINISTRATIVE UNIT ----
rm(list = setdiff(ls(), c("local_directory"))); gc()
# occurences data
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "dataset_occurences.RData")
)

# polygon data
load(file = file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_polygons_key.RData")) # From TetrapodTraits

# shapefile to crop adm unit
load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData"))

# dataset antes de filtrar
data_occurences_precleaned <- open_dataset(
  sources = file.path(
    local_directory,
    "01_data_cleaned",
    "data_occurences_filtered.parquet")
) 

load(file.path(
  local_directory,
  "00_raw_data",
  "trait_data.RData")
) # match taxonomic and trait information with speciesKey

tetrapods_key_species <- trait_data %>%
  arrange(speciesKey, YearOfDescription) %>%
  dplyr::select("speciesKey","scientificName","Class", "Order","Family", 
         "YearOfDescription","BodyLength_mm","ImputedLength","BodyMass_g",
         "ImputedMass","Diu","Noc","Nocturnality", "Fos", "Ter", "Aqu", "Arb", "Aer",
         "ImputedHabitat","Verticality","MajorHabitatSum","ImputedMajorHabitat","RangeSize",
         "HumanDensity", "AssessedStatus", "Latitude", "Elevation", "ETA50K", "AnnuMeanTemp") %>%
  distinct(speciesKey, .keep_all = TRUE) 

# pegando todas as chaves que tem nos dados baixados e filtrados pelos poligonos
keys_with_occ <- data_wallacean_unnested %>%
  dplyr::select(speciesKey) %>%
  distinct() %>%
  pull()
length(keys_with_occ) # 27.395 especies com occ filtradas 

# pegando a chave das especies que tem poligono mas nao estao nas nossas ocorrencias
# aqui eu to pegando todas as especies que nao achei ocorrencia, 
# ou pq sairam devido as ocorrencias cairem fora do poligono de especialista
keys_without_occ <- tetrapods_polygons_key %>% # 31583 com poligono
  filter(!speciesKey %in% keys_with_occ) %>%
  st_drop_geometry() %>%
  # Pegando só as chaves speciesKey delas
  select(speciesKey) %>%
  distinct() %>%
  pull()
# 4188 especies com poligonos mas nao estao nos nossos dados
length(keys_without_occ) 

# Tem occ perdidas por mismatch espacial?
# as especies que tiverem em 'lost_occurences' 
# sao as especies que tem poligono, recuperamos ocorrencias, mas foram excluidas
# porque existe um mismatch entre as duas fontes de dados. O erro pode estar
# tanto nos pontos, quanto nos poligonos, mudanças taxonômicas entre as bases, etc
lost_occurences <- data_occurences_precleaned %>%
  filter(speciesKey %in% keys_without_occ) %>%
  select(speciesKey) %>%
  distinct() %>%
  collect() %>%
  pull()
length(lost_occurences)
# 1064 perdemos pq os pontos de ocorrencia cairam todos fora dos poligonos
# para essas especies nao podemos considerar no completeness pq tem dados delas
# no GBIF, mas por algum motivo, as ocorrencias delas tem mismatch espacial com
# os poligonos. 
count_loss_class <- left_join(data.frame(speciesKey = lost_occurences),
                              tetrapods_key_species, by = "speciesKey")
table(count_loss_class$Class)
# Quem sao as 3124 especies que temos poligono, mas nao temos ocorrencia
count_withoutocc_class <- left_join(
  data.frame(speciesKey = setdiff(keys_without_occ, lost_occurences)),
  tetrapods_key_species, by = "speciesKey")
table(count_withoutocc_class$Class)

tetrapods_polygons_completeness <- tetrapods_polygons_key %>%
  filter(!speciesKey %in% lost_occurences) %>% 
  group_by(speciesKey) %>% 
  summarise(geometry = st_union(geometry))

# Montar a lista
st_crs(tetrapods_polygons_completeness) == st_crs(geographic_shape_data)

# Obter a area total de cada unidade administrativa
admin_units_sf <- geographic_shape_data %>%
  mutate(area_total = st_area(geometry))

# Interseccao das especies com as unidades administrativas
species_admin_intersect <- st_intersection(
  tetrapods_polygons_completeness %>% select(speciesKey, geometry),
  admin_units_sf %>% select(name_en, area_total, geometry)
) %>%
  mutate(
    area_species_in_unit = st_area(geometry),
    cover_percent = as.numeric(area_species_in_unit / area_total * 100)) %>%
  st_set_geometry(NULL)  %>%
  # adicionando dados das especies
  left_join(tetrapods_key_species %>%
              select("speciesKey",
                     "scientificName",
                     "Class", "Order","Family"), by = "speciesKey") 

# Aqui estamos gerando algumas listas de ocorrencias esperadas (poligonos) e
# listas de ocorrencias observadas (pontos de ocorrencia). Essas listas
# serao usadas para calcular o completeness usando diferentes 'cortes' para
# testar a robustez dos nossos dados.

# lista de todas ocorrencias observadas (pontos de ocorrencia)
# lista de todas as ocorrencias dos poligonos possiveis (algumas spp tem no darwin core team, mas nao tem ocorrencia)
# lista de especies por und adm 
# lista de especies que tem occ&poligono nest()
# lista de especies com cobertura > 1% na und adm.
# lista de especies com cobertura > 0.5% na und adm
# lista de especies com cobertura > 0.1% na und adm

# Todas as ocorrencias observadas por unidade administrativa
occ_nested <- data_wallacean_unnested %>%
  dplyr::select(speciesKey, name_en) %>%
  distinct() %>%
  group_by(name_en) %>%
  nest() %>%
  rename(observed_df = data) %>%
  mutate(observerd_rich = map_int(observed_df, ~ n_distinct(.x$speciesKey)))

# Todas as especies que tem poligonos mas podem nao ter ocorrencia
species_admin_all <- species_admin_intersect %>%
  group_by(name_en) %>%
  nest() %>%
  rename(all_species = data) %>%
  mutate(expected_rich = map_int(all_species, ~ n_distinct(.x$speciesKey))) %>%
  left_join(occ_nested, by = "name_en") %>%
  mutate(
    joined_df = map2(
      all_species,
      observed_df,
      ~ left_join(
        .x, 
        if (is.null(.y)) tibble(speciesKey = integer(), observed = logical()) else .y %>% mutate(observed = TRUE),
        by = "speciesKey"
      ) %>%
        mutate(observed = if_else(is.na(observed), FALSE, TRUE))
    )
  ) %>%
  mutate(
    observed_rich_all = map_int(joined_df, ~ sum(.x$observed)),
    completeness_all  = observed_rich_all / expected_rich
  )

# Corte de 1%
species_admin_1perc <- species_admin_intersect %>%
  filter(cover_percent >= 1) %>%
  group_by(name_en) %>%
  nest() %>%
  rename(perc_1 = data) %>%
  mutate(expected_rich = map_int(perc_1, ~ n_distinct(.x$speciesKey))) %>%
  left_join(occ_nested, by = "name_en") %>%
  mutate(
    joined_df = map2(
      perc_1,
      observed_df,
      ~ left_join(
        .x, 
        if (is.null(.y)) tibble(speciesKey = integer(), observed = logical()) else .y %>% mutate(observed = TRUE),
        by = "speciesKey"
      ) %>%
        mutate(observed = if_else(is.na(observed), FALSE, TRUE))
    )
  ) %>%
  mutate(
    observed_rich_all = map_int(joined_df, ~ sum(.x$observed)),
    completeness_all  = observed_rich_all / expected_rich
  )

# Corte de 0.5%
species_admin_05perc <- species_admin_intersect %>%
  filter(cover_percent >= 0.5) %>%
  group_by(name_en) %>%
  nest() %>%
  rename(perc_05 = data) %>%
  mutate(expected_rich = map_int(perc_05, ~ n_distinct(.x$speciesKey))) %>%
  left_join(occ_nested, by = "name_en")  %>%
  mutate(
    joined_df = map2(
      perc_05,
      observed_df,
      ~ left_join(
        .x, 
        if (is.null(.y)) tibble(speciesKey = integer(), observed = logical()) else .y %>% mutate(observed = TRUE),
        by = "speciesKey"
      ) %>%
        mutate(observed = if_else(is.na(observed), FALSE, TRUE))
    )
  ) %>%
  mutate(
    observed_rich_all = map_int(joined_df, ~ sum(.x$observed)),
    completeness_all  = observed_rich_all / expected_rich
  )

# Corte de 0.1%
species_admin_01perc <- species_admin_intersect %>%
  filter(cover_percent >= 0.1) %>%
  group_by(name_en) %>%
  nest() %>%
  rename(perc_01 = data) %>%
  mutate(expected_rich = map_int(perc_01, ~ n_distinct(.x$speciesKey))) %>%
  left_join(occ_nested, by = "name_en") %>%
  mutate(
    joined_df = map2(
      perc_01,
      observed_df,
      ~ left_join(
        .x, 
        if (is.null(.y)) tibble(speciesKey = integer(), observed = logical()) else .y %>% mutate(observed = TRUE),
        by = "speciesKey"
      ) %>%
        mutate(observed = if_else(is.na(observed), FALSE, TRUE))
    )
  ) %>%
  mutate(
    observed_rich_all = map_int(joined_df, ~ sum(.x$observed)),
    completeness_all  = observed_rich_all / expected_rich
  ) 

# Apenas especies com ocorrencia e com poligonos
species_admin_occ_polygons <- species_admin_intersect %>%
  filter(speciesKey %in% keys_with_occ) %>%
  group_by(name_en) %>%
  nest() %>%
  rename(occ_polygons = data) %>%
  mutate(expected_rich = map_int(occ_polygons, ~ n_distinct(.x$speciesKey))) %>%
  left_join(occ_nested, by = "name_en") %>%
  mutate(
    joined_df = map2(
      occ_polygons,
      observed_df,
      ~ left_join(
        .x, 
        if (is.null(.y)) tibble(speciesKey = integer(), observed = logical()) else .y %>% mutate(observed = TRUE),
        by = "speciesKey"
      ) %>%
        mutate(observed = if_else(is.na(observed), FALSE, TRUE))
    )
  ) %>%
  mutate(
    observed_rich_all = map_int(joined_df, ~ sum(.x$observed)),
    completeness_all  = observed_rich_all / expected_rich
  ) 

species_admin_occ_polygons %>%
  filter(observed_rich_all != observerd_rich | is.na(observerd_rich))

# Summarise statistics 
specieslist_count <- species_admin_all %>%
  select(name_en, expected_rich, completeness_all) %>%
  rename(expected_rich_all = expected_rich,
         completeness_all_all = completeness_all) %>%
  left_join(
    species_admin_occ_polygons %>%
      select(name_en, expected_rich, completeness_all) %>%
      rename(expected_rich_occpolygon = expected_rich,
             completeness_all_occpolygon = completeness_all), by = "name_en") %>%
  left_join(
    species_admin_1perc %>%
      select(name_en, expected_rich, completeness_all) %>%
      rename(expected_rich_1 = expected_rich,
            completeness_all_1 = completeness_all), by = "name_en") %>% 
  left_join(
    species_admin_05perc %>%
      select(name_en, expected_rich, completeness_all) %>%
      rename(expected_rich_05 = expected_rich,
             completeness_all_05 = completeness_all), by = "name_en"
  ) %>%
  left_join(
    species_admin_01perc %>%
      select(name_en, expected_rich, completeness_all) %>%
      rename(expected_rich_01 = expected_rich,
             completeness_all_01 = completeness_all), by = "name_en"
  ) 

# Completeness species-level
completeness01 <- species_admin_01perc %>%
  select(name_en, perc_01) %>%
  unnest(perc_01) %>%
  select(speciesKey, name_en) %>%
  group_by(speciesKey) %>%
  summarise(ua_total01 = n_distinct(name_en))

completeness05 <- species_admin_05perc %>%
  select(name_en, perc_05) %>%
  unnest(perc_05) %>%
  select(speciesKey, name_en) %>%
  group_by(speciesKey) %>%
  summarise(ua_total05 = n_distinct(name_en))

completeness1 <- species_admin_1perc %>%
  select(name_en, perc_1) %>%
  unnest(perc_1) %>%
  select(speciesKey, name_en) %>%
  group_by(speciesKey) %>%
  summarise(ua_total1 = n_distinct(name_en))

save(
  specieslist_count,         # data.frame com as contagens por threshold
  occ_nested,                # lista de todas ocorrencias 
  species_admin_intersect,   # lista de todas as ocorrencias dos poligonos 
  species_admin_all,         # lista de especies por und adm nest()
  species_admin_occ_polygons,# lista de especies que tem occ&poligono nest()
  keys_with_occ,             # chaves das especies que tem occ baixadas
  species_admin_1perc,       # lista de especies com cobertura > 1% na und adm.
  species_admin_05perc,      # lista de especies com cobertura > 0.5% na und adm
  species_admin_01perc,      # lista de especies com cobertura > 0.1% na und adm
  completeness01,
  completeness05,
  completeness1,
  file = file.path(
    local_directory,
    "00_raw_data",
    "tetrapods_list_sensitivity.RData"
  )
)
