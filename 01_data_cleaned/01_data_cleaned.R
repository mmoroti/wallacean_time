# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    renv::restore()
    #install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

comparative_plot <- function(data_original, data_filtered,
                             year_column, title = "Comparative Plot") {
  
  # Criar tabelas de frequência para os anos
  df_filter <- as.data.frame(table(data_filtered[[year_column]])) %>%
    rename(Year = Var1)
  df_original <- as.data.frame(table(data_original[[year_column]])) %>%
    rename(Year = Var1)
  
  # Converter a coluna 'Year' para numérico, caso seja fator
  df_filter$Year <- as.numeric(as.character(df_filter$Year))
  df_original$Year <- as.numeric(as.character(df_original$Year))
  
  # Fazer o join e preparar o dataframe comparativo
  df_comparative <- full_join(df_original, df_filter, by = "Year") %>%
    rename(Original = Freq.x, Filtered = Freq.y) %>%
    pivot_longer(cols = c(Original, Filtered), names_to = "Dataset", values_to = "Occurences") %>%
    replace_na(list(Occurences = 0))
  
  # Criar o gráfico comparativo
  ggplot(df_comparative, aes(x = Year, y = Occurences, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 2, alpha = 0.6) +
    geom_line(stat = "smooth", aes(color = Dataset), linewidth = 1.0, linetype = "dashed") + 
    theme_bw() +
    labs(title = title, y = "Frequency of occurrences") +
    scale_x_continuous(breaks = c(1900, 1940, 1980, 2020), limits = c(1900, 2020)) +
    theme(
    axis.title.x = element_text(size = 18),      # Tamanho do título do eixo x
    axis.title.y = element_text(size = 18),      # Tamanho do título do eixo y
    axis.text.x = element_text(size = 16),       # Tamanho dos textos do eixo x
    axis.text.y = element_text(size = 16),       # Tamanho dos textos do eixo y
    plot.title = element_text(size = 18, hjust = 0.5)  # Tamanho do título do gráfico e centralização
  )
}

validar_dados <- function(data) {
  # Lista para armazenar mensagens de validação
  inconsistencias <- list()
  # Verificação de latitude = 0
  if (nrow(data %>% filter(decimalLatitude == 0 & decimalLongitude != 0)) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Latitude igual a 0 com Longitude diferente de 0")
  }
  # Verificação de longitude = 0
  if (nrow(data %>% filter(decimalLatitude != 0 & decimalLongitude == 0)) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Longitude igual a 0 com Latitude diferente de 0")
  }
  # Verificação de latitude e longitude = 0
  if (nrow(data %>% filter(decimalLatitude == 0 & decimalLongitude == 0)) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Latitude e Longitude iguais a 0")
  }
  # Verificação de anos com NA
  if (nrow(data %>% filter(is.na(year))) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Ano com valor NA")
  }
  # Verificação de coordenadas com com NA
  if (nrow(data %>% filter(is.na(decimalLatitude) | is.na(decimalLongitude))) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Ano com valor NA")
  }
  # Verificação de anos igual a 0
  if (nrow(data %>% filter(year == 0)) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Ano igual a 0")
  }
  # Verificação de anos no futuro
  if (nrow(data %>% filter(year > 2024)) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Ano posterior a 2024")
  }
  # Verificação de anos anteriores a 1900
  if (nrow(data %>% filter(year < 1900)) > 0) {
    inconsistencias <- append(inconsistencias, "Inconsistência: Ano anterior a 1900")
  }
  # Retornar resultado
  if (length(inconsistencias) == 0) {
    print("Nenhuma inconsistência encontrada")
  } else {
    print("Inconsistências encontradas:")
    print(inconsistencias)
  }
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
local_directory <- file.path("E:",
                             "datasets_centrais",
                             "wallacean_time") 

# FILTER AND CLEAN ----
## GBIF ----
# load data
rm(list = setdiff(ls(), c("local_directory",
                          "validar_dados",
                          "comparative_plot"))); gc()
load(file.path(
  local_directory,
  "00_raw_data",
  "tetrapodstraits_data.RData")
)

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID) 
nrow(data_tetrapodstraits) 
# 24.591.154 Global species occurrences
table(data_tetrapodstraits$class)

# Amphibia 4745415                         
# Aves 7533355                                     
# Crocodylia 17970
# Mammalia 8218550
# Sphenodontia 420
# Squamata 3810639
# Testudines 26480

data_tetrapods_clean <- data_tetrapodstraits %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "gbif") %>%
  mutate(eventDate = na_if(eventDate, "")) %>%
  # U. ricardinii aparece sem classificacao de classe
  mutate(class = if_else(class == "", "Squamata", class))

nrow(data_tetrapods_clean) # 24.402.844

# how much data is na?
data_tetrapods_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year)),
            missing_date = sum(is.na(eventDate)))
# 9.001.840 sem coordenada
# 3.626.524 sem ano 
# 6.563.455 sem mes
# 7.408.408 sem dia

# exists a difference between "year" and "event date" 
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
#View(data_tetrapods_clean %>%
#       filter(!is.na(eventDate) & is.na(year)))
data_tetrapods_filter <- data_tetrapods_clean %>%
  drop_na(decimalLatitude,
          decimalLongitude,
          year) %>%
  select(-eventDate)

data_tetrapods_filter %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)), # without missing data =)
            missing_long = sum(is.na(decimalLongitude)), 
            missing_day = sum(is.na(day)), # 2.776.369 NA's
            missing_month = sum(is.na(month)), # 2.257.939 NA's
            missing_year = sum(is.na(year))) # without missing data =)
#missing_date = sum(is.na(eventDate))) 

nrow(data_tetrapods_filter) # 14.141.511 Global species occurrences
table(data_tetrapods_filter$class) # occurrences per classes 

# coordsclean
# convert country code from ISO2c to ISO3c
data_tetrapods_filter$countryCode <-  countrycode(
  data_tetrapods_filter$countryCode,
  origin =  'iso2c',
  destination = 'iso3c',
  custom_match = c("XK" = "XKX"))

flags <- clean_coordinates(x = data_tetrapods_filter,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))
# Flagged 678.415 of 14141511 records, EQ = 0.05.

#Exclude problematic records
data_tetrapods_filter_spatialpoints <- data_tetrapods_filter[
  flags$.summary,
] 

nrow(data_tetrapods_filter_spatialpoints) # 13.463.096 occ

# filter coordinates uncertain > 100km
data_tetrapods_filtered <- data_tetrapods_filter_spatialpoints %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
           is.na(coordinateUncertaintyInMeters)) 
nrow(data_tetrapods_filtered) # 13.053.975 

## Human observation non-birds ----
rm(list = setdiff(ls(), c("local_directory",
                          "validar_dados",
                          "comparative_plot",
                          "data_tetrapods_filtered"))); gc()

load(file.path(
  local_directory,
  "00_raw_data",
  "nonbirds_humanobservation_data.RData")
)

# key of each occurence
anyDuplicated(data_nonbirds$gbifID) 
nrow(data_nonbirds) # 38.702.092 Global species occurrences

table(data_nonbirds$class)
# Amphibia 6.484.723                                     
# Crocodylia 92.317
# Mammalia 25.644.889
# Sphenodontia 1.376
# Squamata 5.634.057
# Testudines 844.730

data_nonbirds_clean <- data_nonbirds %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "gbif") %>%
  mutate(eventDate = na_if(eventDate, "")) 

nrow(data_nonbirds_clean) # 38.317.234 Global species occurrences
table(data_nonbirds_clean$class)

# how much data is na?
data_nonbirds_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year)),
            missing_date = sum(is.na(eventDate)))
# 344.879 sem coordenada
# 1.415.954 sem ano 
# 3.887.627 sem mes
# 4.598.425 sem dia

# exists a difference between "year" and "event date" of 7280 occurences
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
#View(data_tetrapods_clean %>%
#       filter(!is.na(eventDate) & is.na(year)))
data_nonbirds_filter <- data_nonbirds_clean %>%
  drop_na(decimalLatitude,
          decimalLongitude,
          year) #%>%
  #select(-eventDate)

data_nonbirds_filter %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)), # without missing data =)
            missing_long = sum(is.na(decimalLongitude)), 
            missing_day = sum(is.na(day)), # 3.080.541 NA's
            missing_month = sum(is.na(month)), # 2.374.207 NA's
            missing_year = sum(is.na(year))) # without missing data =)

nrow(data_nonbirds_filter) # 36.590.997 Global species occurrences
table(data_nonbirds_filter$class) # occurences per classes 

# Coordsclean
# convert country code from ISO2c to ISO3c
data_nonbirds_filter$countryCode <-  countrycode(
  data_nonbirds_filter$countryCode,
  origin =  'iso2c',
  destination = 'iso3c',
  custom_match = c("XK" = "XKX"))

flags <- clean_coordinates(x = data_nonbirds_filter,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))
# Flagged 1381710 of 36590997 records, EQ = 0.04.

#Exclude problematic records
data_nonbirds_filter_spatialpoints <- data_nonbirds_filter[
  flags$.summary,
]
nrow(data_nonbirds_filter_spatialpoints) # 35.209.287 Global species occ

# filter coordinates uncertain > 100km
# filter years from 1900
data_nonbirds_filtered <- data_nonbirds_filter_spatialpoints %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
           is.na(coordinateUncertaintyInMeters)) 

## Human observation birds ----
# Extract ZIP file from eBird
#arquivo_zip <- file.path(local_directory,
#                         "00_raw_data",
#                         "0018107-251025141854904.zip")
#temp_dir <- file.path(local_directory,
#                      "00_raw_data",
#                      "temp_extract")
#dir.create(temp_dir, showWarnings = FALSE)
#unzip(arquivo_zip, exdir = temp_dir)

# save raw .parquet files 
parquet_dir <- file.path(local_directory, 
                         "00_raw_data",
                         "parquet_por_ordem")

# save clean .parquet files
parquet_clean <- file.path(local_directory,
                           "01_data_cleaned",
                           "parquet_clean")

# Obtain human observation data
nested_cols <- c("speciesKey", "class", "order", 
                 "family", "species", "gbifID",
                 "day", "month", "year",
                 "decimalLatitude", "decimalLongitude",
                 "coordinateUncertaintyInMeters", "occurrenceStatus",
                 "taxonRank", "countryCode", "eventDate")

ebird_data <- arrow::open_delim_dataset(
  file.path(temp_dir, "0018107-251025141854904.csv"), 
  parse_options = csv_parse_options(
    delimiter = "\t", 
    newlines_in_values = TRUE),
  read_options = list(
    block_size = 1000000000L, 
    use_threads = TRUE)) %>%
  select(all_of(nested_cols))

# create multiple .parquet files
ebird_data |> 
  group_by(order) |> 
  write_dataset(path = parquet_dir,
                format = "parquet")

tibble(
  arquivos = list.files(parquet_dir, recursive = TRUE),
  tamanho_MB = file.size(file.path(parquet_dir, arquivos)) / 1024^2
) %>% View()

# Cleaning each Order
ordens_disponiveis <- list.dirs(parquet_dir,
                                recursive = FALSE,
                                full.names = FALSE)

ordens_disponiveis <- ordens_disponiveis[
  ordens_disponiveis != "order=Passeriformes"]

# Estatísticas finais
estatisticas <- data.frame(
  order = character(),
  registros_originais = integer(),
  registros_clean = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(ordens_disponiveis)) {
  
  tryCatch({
    cat("Processando", ordens_disponiveis[i])
    
    aves <- open_dataset(file.path(parquet_dir, ordens_disponiveis[i]))
    aves_clean <- aves %>%
      collect()
    
    n_original <- nrow(aves_clean)
    
    # Filter start
    aves_clean <- aves_clean %>%
      filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
      filter(occurrenceStatus != "ABSENT") %>%
      mutate(origin_of_data = "gbif") %>%
      mutate(eventDate = na_if(eventDate, "")) %>%
      drop_na(decimalLatitude,
              decimalLongitude,
              year)
    
    # Remove invalidity coordinates
    aves_clean <- aves_clean %>%
      filter(
        !is.na(decimalLatitude),
        !is.na(decimalLongitude),
        decimalLatitude >= -90 & decimalLatitude <= 90,
        decimalLongitude >= -180 & decimalLongitude <= 180
      ) 
    
    aves_clean$countryCode <-  countrycode(
      aves_clean$countryCode,
      origin =  'iso2c',
      destination = 'iso3c',
      custom_match = c("XK" = "XKX"))
    
    # CoordinateCleaner 
    flags <- clean_coordinates(
      x = aves_clean,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      countries = "countryCode",
      species = "species",
      tests = c("equal", "gbif", "institutions", "seas", "zeros")
    )
    
    aves_clean <- aves_clean[flags$.summary, ] %>%
      filter(
        coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters),
        year < 2026 
      )
    
    aves_clean <- aves_clean %>%
      mutate(
        # duas casas decimais = 1.11 km 
        decimalLatitude = round(decimalLatitude, 2), 
        decimalLongitude = round(decimalLongitude, 2),
        month = ifelse(month > 12 | month < 1, NA, month),
        day = ifelse(day > 31 | day < 1, NA, day)          
      ) %>%
      arrange(speciesKey, year, month, day) 
    
    duplicated_flags <- cc_dupl(
      aves_clean,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "speciesKey",
      value = "flagged"
    )
    
    aves_clean <- aves_clean[
      duplicated_flags == TRUE, ]
    
    # Filter end
    n_final <- nrow(aves_clean)
    
    # Save new parquet file - Ordem
    ordem_nome <- sub("/.*", "", sub("^order=", "", ordens_disponiveis[i]))
    write_parquet(aves_clean, 
                  sink = file.path(parquet_clean,
                                   paste0("clean_", ordem_nome, ".parquet")))
    # How much data did we lose?
    estatisticas <- rbind(estatisticas, data.frame(
      taxa = ordem_nome,
      registros_originais = n_original,
      registros_clean = n_final
    ))
    
  }, error = function(e) {
    cat("\n❌ Erro em", ordem_nome, ":", e$message, "\n")
  })
}

# How much data was cleaned? 
print(estatisticas)
write.table(estatisticas, file.path(
  local_directory, "00_raw_data", "clean_statistics",
  "clean_statistics_order.txt"))

# Passeriformes cannot be processed; let's separate them into families.
aves_passeriformes <- open_dataset(file.path(parquet_dir, 	
                                             "order=Passeriformes/part-0.parquet"))

parquet_dir_family <- file.path(local_directory, 
                                "parquet_por_ordem",
                                "order=Passeriformes",
                                "family")

aves_passeriformes |>
  group_by(family) |> 
  write_dataset(path = parquet_dir_family,
                format = "parquet")

tibble(
  arquivos = list.files(parquet_dir_family, recursive = TRUE),
  tamanho_MB = file.size(file.path(parquet_dir_family, arquivos)) / 1024^2
) %>% View() 

familias_disponiveis <- list.dirs(parquet_dir_family,
                                  recursive = FALSE,
                                  full.names = FALSE) 

# Estatisticas finais
estatisticas <- data.frame(
  ordem = character(),
  registros_originais = integer(),
  registros_clean = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(familias_disponiveis)) {
  
  tryCatch({
    cat("Processando", familias_disponiveis[i])
    
    aves <- open_dataset(file.path(parquet_dir_family, familias_disponiveis[i]))
    aves_clean <- aves %>%
      collect()
    
    n_original <- nrow(aves_clean)
    
    # Filter start
    aves_clean <- aves_clean %>%
      filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
      filter(occurrenceStatus != "ABSENT") %>%
      mutate(origin_of_data = "gbif") %>%
      mutate(eventDate = na_if(eventDate, "")) %>%
      drop_na(decimalLatitude,
              decimalLongitude,
              year)
    
    # Remove invalidity coordinates
    aves_clean <- aves_clean %>%
      filter(
        !is.na(decimalLatitude),
        !is.na(decimalLongitude),
        decimalLatitude >= -90 & decimalLatitude <= 90,
        decimalLongitude >= -180 & decimalLongitude <= 180
      ) 
    
    aves_clean$countryCode <-  countrycode(
      aves_clean$countryCode,
      origin =  'iso2c',
      destination = 'iso3c',
      custom_match = c("XK" = "XKX"))
    
    # CoordinateCleaner 
    flags <- clean_coordinates(
      x = aves_clean,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      countries = "countryCode",
      species = "species",
      tests = c("equal", "gbif", "institutions", "seas", "zeros")
    )
    
    aves_clean <- aves_clean[flags$.summary, ] %>%
      filter(
        coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters),
        year < 2026 #year > 1899 & 
      )
    
    aves_clean <- aves_clean %>%
      mutate(
        # duas casas decimais = 1.11 km 
        decimalLatitude = round(decimalLatitude, 2), 
        decimalLongitude = round(decimalLongitude, 2),
        month = ifelse(month > 12 | month < 1, NA, month),
        day = ifelse(day > 31 | day < 1, NA, day)          
      ) %>%
      arrange(speciesKey, year, month, day) 
    
    duplicated_flags <- cc_dupl(
      aves_clean,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "speciesKey",
      value = "flagged"
    )
    
    aves_clean <- aves_clean[
      duplicated_flags == TRUE, ]
    
    # Filter end
    n_final <- nrow(aves_clean)
    
    # Save new parquet file - Family
    family_nome <- sub("/.*", "", sub("^family=", "", familias_disponiveis[i]))
    write_parquet(aves_clean, 
                  sink = file.path(parquet_clean,
                                   paste0("clean_", family_nome, ".parquet")))
    
    # How much data did we lose?
    estatisticas <- rbind(estatisticas, data.frame(
      taxa = family_nome,
      registros_originais = n_original,
      registros_clean = n_final
    ))
    
  }, error = function(e) {
    cat("\n❌ Erro em",  family_name, ":", e$message, "\n")
  })
}

# How much data was cleaned? 
print(estatisticas)
write.table(estatisticas,file.path(
  local_directory, "00_raw_data", "clean_statistics",
  "clean_statistics_family.txt"))

## BioTIME 2.0v ----
rm(list = setdiff(ls(), c("local_directory",
                          "validar_dados",
                          "comparative_plot",
                          "data_tetrapods_filtered",
                          "data_nonbirds_filtered"))); gc()
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

data_biotime_clean %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_day = sum(is.na(day)),     # 863.686
            missing_month = sum(is.na(month)), # 766.064
            missing_year = sum(is.na(year)))   # without missing data

flags <- clean_coordinates(x = data_biotime_clean,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))
# Flagged 1068661 of 3506963 records, EQ = 0.3.

#Exclude problematic records
data_biotime_filtered <- data_biotime_clean[
  flags$.summary,
]

## speciesLink ----
rm(list = setdiff(ls(), c("local_directory",
                          "validar_dados",
                          "comparative_plot",
                          "data_tetrapods_filtered",
                          "data_nonbirds_filtered",
                          "data_biotime_filtered"))); gc()
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
    day = as.numeric(na_if(day, "")),
    month = as.numeric(na_if(month, "")),
    year = as.numeric(na_if(year, ""))) %>%
  mutate(
    month = ifelse(month > 12 | month < 1, NA, month),
    day = ifelse(day > 31 | day < 1, NA, day)          
  ) %>%
  mutate(coordinateprecision = na_if(coordinateprecision, "")) %>% 
  mutate(year = as.integer(as.character(year))) %>%
  filter(!is.na(year))

data_splink_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),     # 438.648
            missing_month = sum(is.na(month)), # 473.340
            missing_year = sum(is.na(year)))

# invalidity coordinates
data_splink_cleaned <- data_splink_clean %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    decimalLatitude >= -90 & decimalLatitude <= 90,
    decimalLongitude >= -180 & decimalLongitude <= 180
  ) 

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

# filter coordinates uncertain > 100km
# filter years from 1900
data_splink_filtered <- data_splink_filter_spatialpoints %>%
  mutate(coordinateprecision = as.numeric(coordinateprecision)) %>%
  filter(coordinateprecision / 1000 <= 100 | is.na(coordinateprecision)) %>%
  filter(year < 2026)

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

# CHECK DATASETS ----
dados_perdidos <- read.csv2("Figures/cleaned_data.csv") # atualizar manualmente
totais_bases <- c(GBIF = nrow(data_tetrapods_clean),
                  BioTime = nrow(data_biotime_clean),
                  Specieslink = nrow(data_splink_clean))

dados_long <- dados_perdidos %>%
  pivot_longer(cols = c(GBIF, BioTime, Specieslink),
               names_to = "font", 
               values_to = "value") #%>%
               #select(-Total) %>%
               #filter(Filter != "Total temporal filter" & Filter != "Total spatial filter" )

# Calcular a porcentagem dentro de cada "Filter"
dados_long <- dados_long %>%
  mutate(total_specie = case_when(
    font == "GBIF" ~ totais_bases["GBIF"],
    font == "BioTime" ~ totais_bases["BioTime"],
    font == "Specieslink" ~ totais_bases["Specieslink"]
  )) %>%
  mutate(percentage = (value / total_specie) * 100)

plot <- ggplot(dados_long, aes(x = font, y = percentage, fill = Filter)) +
  geom_bar(stat = "identity", position = "fill") +  # Barras empilhadas normalizadas para 100%
  #facet_wrap(~ specie) +  # Um gráfico para cada base
  labs(y = "Proportion of occurrences removed (%)", x = "", fill = "Filtro") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Eixo Y em porcentagem
  coord_flip() +  # Inverte os eixos
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Aumenta o tamanho do texto
    axis.title.x = element_text(size = 16),  # Título do eixo X
    axis.title.y = element_text(size = 16),  # Título do eixo Y
    legend.title = element_blank(),  # Título da legenda
    legend.text = element_text(size = 12),  # Texto da legenda
    legend.justification = c(0.9, 1),  # Justifica a legenda em relação ao seu centro
    legend.position = "bottom",  # Posiciona a legenda em cima
    legend.direction = "horizontal"  # Direção horizontal
  ) +
  guides(fill = guide_legend(ncol = 5))

plot
ggsave("Figures/filter_per_data.png", plot, width = 12, height = 6, dpi = 300)

#--- GBIF
nrow(data_tetrapods_clean) - nrow(data_tetrapods_filtered)
# check difference per species
length(unique(data_tetrapods_clean$speciesKey)) # 9786 spp
length(unique(data_tetrapods_filtered$speciesKey)) # 7897 spp

# check difference per class
count_species_clean <- data_tetrapods_clean %>%
  filter(year >1899) %>%
  distinct(speciesKey, .keep_all = TRUE)
table(count_species_clean$class)

count_species_filter <- data_tetrapods_filtered %>%
  distinct(speciesKey, .keep_all = TRUE)
table(count_species_filter$class)

# Comparative plots of original vs filtered GBIF datasets
gbif_plot <- comparative_plot(
  data_tetrapods_clean %>% filter(year >1899),
  data_tetrapods_filtered,
  year_column = "year",
  title = "GBIF dataset")

# coordenadas ou anos com zero e/ou NA
# algum ano < 1900 ou > 2024
validar_dados(data_tetrapods_filtered)

#data_tetrapods_filtered %>% 
#  summarise(missing_lat = sum(is.na(decimalLatitude)),
#            missing_long = sum(is.na(decimalLongitude)), 
#            missing_year = sum(is.na(year))) # without missing data

#--- BIOTIME
# quantos registros foram perdidos depois do coordscleaned
nrow(data_biotime_clean) - nrow(data_biotime_filtered)

# check difference per class
count_species_clean <- data_biotime_clean %>%
  #filter(year > 1899) %>%
  distinct(speciesKey, .keep_all = TRUE)
table(count_species_clean$class)

count_species_filter <- data_biotime_filtered %>%
  distinct(speciesKey, .keep_all = TRUE)
table(count_species_filter$class)

biotime_plot <- comparative_plot(data_biotime_clean,
                 data_biotime_filtered,
                 "year", title = "BioTime dataset")

# coordenadas ou anos com zero e/ou NA
# algum ano < 1900 ou > 2024
validar_dados(data_biotime_filtered)

#--- SPECIESLINK
nrow(data_splink_clean) - nrow(data_splink_filtered) # 303.476
# the wrongs years difficult to visualize temporal series tendence
data_splink_clean_compare <- data_splink_clean %>%
  filter(year >1899 & year < 2025)

# check difference per class
count_species_clean <- data_splink_clean %>%
  distinct(speciesKey, .keep_all = TRUE)
nrow(count_species_clean) #7330

count_species_filter <- data_splink_filtered %>%
  distinct(speciesKey, .keep_all = TRUE)
nrow(count_species_filter) #6171

splink_plot <- comparative_plot(
  data_splink_clean_compare,
  data_splink_filtered,
  "year", 
  title = "splink dataset")

# coordenadas ou anos com zero e/ou NA
# algum ano < 1900 ou > 2024
validar_dados(data_splink_filtered)

# plot grid
cowplot::plot_grid(gbif_plot, biotime_plot, splink_plot, align = "v", ncol=1)
ggsave("Figures/comparing_filter_data.tiff", dpi = 300, units = "in")

# UNIFYING OCCURENCES DATA ---- 
# load data
rm(list = setdiff(ls(), c("local_directory"))); gc()

# Eu tenho que remover todas as chaves que NAO estao em 
# `species_list_tetrapods_filter` , porque aqui sao todas as especies 
# encontradas no TetrapodTraits, e podem ter saido no data_tetrapods_filtered
# por exemplo, Cristidorsa otai esta no TetrapodTraits, mas todas as ocorrencias
# estao sem data nas observaçoes feitas por especialistas, mas tem dados nas occ
# de observacoes humanas com data. Ou seja, essa chave deve ser mantida. 
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
  sources = parquet_dir
) %>%
  collect()

load(file.path(
  local_directory,
  "01_data_cleaned",
  "datasets_filtered.RData")
)

# BioTIME chaves perdidas
chaves_perdidas_biotime <- setdiff(
  data_biotime_filtered$speciesKey,
  species_list_tetrapods_filter$speciesKey
) 
# Mantenha APENAS as speciesKeys que NÃO estão na lista de perdidas
data_biotime_filtered <- data_biotime_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month))

# speciesLink chaves perdidas
chaves_perdidas_splink <- setdiff(
  data_splink_filtered$speciesKey,
  species_list_tetrapods_filter$speciesKey
) 

data_splink_filtered <- data_splink_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month))

# Non-Birds HumanObsevation chaves perdidas
chaves_perdidas_nonbirdsHO <- setdiff(
  data_nonbirds_filtered$speciesKey,
  species_list_tetrapods_filter$speciesKey
) # nenhuma chave eh perdida, pq todas foram baixadas baseado no speciesKey 

data_nonbirds_filtered <- data_nonbirds_filtered %>% 
  #filter(!(speciesKey %in% chaves_perdidas_nonbirdsHO)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month)) %>%
  mutate(origin_of_data = "HUMAN_OBSERVATION")

# Birds HumanObsevation chaves perdidas
chaves_perdidas_birdsHO <- setdiff(
  data_birds$speciesKey,
  species_list_tetrapods_filter$speciesKey
) # nenhuma chave eh perdida, pq todas foram baixadas baseado no speciesKey

data_birds <- data_birds %>% 
  #filter(!(speciesKey %in% chaves_perdidas_birdsHO)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month)) %>%
  mutate(origin_of_data = "HUMAN_OBSERVATION")

data_occurences_precleaned <- bind_rows(
  data_tetrapods_filtered,
  data_biotime_filtered,
  data_splink_filtered,
  data_nonbirds_filtered,
  data_birds)

nested_cols <- c("speciesKey", "class", "order", "family", "species",
                 "gbifID", "day", "month", "year",
                 "decimalLatitude", "decimalLongitude", "origin_of_data")

data_occurences_precleaned <- data_occurences_precleaned %>% 
  select(all_of(nested_cols))

write_parquet(data_occurences_precleaned,  # dados limpos e integrados
              sink = file.path(
                local_directory,
                "01_data_cleaned",
                "data_occurences_filtered.parquet"))

# Criar dataframe combinando todos os vetores
# TODO da para usar para baixar o registro de occ para mais especies com essas 
# chaves, e buscar + poligonos se necessarios. 
chaves_perdidas_combinadas <- bind_rows(
  data.frame(speciesKey = chaves_perdidas_biotime, Database = "Biotime"),
  data.frame(speciesKey = chaves_perdidas_splink, Database = "splink"),
  #data.frame(speciesKey = chaves_perdidas_nonbirdsHO, Database = "nonbirdsHO"),
  #data.frame(speciesKey = chaves_perdidas_birdsHO, Database = "birdsHO")
) %>%
  arrange(speciesKey, Database)

# 858 chaves ao todo perdidas, estao no BioTime e/ou speciesLink, mas nao
# no TetrapodTraits. Essas podemos conferir para ver se entram ou nao em nossos
# dados.
chaves_perdidas_combinadas %>% 
  distinct(speciesKey) %>% nrow() 

write.table(chaves_perdidas_combinadas, file.path(
  local_directory, "chaves_perdidas.txt"))

# REMOVE DUPLICATES ----
rm(list=ls()); gc() # clean local enviroment

# occurences data
data_occurences_precleaned <- open_dataset(
  sources = file.path(
    local_directory,
    "01_data_cleaned",
    "data_occurences_filtered.parquet")
)

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

nrow(data_occurences_precleaned_duplicate)*100/nrow(data_occurences_precleaned)
# 78% dos dados mantidos
length(unique(data_occurences_precleaned_duplicate$speciesKey)) 
# 28.459 global scale spp

# FILTER POINTS BY RANGE POLYGONS ----
# Shapefile base
load(file = file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_polygons_key.RData")) # From TetrapodTraits

# Sensitivity shapefiles expanded 110km and 220km
load(file = file.path(
       local_directory,
       "00_raw_data",
       "tetrapods_polygons_sensitivity.RData"
     )
)

# deixar apenas as especies que tem poligonos
data_wallacean_knownledge_sa <- data_occurences_precleaned_duplicate %>%
  filter(speciesKey %in% tetrapods_polygons_key$speciesKey)
# deixar apenas os poligonos que tem dados de distribuicao
data_tetrapods_sa <- tetrapods_polygons_key %>% 
  filter(speciesKey %in% data_occurences_precleaned_duplicate$speciesKey) 

# usando os poligonos do tetrapodtraits chegamos a 7893 spp
# usando dados de poligonos de outras bases 6527 spp apenas
# tem que ter o mesmo numero de especies (speciesKey)
length(unique(data_wallacean_knownledge_sa$speciesKey)) 
length(unique(data_tetrapods_sa$speciesKey)) 

# run cc_iucn()
# terra::vect(data_tetrapods_sa) se funcionar vai dar certo
data_tetrapods_sa_transf <- st_cast(data_tetrapods_sa,
 "MULTIPOLYGON")
# Verificar se a conversão funcionou
geometry_types <- st_geometry_type(data_tetrapods_sa_transf)
table(geometry_types)

# check if list_occurences_clean is empty before run
unique_species_keys <- unique(data_wallacean_knownledge_sa$speciesKey)
list_occurences_clean <- data.frame()
species_keys_with_errors <- c()

for (species_key in unique_species_keys) {
  # Usar tryCatch para capturar erros
  tryCatch({
    # Selecionar o polígono correspondente à espécie
    especie_poligono <- data_tetrapods_sa_transf %>%
      filter(speciesKey == species_key)
    
    subset_data <- data_wallacean_knownledge_sa %>%
      filter(speciesKey == species_key)
    
    # Aplicar a função cc_iucn
    range_flags <- cc_iucn(
      x = subset_data,
      range = especie_poligono,
      species = "speciesKey",
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      value = "flagged"
    )
    
    # Filtrar os dados usando os valores retornados por range_flags
    subset_data_filtered <- subset_data[range_flags == TRUE, ]
    
    # Adicionar os dados filtrados à lista final
    list_occurences_clean <- rbind(list_occurences_clean, subset_data_filtered)
    
  }, error = function(e) {
    # Se ocorrer um erro, salvar o species_key no vetor
    species_keys_with_errors <<- c(species_keys_with_errors, species_key)
    # Opcional: Printar o erro ou a species_key para depuração
    message(paste("Erro no speciesKey:", species_key, "-", e$message))
  })
}
# Verificar speciesKey que tiveram erros
print("SpeciesKey com erro:")
print(species_keys_with_errors)
#end_time <- Sys.time()
#print(end_time - start_time) # Time difference of 1h56min

#data_tetrapods_sa_transf <- data_tetrapods_sa %>%
#    filter(speciesKey %in% species_keys_with_errors)

# 405.215 ocorrencias usando dados espaciais do  (7893 spp)
# 269.333 ocorrencias usando poligonos de especialistas (5748 spp)
# 263.628 ocorrencias com poligonos e filtro temporal
nrow(list_occurences_clean) # 124.684.269 occ
nrow(data_wallacean_knownledge_sa)-nrow(list_occurences_clean) # 15.018.549 removed
length(unique(list_occurences_clean$speciesKey)) # 27405 spp 
length(
  unique(
    data_wallacean_knownledge_sa$speciesKey))-length(
      unique(list_occurences_clean$speciesKey)) # 1054 spp removed

# Quantas espécies para cada grupo?
list_species <- list_occurences_clean %>%
 distinct(speciesKey, .keep_all=TRUE)
table(list_species$class)
# Amphibia  4229                         
# Aves      6617                                                   
# Crocodylia 17
# Mammalia 4139
# Sphenodontia 1
# Squamata 6787
# Testudines 257

write_parquet(list_occurences_clean,  # dados limpos e integrados
              sink = file.path(
                local_directory, 
                "01_data_cleaned",
                "list_occurences_clean.parquet"))

# Compare
plot_sujo <- ggplot() +
  coord_fixed() +
  borders("world", colour = "gray50", fill = "gray50") +
  geom_point(data = data_wallacean_knownledge_sa,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

plot_clean <- ggplot() +
  coord_fixed() +
  borders("world", colour = "gray50", fill = "gray50") +
  geom_point(data = list_occurences_clean,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

# EXTRACTING GEOLOCATIONS ----
parquet_dir <- file.path(local_directory, 
                         "01_data_cleaned",
                         "list_occurences_clean.parquet")
list_occurences_clean <- open_dataset(
  sources = parquet_dir
) %>%
  collect()

load(file = file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_polygons_key.RData")) # From TetrapodTraits

# deixar apenas os poligonos que tem dados de distribuicao
data_tetrapods_sa <- tetrapods_polygons_key %>% 
  filter(speciesKey %in% list_occurences_clean$speciesKey) 

# run cc_iucn()
# terra::vect(data_tetrapods_sa) se funcionar vai dar certo
data_tetrapods_sa_transf <- st_cast(data_tetrapods_sa,
                                    "MULTIPOLYGON")
# Verificar se a conversão funcionou
geometry_types <- st_geometry_type(data_tetrapods_sa_transf)
table(geometry_types)

rm(list = setdiff(ls(), c("local_directory",
                          "list_occurences_clean",
                          "data_tetrapods_sa"
))); gc()

# shapefile to crop adm unit
load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData"))

# Verificar validade das geometrias
sf_use_s2(FALSE)
st_is_valid(data_tetrapods_sa) %>% table()
st_is_valid(geographic_shape_data) %>% table()

# extracting geolocation per specie
data_occurences <- st_as_sf(list_occurences_clean,
                            coords = c("decimalLongitude", "decimalLatitude"),
                            crs = 4326) %>%    # WGS84
  st_transform(st_crs(geographic_shape_data))  # equal area

nrow(data_occurences) == nrow(list_occurences_clean) 

data_occurences <- data_occurences %>%
  mutate(ID = 1:nrow(data_occurences))
any(duplicated(data_occurences$ID))

# validaty of polygons 
any(!st_is_valid(geographic_shape_data))
any(!st_is_valid(data_occurences))

geographic_shape_data_adj <- geographic_shape_data %>%
  select(name_en, geometry)

data_occurences_geo <- st_join(data_occurences,
  geographic_shape_data_adj, 
  join = st_intersects,
  left = TRUE) %>% 
  filter(!is.na(name_en)) %>% # 'col' filtra NA's (buracos) nos poligonos
  group_by(ID) %>%
  filter(n() == 1) %>% # remove registros em polígonos sobrepostos
  ungroup()

any(duplicated(data_occurences_geo$ID))
# 907.533 occ registros fora do poligonos (oceano, areas em litigio, etc)
nrow(data_occurences)-nrow(data_occurences_geo) # 907.533 occ removidas 
names(data_occurences_geo)

data_occurences_geometry <- data_occurences_geo %>%
  select(speciesKey, species, day, month, year, origin_of_data, name_en, ID) 

data_occurences_units <- data_occurences_geometry %>%
  st_drop_geometry()

# SPECIES LIST PER ADMINASTRIVE UNIT ----
interseccao_sf <- st_intersection(data_tetrapods_sa, geographic_shape_data)

list_per_admunit <- interseccao_sf %>%
  st_drop_geometry()

list_per_admunit$name_en %>% table()

save(list_occurences_clean, # sem duplicatas e filtrados por poligonos
     data_occurences_units, # com as unidades administrativas
     list_per_admunit,      # lista de especies por und adm
     file = file.path(
       local_directory,
       "01_data_cleaned",
       "data_occurences_geolocation.RData")
)

# NESTED DATAFRAME WITH BIOLOGICAL TRAITS ----
rm(list = setdiff(ls(), c("local_directory",
                          "list_per_admunit"
))); gc()

load(file = file.path(
  local_directory,
  "01_data_cleaned",
  "data_occurences_geolocation.RData")
)

load(file.path(
  local_directory,
  "00_raw_data",
  "trait_data.RData")
) # match taxonomic and trait information with speciesKey

tetrapods_key_species <- trait_data %>%
  select("speciesKey","scientificName","Class", "Order","Family", 
         "YearOfDescription","BodyLength_mm","ImputedLength","BodyMass_g",
         "ImputedMass","Diu","Noc","Nocturnality", "Fos", "Ter", "Aqu", "Arb", "Aer",
         "ImputedHabitat","Verticality","MajorHabitatSum","ImputedMajorHabitat","RangeSize",
         "HumanDensity", "AssessedStatus", "Latitude", "Elevation", "ETA50K", "AnnuMeanTemp") %>%
  distinct(speciesKey, .keep_all = TRUE) # 117 speciesKey duplicated
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

# para o mundo sem human observation
#Amphibia  Aves   Mammalia Reptilia 
# 4942     8250     4459     7050 
# para o mundo com human observation
#Amphibia  Aves   Mammalia Reptilia 
# 5553     8941     4794     8107 

anyDuplicated(data_wallacean_nested$speciesKey)

# unnest data
data_wallacean_unnested <- data_wallacean_nested %>%
  unnest(cols = c(event_table))
nrow(data_wallacean_unnested) 

# 3.144.303 occ global
# 123.776.736 occ global com Human Observation
save(data_wallacean_nested,
     data_wallacean_unnested,
     list_per_admunit,
     file = file.path(
       "01_data_cleaned",
       "dataset_occurences.RData")
)

# Explore data ----
test <- data_wallacean_unnested %>%
  filter(year > 1200) %>%
  group_by(speciesKey) %>%
  arrange(year) %>%
  slice(1) %>%  # Pega a primeira linha de cada grupo (ano mais antigo)
  ungroup()

test %>%
  select(scientificName, year, YearOfDescription) %>%
  mutate(diff = year - YearOfDescription ) %>%
  View()

data_wallacean_unnested %>%
  filter(scientificName == "Chiasmocleis albopunctata (Boettger, 1885)") %>%
  View()

tetrapods_polygons_key %>%
  filter(verbatim_name == "Aratinga solstitialis")

# Dados do polígono (assumindo que tetrapods_polygons_key já está carregado)
polygon_data <- tetrapods_polygons_key %>%
  filter(verbatim_name == "Aratinga solstitialis")

# 1. Mapa simples com ggplot2
ggplot() +
  # Adicionar polígono da espécie
  geom_sf(data = polygon_data, 
          fill = "red", 
          alpha = 0.3, 
          color = "darkred", 
          linewidth = 1) +
  # Título
  ggtitle("Distribuição de Aratinga solstitialis") +
  theme_minimal()

# 2. Mapa com fundo do mundo
ggplot() +
  # Adicionar mapa do mundo de fundo
  borders("world", 
          colour = "gray50", 
          fill = "lightgray") +
  # Adicionar polígono da espécie
  geom_sf(data = polygon_data, 
          fill = "red", 
          alpha = 0.5, 
          color = "darkred", 
          linewidth = 1,
          inherit.aes = FALSE) +
  # Limitar ao bounding box do polígono + margem
  coord_sf(xlim = c(-65, -55),  # Extensão leste-oeste
           ylim = c(-5, 7),      # Extensão norte-sul
           expand = TRUE) +
  # Título e labels
  ggtitle("Distribuição de Aratinga solstitialis") +
  labs(subtitle = "Ararajuba - Papagaio-amarelo") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
