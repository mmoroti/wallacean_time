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
  "here",
  "CoordinateCleaner",
  "sf",
  "countrycode"
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
rm(list = setdiff(ls(), "local_directory")); gc()

load(file.path(
  local_directory,
  "tetrapodstraits_data.RData")
)

# conferir as chaves dos dados baixados
# conferir % das colunas var. resp e preditoras
#head(data_tetrapodstraits)
#table(data_tetrapodstraits$occurrenceStatus) # remover ausencias
#table(data_tetrapodstraits$taxonRank)

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID) 
nrow(data_tetrapodstraits) 
# 24.441.908 occurences species Global
# 2.501.450 occurences species America do Sul
table(data_tetrapodstraits$class)

data_tetrapods_clean <- data_tetrapodstraits %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "gbif") %>%
  mutate(eventDate = na_if(eventDate, "")) %>%
  # U. ricardinii aparece sem classificacao de classe
  mutate(class = if_else(class == "", "Squamata", class))

nrow(data_tetrapods_clean) 
# 2.478.512 occurences species America do Sul
# 24.263.774

table(data_tetrapods_clean$class)

# how much data is na?
data_tetrapods_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year)),
            missing_date = sum(is.na(eventDate)))
# 8.916.385 sem coordenada
# 3.568.206 sem ano 

# exists a difference between "year" and "event date" of 7280 occurences
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
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), # 112.301
            missing_day = sum(is.na(day)), # 62.065 NA's
            missing_month = sum(is.na(month)), # 62.065 NA's
            missing_year = sum(is.na(year))) # without missing data =)
            #missing_date = sum(is.na(eventDate))) 

nrow(data_tetrapods_filter) 
# 1.153.622 total occurences america do sul
# 14.108.708 global
table(data_tetrapods_filter$class) # occurences per classes 

# coordsclean
# convert country code from ISO2c to ISO3c
data_tetrapods_filter$countryCode <-  countrycode(
  data_tetrapods_filter$countryCode,
  origin =  'iso2c',
  destination = 'iso3c')

flags <- clean_coordinates(x = data_tetrapods_filter,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros"))

#Flagged 674578 of 14108708 records, EQ = 0.05.
#Exclude problematic records
data_tetrapods_filter_spatialpoints <- data_tetrapods_filter[
  flags$.summary,
  ]

#The flagged records
# dat_fl <- data_teste[!flags$.summary,]

# Remove records that are temporal outliers
#flags_temporal <- cf_age(x = data_tetrapods_filter_spatialpoints,
#                lon = "decimalLongitude",
#                lat = "decimalLatitude",
#                taxon = "species", 
#                min_age = "year",
#                max_age = "year",
#                value = "flagged")
#
#data_tetrapods_filter_spatial_temp <- data_tetrapods_filter_spatialpoints[
#  flags_temporal,
#  ]

# Remove records with low coordinate precision
#data_tetrapods_filter_spatial_temp %>% 
#  mutate(Uncertainty = coordinateUncertaintyInMeters / 1000) %>% 
#  ggplot(aes(x = Uncertainty)) + 
#  geom_histogram() +
#  xlab("Coordinate uncertainty in meters") +
#  theme_bw()

# filter coordinates uncertain > 100km
# filter years from 1900
data_tetrapods_filtered <- data_tetrapods_filter_spatialpoints %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(year >1899)

## BioTIME 2.0v ----
#start_time <- Sys.time()
load(file.path(
  local_directory,
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
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year))) # without missing data

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

min(data_biotime_filtered$year)
max(data_biotime_filtered$year)

## speciesLink ----
load(file.path(
  local_directory,
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
    day = na_if(day, ""),
    month = na_if(month, ""),
    year = na_if(year, "")) %>%
  mutate(coordinateprecision = na_if(coordinateprecision, "")) %>% 
  mutate(year = as.integer(as.character(year))) %>%
  filter(!is.na(year))

data_splink_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_day = sum(is.na(day)),
            missing_month = sum(is.na(month)),
            missing_year = sum(is.na(year)))

# invalidity coordinates
vetor <- c(19875, 20466, 20880, 21001, 28842, 30252, 33789, 34402, 34404, 34407, 
            34409, 34411, 34413, 34415, 34417, 34418, 34474, 34493, 34495, 34497, 
            34499, 34533, 34534, 50880, 169018, 170384, 170499, 171015, 171016, 
            211761, 211782, 211828, 317952, 344767, 410361, 411544, 411809, 411936, 
            412107, 412108, 413889, 414494, 422455, 424351, 424529, 424546, 426406, 
            426407, 428892, 429413, 429755, 429868, 430062, 445619, 445625, 445711, 
            445712, 445713, 445714, 445715, 445716, 445717, 445718, 445719, 445720, 
            445721, 445722, 445723, 445724, 445725, 445726, 583577)

data_splink_cleaned <-  data_splink_clean[-vetor, ]

flags <- clean_coordinates(x = data_splink_cleaned,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("equal","gbif",
                                     "institutions",
                                     "seas","zeros")) 
# Flagged 225776 of 597492 records, EQ = 0.38.
#Exclude problematic records
data_splink_filter_spatialpoints <- data_splink_cleaned[
  flags$.summary,
]

#The flagged records
# dat_fl <- data_teste[!flags$.summary,]

## Remove records that are temporal outliers
#flags_temporal <- cf_age(x = data_splink_filter_spatialpoints,
#                         lon = "decimalLongitude",
#                         lat = "decimalLatitude",
#                         taxon = "species",
#                         min_age = "year",
#                         max_age = "year",
#                         value = "flagged") # Flagged 16743 records
#
#data_splink_filter_spatial_temp <- data_splink_filter_spatialpoints[
#  flags_temporal,
#]

# filter coordinates uncertain > 100km
# filter years from 1900
data_splink_filtered <- data_splink_filter_spatialpoints %>%
  mutate(coordinateprecision = as.numeric(coordinateprecision)) %>%
  filter(coordinateprecision / 1000 <= 100 | is.na(coordinateprecision)) %>%
  filter(year > 1899 & year < 2025)

nrow(data_splink_filter_spatialpoints)-nrow(data_splink_filtered)

end_time <- Sys.time()
print(end_time - start_time) # 7min57s
# por todas as bases: Time difference of 39.19346 mins

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

# If data its ok, save data
save(
  data_tetrapods_filtered,
  data_biotime_filtered,
  data_splink_filtered, 
  file = file.path(
    "01_data_cleaned",
    "datasets_filtered.RData")
)

# UNIFYING OCCURENCES DATA ---- 
# TODO aqui vai ser necessario checar quando usarmos dados globais
# a maioria sao de especies fora da america do sul
# por isso estamos perdendo tantos dados aqui
# quando for para fazer para o mundo sera necessario reconferir
#key_losing <- left_join(data.frame(speciesKey = chaves), 
#                        data_biotime_clean_v1,
#                        by = "speciesKey")
# load data
rm(list=ls()); gc() # clean local enviroment
load(file.path(
  "01_data_cleaned",
  "datasets_filtered.RData")
)
# TODO da para usar para baixar o registro de occ para mais especies com essas 
# chaves, e buscar + poligonos se necessarios. 
chaves_perdidas_biotime <- setdiff(
  data_biotime_filtered$speciesKey,
  data_tetrapods_filtered$speciesKey
) 

data_biotime_filtered_sa <- data_biotime_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month))

chaves_perdidas_splink <- setdiff(
  data_splink_filtered$speciesKey,
  data_tetrapods_filtered$speciesKey
) 

data_splink_filtered_sa <- data_splink_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink)) %>%
  mutate(day = as.integer(day),
         month = as.integer(month))

glimpse(data_splink_filtered_sa)
glimpse(data_biotime_filtered_sa)

data_occurences_precleaned <- bind_rows(
  data_tetrapods_filtered,
  data_biotime_filtered_sa,
  data_splink_filtered_sa)

nested_cols <- c("speciesKey", "class", "order", "family", "species",
  "gbifID", "day", "month", "year","decimalLatitude", "decimalLongitude", "origin_of_data")

data_occurences_precleaned <- data_occurences_precleaned %>% 
  select(all_of(nested_cols))

biotime_species <- unique(data_occurences_precleaned$speciesKey[data_occurences_precleaned$origin_of_data == "biotime"])
splink_species <- unique(data_occurences_precleaned$speciesKey[data_occurences_precleaned$origin_of_data == "splink"])
gbif_species <- unique(data_occurences_precleaned$speciesKey[data_occurences_precleaned$origin_of_data == "gbif"])
# spp em biotime que não estão no GBIF
setdiff(biotime_species, gbif_species)
# splink que não estão no GBIF
setdiff(splink_species, gbif_species)

save(data_occurences_precleaned,  # dados limpos e integrados
     file = file.path(
       "01_data_cleaned",
       "data_occurences_filtered.RData"))

# REMOVE DUPLICATES ----
rm(list=ls()); gc() # clean local enviroment
start_time <- Sys.time()
# occurences data
load(file = file.path(
  "01_data_cleaned",
  "data_occurences_filtered.RData"))

# plot points without polygons filter
#ggplot() +
#  coord_fixed() +
#  borders("world", colour = "gray50", fill = "gray50") +
#  geom_point(data = data_occurences_precleaned,
#             aes(x = decimalLongitude, y = decimalLatitude),
#             colour = "darkred",
#             size = 0.5) +
#  theme_bw()

data_occurences_precleaned <- data_occurences_precleaned %>%
  mutate(
    # duas casas decimais = 1.11 km 
    decimalLatitude = round(decimalLatitude, 2), 
    decimalLongitude = round(decimalLongitude, 2),
    month = ifelse(month > 12 | month < 1, NA, month),
    day = ifelse(day > 31 | day < 1, NA, day)          
  ) %>%
  arrange(speciesKey, year, month, day) 

# Remove duplicates
duplicated_flags <- cc_dupl(
  data_occurences_precleaned,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "speciesKey",
  value = "flagged"
) # TODO conferir ordenacao se ta certinho, mas ja havia feito isso

data_occurences_precleaned_duplicate <- data_occurences_precleaned[
  duplicated_flags == TRUE, ]

nrow(data_occurences_precleaned_duplicate)*100/nrow(data_occurences_precleaned)
# View(data_occurences_precleaned_duplicate)
length(unique(data_occurences_precleaned_duplicate$speciesKey)) 
# 7691 south america spp
# 25928 global scale spp

# FILTER POINTS BY RANGE POLYGONS ----
# shapefile
load(file = file.path(
  "00_raw_data",
  "tetrapods_polygons_key.RData")) # From TetrapodTraits

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
length(unique(data_tetrapods_sa$speciesKey)) # perdendo 3 spp sem poligono

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
nrow(list_occurences_clean) 
nrow(data_wallacean_knownledge_sa)-nrow(list_occurences_clean) 
length(unique(list_occurences_clean$speciesKey)) 
length(unique(data_wallacean_knownledge_sa$speciesKey))

View(data_wallacean_knownledge_sa)
# Quantas espécies para cada grupo?
list_species <- list_occurences_clean %>%
 distinct(speciesKey, .keep_all=TRUE)
table(list_species$class)

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
# shapefile to crop adm unit
south_america <- st_read(file.path(
  "Shapefiles",
  "GEOGRAPHIC_SHAPEFILE",
  "south_america_br_states.shp"
))

load(file.path("00_raw_data",
       "geographic_shape_data.RData"))

# Verificar validade das geometrias
sf_use_s2(FALSE)
st_is_valid(data_tetrapods_sa) %>% table()
st_is_valid(geographic_shape_data) %>% table()

# extracting geolocation per specie
data_occurences <- st_as_sf(list_occurences_clean,
                   coords = c("decimalLongitude", "decimalLatitude"),
                   crs = st_crs(geographic_shape_data))

# com os poligonos de especialistas tinhamos 269333
# usando o tetrapodtraits aumentamos para 405215 occ
nrow(data_occurences) 
nrow(list_occurences_clean) 

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
# 62.956 occ registros fora do poligonos (oceano, areas em litigio, etc)
nrow(data_occurences)-nrow(data_occurences_geo) # 21821 occ removidas 
names(data_occurences_geo)

data_occurences_geometry <- data_occurences_geo %>%
  select(speciesKey, species,day,month, year, origin_of_data, name_en, ID) #%>%
  #rename(codeAdmUnit = adm0_a3, adm_unit = name_en)

data_occurences_units <- data_occurences_geometry %>%
  st_drop_geometry()

#mapview(geographic_shape_data,
# col.regions = "lightblue", alpha.regions = 0.4, color = "black", layer.name = "Regiões") +
#  mapview(data_occurences_geo, 
#  color = "red", cex = 2, layer.name = "Ocorrências", popup = popupTable(data_occurences_geo, zcol = "name_en"))

# SPECIES LIST PER ADMINASTRIVE UNIT ----
interseccao_sf <- st_intersection(data_tetrapods_sa, geographic_shape_data)

View(geographic_shape_data)
list_per_admunit <- interseccao_sf %>%
  st_drop_geometry()

list_per_admunit$name_en %>% table()

save(list_occurences_clean, # sem duplicatas e filtrados por poligonos
     data_occurences_units, # com as unidades administrativas
     list_per_admunit,      # lista de especies por und adm
     file = file.path(
       "01_data_cleaned",
       "data_occurences_geolocation.RData")
)

# NESTED DATAFRAME WITH BIOLOGICAL TRAITS ----
rm(list=ls()); gc() # clean local enviroment

load(file = file.path(
       "01_data_cleaned",
       "data_occurences_geolocation.RData")
)

load(file.path(
  "00_raw_data",
  "trait_data.RData")
) # match taxonomic and trait information with speciesKey

tetrapods_key_species <- trait_data %>%
  select("speciesKey","scientificName","Class", "Order","Family", 
  "YearOfDescription","BodyLength_mm","ImputedLength","BodyMass_g",
  "ImputedMass","Diu","Noc","Nocturnality", "Fos", "Ter", "Aqu", "Arb", "Aer",
  "ImputedHabitat","Verticality","MajorHabitatSum","ImputedMajorHabitat","RangeSize",
   "HumanDensity", "AssessedStatus") %>%
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
  mutate(count_events = NA) %>%
  relocate(event_table,count_events, .after = Family) %>%
  arrange(Class, Order, scientificName)

View(data_wallacean_nested)
# check data
table(data_wallacean_nested$Class)
# com os poligonos de especialistas
# Amphibia     Aves   Mammalia Reptilia
#     797      2764      898     1281
# com o tetrapodtraits
# Amphibia     Aves Mammalia Reptilia
#    1825     2829     1024     1367
# para o mundo
#Amphibia  Aves   Mammalia Reptilia 
# 4942     8250     4459     7050 

table(data_wallacean_nested$Order)
anyDuplicated(data_wallacean_nested$speciesKey)

# adicionar contagem de unidades administrativas diferentes
#for (i in 1:nrow(data_wallacean_nested)) {
#  data_wallacean_nested[i,"count_events"] <- nrow(
#    data_wallacean_nested$event_table[[i]])
#  
#  data_wallacean_nested$event_table[[i]] <- 
#    data_wallacean_nested$event_table[[i]][
#      order(data_wallacean_nested$event_table[[i]]$year), 
#    ]
#}

# unnest data
data_wallacean_unnested <- data_wallacean_nested %>%
  unnest(cols = c(event_table))
nrow(data_wallacean_unnested) 
# 342.259  occ south america
# 3.144.303 occ global

save(data_wallacean_nested,
    data_wallacean_unnested,
    file = here(
      "01_data_cleaned",
      "dataset_occurences.RData")
)
