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

# GBIF FILTER AND CLEAN ----
# load data
start_time <- Sys.time()
load(file.path(
  "00_raw_data",
  "tetrapodstraits_data.RData")
)

# conferir as chaves dos dados baixados
# conferir % das colunas var. resp e preditoras
#head(data_tetrapodstraits)
#table(data_tetrapodstraits$occurrenceStatus) # remover ausencias
#table(data_tetrapodstraits$taxonRank)

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID) 
nrow(data_tetrapodstraits) # 2513477 occurences species
table(data_tetrapodstraits$class)

data_tetrapods_clean <- data_tetrapodstraits %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "gbif") %>%
  mutate(eventDate = na_if(eventDate, "")) %>%
  mutate(class = if_else(class == "", "Squamata", class)) # TODO REVER quando usar dados globais

nrow(data_tetrapods_clean) # 2.489.856 occurences species
table(data_tetrapods_clean$class)

# how much data is na?
data_tetrapods_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
            missing_year = sum(is.na(year)),
            missing_date = sum(is.na(eventDate)))

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
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data =)
            #missing_date = sum(is.na(eventDate))) 

nrow(data_tetrapods_filter) # 1153622 total occurences 
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
                           tests = c("capitals", "centroids", "equal",
                            "gbif","institutions","seas","zeros", "countries")) 
                            # Flagged 58867 of 1153622 records, EQ = 0.05.

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

#end_time <- Sys.time()
#print(end_time - start_time) # Time difference of 28 mins

# BIOTIME FILTER AND CLEAN ----
#start_time <- Sys.time()
load(file.path(
  "00_raw_data",
  "biotime_data.RData")
) 

data_biotime_clean <- biotime_data_key_precleaned %>%
  select("speciesKey", "class", "order", "family", "GENUS_SPECIES", 
         "LATITUDE", "LONGITUDE", "YEAR", "STUDY_ID") %>%
  rename(decimalLatitude = LATITUDE, 
         decimalLongitude = LONGITUDE,
         species = GENUS_SPECIES,
         year = YEAR) %>%
  mutate(origin_of_data = "biotime")

data_biotime_clean %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

flags <- clean_coordinates(x = data_biotime_clean,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("capitals", "centroids", "equal",
                          "gbif","institutions","seas", "zeros")) 
                          # Flagged 397771 of 1982430 records, EQ = 0.2.

#Exclude problematic records
data_biotime_filtered <- data_biotime_clean[
  flags$.summary,
]

#The flagged records
# dat_fl <- data_teste[!flags$.summary,]
# Remove records that are temporal outliers
#flags_temporal <- cf_age(x = data_biotime_filter_spatialpoints,
#                         lon = "decimalLongitude",
#                         lat = "decimalLatitude",
#                         taxon = "species",
#                         min_age = "year",
#                         max_age = "year",
#                         value = "flagged")
#
#data_biotime_filtered <- data_biotime_filter_spatialpoints[
#  flags_temporal,
#] 
min(data_biotime_filtered$year)
max(data_biotime_filtered$year)

#end_time <- Sys.time()
#print(end_time - start_time)

# SPLINK FILTER AND CLEAN ----
#start_time <- Sys.time()
load(file.path(
  "00_raw_data",
  "splinks_data.RData")
) 

# names(splink_data_key_precleaned)
data_splink_clean <- splink_data_key_precleaned %>%
  select("speciesKey", "scientificname",
         "latitude", "longitude", "yearcollected", "coordinateprecision") %>%
  rename(decimalLatitude = latitude, 
         decimalLongitude = longitude,
         species = scientificname,
         year = yearcollected) %>%
  mutate(origin_of_data = "splink") %>%
  mutate(year = na_if(year, "")) %>%
  mutate(coordinateprecision = na_if(coordinateprecision, "")) %>% 
  mutate(year = as.integer(as.character(year))) %>%
  filter(!is.na(year))

data_splink_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)),
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
                           tests = c("capitals", "centroids", "equal",
                            "gbif","institutions","seas","zeros")) # EQ = 0.39.

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
  filter(coordinateprecision / 1000 <= 100 | is.na(coordinateprecision)) %>%
  filter(year > 1899 & year < 2025)

nrow(data_splink_filtered) # 294088

end_time <- Sys.time()
print(end_time - start_time) # 7min57s
# por todas as bases: Time difference of 39.19346 mins

# CHECK DATASETS GBIF, BIOTIME AND SPECIESLINK ----
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

chaves_perdidas_biotime <- setdiff(
  data_biotime_filtered$speciesKey,
  data_tetrapods_filtered$speciesKey
) 

data_biotime_filtered_sa <- data_biotime_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime)) 

chaves_perdidas_splink <- setdiff(
  data_splink_filtered$speciesKey,
  data_tetrapods_filtered$speciesKey
) 

data_splink_filtered_sa <- data_splink_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink))

glimpse(data_splink_filtered_sa)
glimpse(data_biotime_filtered_sa)

data_occurences_precleaned <- bind_rows(
  data_tetrapods_filtered,
  data_biotime_filtered_sa,
  data_splink_filtered_sa)

nested_cols <- c("speciesKey", "class", "order", "family", "species",
  "gbifID", "year","decimalLatitude", "decimalLongitude", "origin_of_data")

data_occurences_precleaned <- data_occurences_precleaned %>% 
  select(all_of(nested_cols))

save(data_occurences_precleaned,  # dados limpos
     file = file.path(
       "01_data_cleaned",
       "data_occurences_filtered.RData"))
# quantidade de registros de ocorrencia por origem dos dados
table(data_occurences_precleaned$origin_of_data)

# REMOVE DUPLICATES AND FILTER POINTS BY RANGE POLYGONS ----
rm(list=ls()); gc() # clean local enviroment
start_time <- Sys.time()
# shapefile
load(file = file.path(
  "00_raw_data",
  "shapefiles_data.RData"))
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

# Remove duplicates
duplicated_flags <- cc_dupl(
  data_occurences_precleaned,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "speciesKey",
  value = "flagged"
) # TODO checar se estamos mantendo o ano mais antigo

data_occurences_precleaned_duplicate <- data_occurences_precleaned[
  duplicated_flags == TRUE, ]

nrow(data_occurences_precleaned_duplicate)*100/nrow(data_occurences_precleaned)
# View(data_occurences_precleaned_duplicate)
length(unique(data_occurences_precleaned_duplicate$speciesKey)) # 7897 spp
# deixar apenas as especies que tem poligonos
data_wallacean_knownledge_sa <- data_occurences_precleaned_duplicate %>%
  filter(speciesKey %in% tetrapod_shapefile$speciesKey)
# deixar apenas os poligonos que tem dados de distribuicao
data_tetrapods_sa <- tetrapod_shapefile %>% 
  filter(speciesKey %in% data_occurences_precleaned_duplicate$speciesKey)
# tem que ter o mesmo numero de especies (speciesKey)
length(unique(data_wallacean_knownledge_sa$speciesKey)) # 6527 spp
length(unique(data_tetrapods_sa$speciesKey)) # 6527 spp

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
end_time <- Sys.time()
print(end_time - start_time) # Time difference of 1h56min

#data_tetrapods_sa_transf <- data_tetrapods_sa %>%
#    filter(speciesKey %in% species_keys_with_errors)

nrow(list_occurences_clean) # 269333 (com filtro temporal 263628)
nrow(data_wallacean_knownledge_sa)-nrow(list_occurences_clean) # 111096 excluidos pelo poligono
length(unique(list_occurences_clean$speciesKey)) # 5748 spp
length(unique(data_wallacean_knownledge_sa$speciesKey)) # 6527 spp

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

# Save 
save(list_occurences_clean,
     species_keys_with_errors,
     file = here(
       "01_data_cleaned",
       "data_occurences_clean.RData")
)
# NESTED DATAFRAME ----
rm(list=ls()); gc() # clean local enviroment
load(file = here(
       "01_data_cleaned",
       "data_occurences_clean.RData")
)
load(file.path(
  "00_raw_data",
  "trait_data.RData")
) # match taxonomic and trait information with speciesKey

tetrapods_key_species <- trait_data %>%
  select("speciesKey","scientificName","Class", "Order","Family", 
  "YearOfDescription","BodyLength_mm","ImputedLength","BodyMass_g",
  "ImputedMass","Diu","Noc","Nocturnality", "Fos", "Ter", "Aqu", "Arb", "Aer",
  "ImputedHabitat","MajorHabitatSum","ImputedMajorHabitat","RangeSize",
   "HumanDensity", "AssessedStatus") %>%
  distinct(speciesKey, .keep_all = TRUE) # 117 speciesKey duplicated
anyDuplicated(tetrapods_key_species$speciesKey) # ok

data_tetrapods_nested <- list_occurences_clean %>%
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
  arrange(Class, Order, Family)

head(data_wallacean_nested)

# check data
table(data_wallacean_nested$Class)
# Amphibia     Aves   Mammalia Reptilia
#     789      2765      893     1279
table(data_wallacean_nested$Order)

anyDuplicated(data_wallacean_nested$speciesKey)

for (i in 1:nrow(data_wallacean_nested)) {
  data_wallacean_nested[i,"count_events"] <- nrow(
    data_wallacean_nested$event_table[[i]])
  
  data_wallacean_nested$event_table[[i]] <- 
    data_wallacean_nested$event_table[[i]][
      order(data_wallacean_nested$event_table[[i]]$year), 
    ]
}

save(data_wallacean_nested, 
     file = here(
       "01_data_cleaned",
       "data_wallacean_nested.RData")
)

# Explorar os dados
ggplot(data_wallacean_nested, aes(x = log(RangeSize), y = log(count_events))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Adiciona uma linha de tendência
  labs(x = "Comprimento do Corpo (mm)", y = "Número de Registros (count_events)", title = "Relação entre Comprimento do Corpo e Registros") +
  theme_minimal()

ggplot(data_wallacean_nested, aes(x = AssessedStatus, y = log(count_events))) +
  geom_boxplot(fill = "lightblue", colour = "darkblue") +
  theme_minimal() +
  labs(title = "Comparação de quantidade de registros por categoria da IUCN",
       x = "Categoria da IUCN",
       y = "Quantidade de registros (event_counts)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
