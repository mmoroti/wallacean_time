# TODO adaptar para usar o ambiente virtual
# Functions and packages ----
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
    labs(title = title, y = "Frequency of occurrences")
}

install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    renv::restore()
    #install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "cowplot",
  "tidyverse",
  "here"
)

# executing install & load for each package
sapply(package_vec, install.load.package)

# LOADING DATASETS ----
rm(list=ls()); gc() # clean local enviroment

# load data
load(file.path(
  "00_raw_data",
  "tetrapodstraits_data.RData")
) 
load(file.path(
  "00_raw_data",
  "biotime_data.RData")
)
load(file.path(
  "00_raw_data",
  "splinks_data.RData")
)

# GBIF FILTER AND CLEAN ----
# conferir as chaves dos dados baixados
# conferir % das colunas var. resp e preditoras
head(data_tetrapodstraits)
table(data_tetrapodstraits$occurrenceStatus) # remover ausencias
table(data_tetrapodstraits$taxonRank)

View(data_tetrapodstraits %>% 
       filter(taxonRank == "UNRANKED"))

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID) 

data_tetrapods_clean <- data_tetrapodstraits %>%
  select("speciesKey", "class", "order", "family", "species",
         #"countryCode", "locality", "stateProvince",
         # "coordinateUncertaintyInMeters", "coordinatePrecision",
         "decimalLatitude", "decimalLongitude", "eventDate", "year",
         "taxonRank", "occurrenceStatus", "gbifID") %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  select(-"occurrenceStatus", -"taxonRank") %>%
  mutate(origin_of_data = "gbif") %>%
  mutate(eventDate = na_if(eventDate, ""))

nrow(data_tetrapods_clean) # 2.489.856 occurences species
table(data_tetrapods_clean$class)

# TODO rever quando usar dados globais
#View(data_tetrapods_clean_v1 %>% 
#       filter(class == "")) # todas que nao tem classe correspondem ao mesmo sp
data_tetrapods_clean <- data_tetrapods_clean %>%
  mutate(class = if_else(class == "", "Squamata", class)) 

table(data_tetrapods_clean$class)
glimpse(data_tetrapods_clean)
names(data_tetrapods_clean)

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

# per species
length(unique(data_tetrapods_filter$speciesKey)) # 8122 spp
count_species <- data_tetrapods_filter %>%
  distinct(species, .keep_all = TRUE)
table(count_species$class)

# Comparative plots of original vs filtered GBIF datasets
comparative_plot(data_tetrapods_clean,
                 data_tetrapods_filter,
                 "year",
                 title = "GBIF dataset")

# BIOTIME FILTER AND CLEAN ----
data_biotime_clean <- biotime_data_key_precleaned %>%
  select("speciesKey", "class", "order", "family", "GENUS_SPECIES", 
         "LATITUDE", "LONGITUDE", "YEAR", "STUDY_ID") %>%
  rename(decimalLatitude = LATITUDE, 
         decimalLongitude = LONGITUDE,
         species = GENUS_SPECIES,
         year = YEAR) %>%
  mutate(origin_of_data = "biotime")

names(data_biotime_clean)

data_biotime_clean %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

summary(data_biotime_clean$decimalLatitude)

# plot time series of biotime
df_filter <- as.data.frame(table(data_biotime_clean$year)) %>%
  rename(Year = Var1) %>%
  mutate(Year = as.numeric(as.character(Year)))

# Criar o gráfico comparativo
ggplot(df_filter, aes(x = Year, y = Freq)) +
  geom_bar(stat = "identity", position = "dodge", width = 2) +
  geom_line(stat = "smooth", linewidth = 1.0, linetype = "dashed") + 
  theme_bw() +
  labs(title = "BioTime", y = "Frequency of occurrences")

# nenhuma coordenada zerada
#nrow(data_biotime_clean %>%
#       filter(decimalLatitude == 0 | decimalLongitude == 0))

# SPLINK FILTER AND CLEAN ----
names(splink_data_key_precleaned)

data_splink_clean <- splink_data_key_precleaned %>%
  select("speciesKey", "taxonclass", "ordem", "family", "scientificname", 
         "latitude", "longitude", "yearcollected") %>%
  rename(decimalLatitude = latitude, 
         decimalLongitude = longitude,
         class = taxonclass,
         order = ordem,
         species = scientificname,
         year = yearcollected) %>%
  mutate(origin_of_data = "splink") %>%
  mutate(year = na_if(year, "")) %>%
  mutate(year = as.integer(as.character(year))) 
  
names(data_splink_clean)
nrow(data_splink_clean) # 660.781 ocorrencias

# alguns casos tem 0 em uma das coordenadas
# latitude com 0
nrow(data_splink_clean %>%
       filter(decimalLatitude == 0 & decimalLongitude != 0)) 
# longitude com 0
nrow(data_splink_clean %>%
       filter(decimalLatitude != 0 & decimalLongitude == 0)) 

nrow(data_splink_clean %>%
  filter(decimalLatitude == 0 & decimalLongitude == 0)) # coordenadas com zero
# ocorrencias de coordenadas e longitudes iguais a zero = 243.006

nrow(data_splink_clean %>% 
  filter(is.na(year))) # 63.167 anos com NA

nrow(data_splink_clean %>%
       filter(year == 0)) # 605 anos com zero

nrow(data_splink_clean %>%
       filter(year > 2024)) # 6 registros depois de 2024

nrow(data_splink_clean %>%
       filter(year < 1756)) # 1017 nomes antes de linneu, muitos 1111 e 0

nrow(data_splink_clean %>%
  filter(year > 1755 & year <2025)) # 596.591

# remove coordenadas com 0 na lat e na long
# filtra dados sem ano
data_splink_filter <- data_splink_clean %>%
  filter(!(decimalLatitude == 0 & decimalLongitude == 0)) %>% 
  filter(!is.na(year)) %>%
  filter(year > 1755 & year <2025)

nrow(data_splink_clean) - nrow(data_splink_filter) # 269.489 occ removed

data_splink_filter %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

# the wrongs years difficult to visualize temporal series tendence 
data_splink_clean_compare <- data_splink_clean %>%
  filter(!is.na(year)) %>%
  filter(year > 1755 & year < 2025)

comparative_plot(data_splink_clean_compare, data_splink_filter,
                "year", title = "splink dataset")

save(
  data_tetrapods_filter,
  data_biotime_clean, 
  data_splink_filter,
  file = file.path(
    "01_data_cleaned",
    "occurences_filter_data.RData")
  )

# UNIFYING OCCURENCES DATA ---- 
# TODO aqui vai ser necessario checar quando usarmos dados globais
# a maioria sao de especies fora da america do sul
# por isso estamos perdendo tantos dados aqui
# quando for para fazer para o mundo sera necessario reconferir
#key_losing <- left_join(data.frame(speciesKey = chaves), 
#                        data_biotime_clean_v1,
#                        by = "speciesKey")
#
#length(unique(data_tetrapods_filter$speciesKey))
#length(unique(data_biotime_clean_v1$speciesKey))
#length(chaves)
#nrow(key_losing) # 1022095 occurences
#View(key_losing %>%
#  distinct(speciesKey, .keep_all = TRUE))
rm(list=ls()); gc() # clean local enviroment

# load data
load(file.path(
  "01_data_cleaned",
  "occurences_filter_data.RData")
) 

chaves_perdidas_biotime <- setdiff(
  data_biotime_clean$speciesKey,
  data_tetrapods_filter$speciesKey
) 

data_biotime_clean_sa <- data_biotime_clean %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime)) 

chaves_perdidas_splink <- setdiff(
  data_splink_filter$speciesKey,
  data_tetrapods_filter$speciesKey
) 

data_splink_clean_sa <- data_splink_filter %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink))

glimpse(data_splink_clean_sa)
glimpse(data_biotime_clean_sa)

data_wallacean_knownledge <- bind_rows(
  data_tetrapods_filter,
  data_biotime_clean_sa,
  data_splink_clean_sa)

table(data_wallacean_knownledge$origin_of_data)

# plot time series
df_count <- data_wallacean_knownledge %>%
  group_by(year, origin_of_data) %>%
  summarise(Occurences = n()) %>%
  ungroup() %>%
  filter(year >1899)

df_total <- df_count %>%
  group_by(year) %>%
  summarise(Total_Occurences = sum(Occurences)) %>%
  filter(year >1899)

ggplot() +
  geom_bar(data = df_count, aes(x = year, y = Occurences, fill = origin_of_data), 
           stat = "identity", position = "stack", width = 0.7, alpha = 0.4) +  # Barras empilhadas
  geom_line(data = df_count, aes(x = year, y = Occurences, color = origin_of_data, group = origin_of_data), 
            size = 1.0, linetype = "dashed", stat = "smooth") +  # Linhas para cada origem de dados
  #geom_smooth(data = df_total, aes(x = year, y = Total_Occurences), 
  #            method = "loess", se = FALSE, color = "black", size = 1.0, linetype = "dashed") +
  theme_bw() +
  labs(title = "Série Temporal de ocorrências por origem dos dados", 
       x = "Ano", 
       y = "Número de Ocorrências",
       fill = "Origem dos Dados",
       color = "Origem dos Dados") +  # Sincroniza legendas
  scale_fill_brewer(palette = "Set1") +  # Paleta para barras
  scale_color_brewer(palette = "Set1")

nrow(data_wallacean_knownledge) # 2497990 total occurences 
table(data_wallacean_knownledge$class) # occurences per classes 

# per species
length(unique(data_wallacean_knownledge$speciesKey)) # 8122 spp
count_species <- data_tetrapods_filter %>%
  distinct(species, .keep_all = TRUE)
table(count_species$class)

#---- COORDSCLEAN
# TODO executar o coordsclean e limpar coordenadas erradas

# NESTED DATAFRAME
data_tetrapods_nested <- data_wallacean_knownledge %>%
  group_by(speciesKey) %>%  # Agrupa pelo speciesKey
  nest() %>%
  rename(event_table = data) 

anyDuplicated(data_tetrapods_nested$speciesKey) # ok

tetrapods_key_species <- data_tetrapods_filter %>% 
  select(speciesKey, species, class, order, family) %>%
  distinct(speciesKey, .keep_all = TRUE)

data_wallacean_nested_v1 <- left_join(
  data_tetrapods_nested,
  tetrapods_key_species,
  by = "speciesKey"
)

data_wallacean_nested <-  data_wallacean_nested_v1 %>% 
  relocate(event_table, .after = family) %>%
  arrange(class, order, family) %>%
  mutate(count_events = NA)

head(data_wallacean_nested)

# check data
table(data_wallacean_nested$class)
# Amphibia   Aves     Crocodylia   Mammalia   Squamata Testudines 
# 2065       3230          8       1183       1571         65
anyDuplicated(data_wallacean_nested$species) # nenhum nome repetido
anyDuplicated(data_wallacean_nested$speciesKey)

# TODO adicionar a contagem por origem dos dados por especie
for (i in 1:nrow(data_wallacean_nested)) {
  data_wallacean_nested[i,"count_events"] <- nrow(
    data_wallacean_nested$event_table[[i]])
  
  data_wallacean_nested$event_table[[i]] <- 
    data_wallacean_nested$event_table[[i]][
      order(data_wallacean_nested$event_table[[i]]$year), 
    ]
}

data_wallacean_nested <- data_wallacean_nested %>% 
       relocate(count_events, .before = event_table) 

save(data_wallacean_nested, 
     file = here(
       "01_data_cleaned",
       "data_wallacean_nested.RData")
)
