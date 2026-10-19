library(tidyverse)
library(sf)
library(arrow)

rm(list=ls()); gc() # clean local enviroment

# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# Load data ----
# occurences data
load(
    file = file.path(
      local_directory,
      "01_data_cleaned",
      "dataset_occurences_220.RData")
)

rm(list = setdiff(ls(), c("local_directory",
                          "list_per_admunit",
                          "data_wallacean_unnested"))); gc()
# trait data
load(file.path(
  local_directory,
  "00_raw_data",
  "trait_data.RData")
)

# socieconomic data
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_temporal.RData")
)
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_humandensity_temporal.RData")
)

# species list
load(
  file.path(
    local_directory,
    "00_raw_data",
    "tetrapods_list_sensitivity_220.RData"
    )
  )

# Completeness list ----
# Amphibia
rich_admunit_expected_amphibia <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Amphibia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_amphibia <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Amphibia") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_amphibia <- left_join(
  rich_admunit_expected_amphibia,
  rich_admunit_observed_amphibia,
  by = "name_en"
)
# Squamata
rich_admunit_expected_reptilia <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Reptilia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_reptilia <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Reptilia") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_reptilia <- left_join(
  rich_admunit_expected_reptilia,
  rich_admunit_observed_reptilia,
  by = "name_en"
)
# Aves
rich_admunit_expected_aves <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Aves") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_aves <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Aves") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_aves <- left_join(
  rich_admunit_expected_aves,
  rich_admunit_observed_aves,
  by = "name_en"
)
# Mammals
rich_admunit_expected_mammalia <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Mammalia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_mammalia <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Mammalia") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_mammalia <- left_join(
  rich_admunit_expected_mammalia,
  rich_admunit_observed_mammalia,
  by = "name_en"
)

save(
  richness_completude_amphibia,
  richness_completude_reptilia,
  richness_completude_aves,
  richness_completude_mammalia,
  file = file.path(
   local_directory,
   "02_data_analysis",
   "richness_completude_220.RData")
)

# Species count per adm. unit ----
list_per_admunit <- species_admin_all %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  select(name_en, all_species) %>%
  unnest(all_species)

sp_count <- list_per_admunit %>% 
  group_by(speciesKey) %>%
  mutate(ua_total = n()) %>%
  distinct(speciesKey, .keep_all = TRUE) %>%
  select(speciesKey, ua_total)
any(is.na(data_wallacean_unnested$year))
# na analise de sensibilidade, checar 2421643

# Create time-to-event table ----
data_wallacean_unnested_modified <- data_wallacean_unnested %>%
  #filter(speciesKey %in% c("2426805")) %>%
  #filter(speciesKey %in% c("2433011", "2423581", "2472164","2426805")) %>% 
  filter(year > 1500) %>% 
  # filter only adm unit with GDP information
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  # spp so pode ser encontrada depois de ja ter sido descrita
  mutate(year_modified = if_else(year < YearOfDescription, YearOfDescription, year)) %>%
  distinct(speciesKey, year_modified, name_en, .keep_all = TRUE) %>%
  relocate(year, YearOfDescription, year_modified, .after = species) %>%
  mutate(
    # Substitui valores ausentes com o mínimo possível válido
    year_modified = as.integer(year_modified),
    month = if_else(is.na(month), 1L, as.integer(month)),
    day = if_else(is.na(day), 1L, as.integer(day)),
    # cria uma data completa
    date = make_datetime(year_modified, month, day),
    date = as.POSIXct(date, tz = "UTC"),
    date = date + seconds(1) # TODO CHECAR SPP QUE SOH TEM ANO NO PRIMEIRO REGISTRO
  ) %>%
  # completeness sensitivity 
  left_join(sp_count, by = "speciesKey") %>% # adm. unit count
  left_join(completeness01, by = "speciesKey") %>%
  left_join(completeness05, by = "speciesKey") %>%
  left_join(completeness1, by = "speciesKey") %>%
  # Join socieconomic data with occurences
  # TODO: aqui tem uma premissa, nao necessariamente a ocorrencia
  # representa o ano que o registro aconteceu. devido a colocarmos o ano de ocorrencia
  # sempre depois da descricao. entao se usarmos "year", teremos NA's porque
  # mudamos o ano de ocorrencia para ser sempre depois da descricao (year_modified).
  left_join(data_socieconomic_temporal, by = c("name_en", "year_modified" = "year")) %>%
  # join human density over time
  left_join(data_humandensity_temporal %>%
              select(speciesKey, year, HumanDensityTime),
            by = c("speciesKey", "year_modified" = "year")) %>%
  select(speciesKey, scientificName, Class, Order, Family, 
         name_en, date, year_modified, log_gdp_final, gdp_rank, gdp_percentile,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus,
         ua_total, ua_total01, ua_total05, ua_total1,
         YearOfDescription, RangeSize, origin_of_data)

data_wallacean_unnested_modified <- data_wallacean_unnested_modified %>%
  group_by(speciesKey) %>%
  #mutate(date = as.Date(date)) %>%
  arrange(date, .by_group = TRUE) %>%
  # REMOVE COMBINATIONS OF SPECIESKEY-ADMUNIT-DATE
  distinct(name_en, date, .keep_all = TRUE) %>%
  # Adiciona secundos sequenciais (0, 1, 2...) para cada data repetida
  # TODO uma premissa chatinha porque algo pode ter vindo antes,
  # como a escala da analise é em anos, acho que acabara ficando diluido
  # ao diferenciar as ocorrencias em segundos!
  mutate(
    date_corrigida = {
      datas <- date
      if (length(datas) > 1) {
        for (i in 2:length(datas)) {
          if (datas[i] <= datas[i - 1]) {
            datas[i] <- datas[i - 1] + seconds(1)
          }
        }
      }
      datas
    }
  ) %>%
  mutate(date = date_corrigida) %>%
  select(-date_corrigida) %>%
  ungroup()

# Create time to event table 
df_wallacean_time <- data_wallacean_unnested_modified %>%
  #filter(ua_total > 2) %>% # pelo menos tres unidades administrativas
  group_by(speciesKey) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    # se o ano da descrição for anterior a 1900, inicia em 1900
    #ano_descricao = if_else(YearOfDescription < 1900, 1900, YearOfDescription),
    # a data do stop corresponde a data da ocorrencia (evento ou censura)
    #t.stop.date = as.Date(date),
    t.stop.date = as.POSIXct(date, tz = "UTC"),
    # a data de start corresponde a ultima data de ocorrencia
    t.start.date = lag(t.stop.date),
    # o primeiro e o ultimo valor sempre terminam com NA, por isso para cada
    # especie, a data de inicio corresponde a descricao da especie - ou seja - 
    # inicia a possibilidade de geolocacao
    t.start.date = if_else(is.na(t.start.date), make_date(YearOfDescription, 1,1), t.start.date),
    t.stop.date = if_else(is.na(t.stop.date), make_date(2026, 1,1), t.stop.date),
    # event
    event = if_else(!duplicated(name_en), 1L, 0L),
    ua_acumulada = cumsum(!duplicated(name_en)),
    # se 1L em event, a wallacecompletude eh junto do ultimo evento
    # se 0L em event, a wallacecompletude eh depois do ultimo evento
    WallaceCompletude = if_else(event == 1L & ua_acumulada == ua_total, 1L, 0L),
    WallaceCompletude01 = if_else(event == 1L & ua_acumulada == ua_total01, 1L, 0L),
    WallaceCompletude05 = if_else(event == 1L & ua_acumulada == ua_total05, 1L, 0L),
    WallaceCompletude1 = if_else(event == 1L & ua_acumulada == ua_total1, 1L, 0L)) %>%
    #Wallace75 = if_else(event == 1L & ua_acumulada >= ua_total * 0.75, 1L, 0L),
    #Wallace50 = if_else(event == 1L & ua_acumulada >= ua_total * 0.5, 1L, 0L)) %>%
  # cria o t.start e o t.stop em dias e transforma em anos
  mutate(
    duracao_evento = as.numeric(difftime(t.stop.date, t.start.date, units = "secs")),
    duracao_evento = if_else(duracao_evento <= 0, 1, duracao_evento),
    t.start = cumsum(lag(duracao_evento, default = 0)) / (86400 * 365.25),
    t.stop  = t.start + duracao_evento / (86400 * 365.25)
  ) %>% 
  mutate(
    t.start.cal = as.numeric(difftime(t.start.date,
                                      as.POSIXct("1756-01-01 00:00:00", tz="UTC"),
                                      units = "secs")) / (365.25 * 24 * 3600),
    t.stop.cal  = as.numeric(difftime(t.stop.date,
                                      as.POSIXct("1756-01-01 00:00:00", tz="UTC"), 
                                      units = "secs")) / (365.25 * 24 * 3600)
  ) %>%
  mutate(
    t.start.year = t.start.cal + 1756,
    t.stop.year  = t.stop.cal  + 1756
  )

df_wallacean_100 <- df_wallacean_time %>%
  filter(!is.na(WallaceCompletude)) %>% # spp sem poligono que entraram aqui
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(WallaceCompletude == 1)) min(linha[WallaceCompletude == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, date,
         t.start.date, t.stop.date, t.start, t.stop, t.start.year, t.stop.year,
         event, WallaceCompletude, log_gdp_final, gdp_rank, gdp_percentile,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus, RangeSize, 
         ua_total, ua_total01, ua_total05, ua_total1, origin_of_data, YearOfDescription) # 24969 spp

df_wallacean_01 <- df_wallacean_time %>%
  filter(!is.na(WallaceCompletude01)) %>% # spp sem poligono que entraram aqui
  # remove com cobertura pequena dos poligonos
  select(-WallaceCompletude, -WallaceCompletude05, -WallaceCompletude1) %>%
  remove_missing() %>%
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(WallaceCompletude01 == 1)) min(linha[WallaceCompletude01 == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, date,
         t.start.date, t.stop.date, t.start, t.stop, t.start.year, t.stop.year,
         event, WallaceCompletude01, log_gdp_final, gdp_rank, gdp_percentile,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus, RangeSize, 
         ua_total01, origin_of_data, YearOfDescription)

df_wallacean_05 <- df_wallacean_time %>%
  filter(!is.na(WallaceCompletude05)) %>% # spp sem poligono que entraram aqui
  # remove com cobertura pequena dos poligonos
  select(-WallaceCompletude, -WallaceCompletude01, -WallaceCompletude1) %>%
  remove_missing() %>%
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(WallaceCompletude05 == 1)) min(linha[WallaceCompletude05 == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, date, 
         t.start.date, t.stop.date, t.start,t.stop, t.start.year, t.stop.year,
         event, WallaceCompletude05, log_gdp_final, gdp_rank, gdp_percentile,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus, RangeSize, 
         ua_total05, origin_of_data, YearOfDescription) 

df_wallacean_1 <- df_wallacean_time %>%
  filter(!is.na(WallaceCompletude1)) %>% # spp sem poligono que entraram aqui
  # remove com cobertura pequena dos poligonos
  select(-WallaceCompletude, -WallaceCompletude01, -WallaceCompletude05) %>%
  remove_missing() %>%
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(WallaceCompletude1 == 1)) min(linha[WallaceCompletude1 == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, date, 
         t.start.date, t.stop.date, t.start,t.stop, t.start.year, t.stop.year,
         event, WallaceCompletude1, log_gdp_final, gdp_rank, gdp_percentile,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus, RangeSize, 
         ua_total1, origin_of_data, YearOfDescription) 

# TODO Precisa pensar em validacoes para atender o modelo de evento
df_wallacean_time %>%
  filter(t.start >= t.stop) %>% View()

table(df_wallacean_time$WallaceCompletude) 
# 1988 completeness sem threshold
# 2546 completeness com 0.1%
# 3000 completeness com 0.5%
# 3294 completeness com 1%

table(df_wallacean_time$WallaceCompletude01) 
table(df_wallacean_time$WallaceCompletude05) 
table(df_wallacean_time$WallaceCompletude1) 


# bons exemplos para checar problemas
# 2433011
# 2423581
# 2472164

save(df_wallacean_100,
     df_wallacean_01,
     df_wallacean_05,
     df_wallacean_1,
    file = file.path(
      local_directory,
      "02_data_analysis",
      "data_wallacean_time_220.RData")
)
