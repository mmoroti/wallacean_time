library(tidyverse)

rm(list=ls()); gc() # clean local enviroment
# TODO:
# Separar por grupos
# incluir dados socieconomicos

# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# Load data ----
# socieconomic data
#load(
#  file = file.path(
#    "01_data_cleaned",
#    "data_socieconomic_cleaned.RData")
#)
                                        
# occurences data
load(
    file = file.path(
      local_directory,
      "dataset_occurences.RData")
)
rm(list = setdiff(ls(), c("local_directory",
                          "list_per_admunit",
                          "data_wallacean_unnested"))); gc()

# Species count per adm. unit ----
sp_count <- list_per_admunit %>% 
  group_by(speciesKey) %>%
  mutate(ua_total = n()) %>%
  distinct(speciesKey, .keep_all = TRUE) %>%
  select(speciesKey, ua_total)
any(is.na(data_wallacean_unnested$year))

# Gerar uma lista de quantos registros sao esperados pelos poligonos
# nas unidades administrativas pelo numero de ocorrencias disponiveis na base
# Global
data_key <- data_wallacean_unnested %>%
  distinct(speciesKey, .keep_all = TRUE) %>%
  select(Class, Order, speciesKey)

list_species <- left_join(list_per_admunit,
          data_key, 
          by = "speciesKey")
# Amphibia
rich_admunit_expected_amphibia <- list_species %>%
  filter(Class == "Amphibia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_amphibia <- data_wallacean_unnested %>%
  filter(Class == "Amphibia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_amphibia <- left_join(
  rich_admunit_expected_amphibia,
  rich_admunit_observed_amphibia,
  by = "name_en"
)
# Squamata
rich_admunit_expected_reptilia <- list_species %>%
  filter(Class == "Reptilia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_reptilia <- data_wallacean_unnested %>%
  filter(Class == "Reptilia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_reptilia <- left_join(
  rich_admunit_expected_reptilia,
  rich_admunit_observed_reptilia,
  by = "name_en"
)
# Squamata
rich_admunit_expected_aves <- list_species %>%
  filter(Class == "Aves") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_aves <- data_wallacean_unnested %>%
  filter(Class == "Aves") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_aves <- left_join(
  rich_admunit_expected_aves,
  rich_admunit_observed_aves,
  by = "name_en"
)
# Mammals
rich_admunit_expected_mammalia <- list_species %>%
  filter(Class == "Mammalia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_mammalia <- data_wallacean_unnested %>%
  filter(Class == "Mammalia") %>%
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
  file = "02_data_analysis/richness_completude.RData"
)

# Create time-to-event table ----
data_wallacean_unnested_modified <- data_wallacean_unnested %>%
  #filter(speciesKey %in% c("2426805")) %>%
  #filter(speciesKey %in% c("2433011", "2423581", "2472164","2426805")) %>% 
  filter(year > 1500) %>% 
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
  left_join(sp_count, by = "speciesKey") %>% # adm. unit count
  # Join socieconomic data with occurences
  #left_join(dados_socieconomic_merge, by = c("sigla_admu", "year")) %>%
  select(speciesKey, scientificName, Class, Order, Family, 
         name_en, date, year_modified,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensity, Latitude, Elevation, ETA50K, AnnuMeanTemp,
         ua_total, YearOfDescription, RangeSize, origin_of_data)

data_wallacean_unnested_modified <- data_wallacean_unnested_modified %>%
  group_by(speciesKey) %>%
  #mutate(date = as.Date(date)) %>%
  arrange(date, .by_group = TRUE) %>%
  distinct(name_en, date, .keep_all = TRUE) %>%
  # Adiciona dias sequenciais (0, 1, 2...) para cada data repetida
  # TODO uma premissa chatinha porque algo pode ter vindo antes,
  # como a escala da analise é em anos, acho que acabara ficando diluido
  # mudar a escala em dias!
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
  filter(ua_total > 2) %>% # pelo menos tres unidades administrativas
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
    Wallace75 = if_else(event == 1L & ua_acumulada >= ua_total * 0.75, 1L, 0L),
    Wallace50 = if_else(event == 1L & ua_acumulada >= ua_total * 0.5, 1L, 0L)) %>%
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
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(WallaceCompletude == 1)) min(linha[WallaceCompletude == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, date,
         t.start.date, t.stop.date, t.start, t.stop, t.start.year, t.stop.year,
         event, WallaceCompletude,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensity, Latitude, Elevation, ETA50K, AnnuMeanTemp, RangeSize, 
         ua_total, origin_of_data, YearOfDescription) 

df_wallacean_75 <- df_wallacean_time %>%
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(Wallace75 == 1)) min(linha[Wallace75 == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, 
         date, t.start.date, t.stop.date, t.start,t.stop, t.start.year, t.stop.year,
         event, Wallace75,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensity, Latitude, Elevation, ETA50K, AnnuMeanTemp, RangeSize, 
         ua_total, origin_of_data, YearOfDescription) 

df_wallacean_50 <- df_wallacean_time %>%
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(Wallace50 == 1)) min(linha[Wallace50 == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, Order, Family, scientificName, name_en, date, 
         t.start.date, t.stop.date, t.start,t.stop, t.start.year, t.stop.year,
         event, Wallace50,
         BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
         HumanDensity, Latitude, Elevation, ETA50K, AnnuMeanTemp, RangeSize, 
         ua_total, origin_of_data, YearOfDescription) 

# TODO precisa checar isso!
# TODO Precisa pensar em validacoes para atender o modelo de evento
df_wallacean_time %>%
  filter(t.start >= t.stop) %>% View()

# bons exemplos para checar problemas
# 2433011
# 2423581
# 2472164
save(df_wallacean_100,
     df_wallacean_75,
     df_wallacean_50,
    file = file.path(
      "02_data_analysis",
      "data_wallacean_time_amphibia.RData")
)
