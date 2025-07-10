library(tidyverse)

rm(list=ls()); gc() # clean local enviroment
# socieconomic data
load(
  file = file.path(
    "01_data_cleaned",
    "data_socieconomic_cleaned.RData")
)
# occurences data
load(
    file = file.path(
      "01_data_cleaned",
      "dataset_occurences.RData")
)
# list per adm unit
load(file = file.path(
       "01_data_cleaned",
       "data_occurences_geolocation.RData")
)

# Gerando lista de quantas und administrativas cada especie pode ocupar
sp_count <- list_per_admunit %>% 
  group_by(speciesKey) %>%
  mutate(ua_total = n()) %>%
  distinct(speciesKey, .keep_all = TRUE) %>%
  select(speciesKey, ua_total)

# Change data to discovery
data_wallacean_unnested_modified <- data_wallacean_unnested %>%
  # especie so pode ser encontrada depois de ja ter sido descrita
  mutate(year_modified = if_else(year < YearOfDescription, YearOfDescription, year)) %>%
  relocate(year, YearOfDescription, year_modified, .after = species) %>%
  mutate(
    # Substitui valores ausentes com o mínimo possível válido
    year_modified = as.integer(year_modified),
    month = if_else(is.na(month), 1L, as.integer(month)),
    day = if_else(is.na(day), 1L, as.integer(day)),
    
    # cria uma data completa
    date = make_date(year_modified, month, day),
    date = date + days(1) # TODO CHECAR SPP QUE SOH TEM ANO NO PRIMEIRO REGISTRO
  ) %>%
  # TODO CHECAR
  # removendo datas repetidas para specieskey,datas e admu repetidas
  arrange(speciesKey, date) %>%
  filter(!duplicated(paste(speciesKey, date, sigla_admu))) 

# Join socieconomic data with occurences
data_wallacean <- data_wallacean_unnested_modified %>%
  left_join(sp_count, by = "speciesKey") %>% # adm. unit count
  left_join(dados_socieconomic_merge, by = c("sigla_admu", "year")) %>% #socioeconomic data
  select(speciesKey, scientificName, Class, sigla_admu, date, 
  BodyLength_mm, Nocturnality, Verticality, ua_total, YearOfDescription,
  'pop_density(person/km2)','pib_dollar_current') %>% # id, time, cov., adm. unit count
  rename('pop_density' = 'pop_density(person/km2)') 

#df_eventos <- data_wallacean %>%
#  filter(ua_total > 2) %>% # spp esperadas em + de 2 unid. adm.
#  group_by(speciesKey) %>%
#  mutate(
#    year_relative = year_modified - min(year_modified),
#    t.start = as.numeric(year_relative),
#    t.stop  = lead(t.start),
#    t.stop = ifelse(is.na(t.stop), 2025 - min(year_modified), t.stop),
#    event   = ifelse(!duplicated(sigla_admu), 1, 0), 
#    ua_acumulada = cumsum(!duplicated(sigla_admu)),
#    WallaceCompletude = ifelse(event == 1 & ua_acumulada == ua_total, 1, 0),
    #Wallace50 = ifelse(event == 1 & ua_acumulada => ua_total/2, 1, 0),
#    event = as.integer(event),
#    WallaceCompletude = as.integer(WallaceCompletude)
#  ) %>%
#  ungroup() %>%
#  select(speciesKey, scientificName, Diu, Noc, BodyLength_mm,
#         t.start, t.stop, event, WallaceCompletude) #%>%
# #filter(speciesKey %in% c("2424006", "2428622", "2422148"))

# Create time to event table 
df_wallacean_time <- data_wallacean %>%
  filter(ua_total > 2) %>% # pelo menos duas unidades administrativas
  group_by(speciesKey) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
  # se o ano da descrição for anterior a 1900, inicia em 1900
    ano_descricao = if_else(YearOfDescription < 1900, 1900, YearOfDescription),
  # a data do stop corresponde a data da ocorrencia (evento ou censura)
    t.stop.date = as.Date(date),
  # a data de start corresponde a ultima data de ocorrencia
    t.start.date = lag(t.stop.date),
  # o primeiro e o ultimo valor sempre terminam com NA, por isso para cada
  # especie, a data de inicio corresponde a descricao da especie - ou seja - 
  # inicia a possibilidade de geolocacao
    t.start.date = if_else(is.na(t.start.date), make_date(ano_descricao, 1,1), t.start.date),
    t.stop.date = if_else(is.na(t.stop.date), make_date(2026, 1,1), t.stop.date)) %>%
  # cria o t.start e o t.stop em dias
  mutate(
    duracao_evento = as.numeric(difftime(t.stop.date, t.start.date, units = "days")),
    t.start = cumsum(lag(duracao_evento, default = 0))/365.25,
    t.stop  = t.start + duracao_evento/365.25,
    event = if_else(!duplicated(sigla_admu), 1L, 0L),
    ua_acumulada = cumsum(!duplicated(sigla_admu)),
    # se 1L em event, a wallacecompletude eh junto do ultimo evento
    # se 0L em event, a wallacecompletude eh depois do ultimo evento
    WallaceCompletude = if_else(event == 0L & ua_acumulada == ua_total, 1L, 0L) 
  ) %>%
  # remove todos eventos apos completude 
  mutate(
    linha = row_number(),
    pos_completude = if (any(WallaceCompletude == 1)) min(linha[WallaceCompletude == 1]) else Inf
  ) %>%
  filter(linha <= pos_completude) %>%
  ungroup() %>%
  select(speciesKey, Class, scientificName, 
         t.start.date, t.stop.date, t.start,t.stop, event, WallaceCompletude,
         Nocturnality, Verticality, BodyLength_mm) 

# TODO precisa checar isso!
df_wallacean_time %>%
  filter(t.start == t.stop & event == 1) %>%
  View()

# um bom exemplo para checar problemas
# 2433011 

save(df_wallacean_time,
    file = file.path(
      "02_data_analysis",
      "data_wallacean_time.RData")
)