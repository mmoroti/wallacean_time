library(sf)
library(tidyverse)

south_america <- st_read(file.path(
  "Shapefiles",
  "GEOGRAPHIC_SHAPEFILE",
  "south_america_br_states.shp"
)) %>%
  mutate(sigla_admu = case_when(
    sigla_admu == "URU" ~ "URY",
    sigla_admu == "PAR" ~ "PRY",
    TRUE ~ sigla_admu
  ))

adm_unit_list <- south_america %>%
  select(adm_unit,sigla_admu) %>%
  st_drop_geometry()

# unidades administrativas nao encontradas
setdiff(adm_unit_list$sigla_admu,
   data_socieconomic$iso3c)

# Adicionando dados dos estados do Brasil IBGE
ibge_state_codes <- data.frame(
  tcode = c(12, 27, 16, 13, 29, 23, 53, 32, 52,
             21, 51, 50, 31, 15, 25, 41, 26, 22,
              24, 43, 33, 11, 14, 42, 35, 28, 17),
  state = c("AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES",
              "GO", "MA", "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI",
              "RN", "RS", "RJ", "RO", "RR", "SC", "SP", "SE", "TO")
)

pib_estadual <- ipeadata("PIBPMCE") %>%
  filter(uname == "States") %>%
  mutate(year = year(date)) %>%
  select(tcode, value, year) %>%
  left_join(ibge_state_codes, by = "tcode") %>%
  arrange(year, state)

# Os dados do World Bank geralmente estão em dólares americanos (USD). 
# Para garantir comparabilidade, será necessário converter os valores de PIB
#  estadual (em reais) para dólares. Isso pode ser feito utilizando a taxa
#  de câmbio nominal anual média (R$ / USD) correspondente a cada ano.
# taxa de câmbio comercial - venda - média anual para cada ano
taxa_cambio <- ipeadata("BM_ERV") %>%
  mutate(year = year(date)) %>%
  select(year, value) %>%
  rename("taxa_cambio" = "value")

# consultar worldbank
# Para maior precisão, busque ajustar os dados para PPP se estiver disponível.
pib_estadual_dollar <- pib_estadual %>%
  left_join(taxa_cambio, by = "year") %>%
  mutate(pib_dollar_current = value/taxa_cambio) %>%
  select(state, year, pib_dollar_current) %>%
  rename(iso3c = "state")

# Agora você pode fazer o merge com o dataset socieconômico
dados_socieconomic_merge <- bind_rows(
  data_socieconomic,
  pib_estadual_dollar) %>%
  select("iso3c","year", "pop_density(person/km2)", "pib_dollar_current") 

# unidades administrativas nao encontradas
setdiff(adm_unit_list$sigla_admu,
   dados_socieconomic_merge$iso3c) 

# Save in 01_data_cleaned
save(
  dados_socieconomic_merge,
  file = file.path(
    "01_data_cleaned",
    "data_socieconomic_cleaned.RData")
)
