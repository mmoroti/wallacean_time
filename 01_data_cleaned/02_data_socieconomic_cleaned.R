library(sf)
library(tidyverse)

# preencher os anos intermediarios com base na variacao entre decadas
fill_population_data <- function(df) {
  # Ordena por UF e ano
  df <- df %>%
    arrange(UF, Year)
  
  # Criação de um novo dataframe para os anos preenchidos
  df_filled <- df %>%
    group_by(UF) %>%
    # Para cada par de anos consecutivos com dados, cria os anos intermediários
    do({
      # Extrair o subdataframe para um grupo de UF específico
      df_group <- .
      
      # Encontrar os pares de anos onde temos valores conhecidos
      completed_years <- df_group %>% filter(!is.na(`pop_density(person/km2)`))
      
      # Lista para armazenar os dados preenchidos
      filled_rows <- list()
      
      # Preencher para cada par de anos consecutivos com valores conhecidos
      for(i in 1:(nrow(completed_years) - 1)) {
        year_start <- completed_years$Year[i]
        year_end <- completed_years$Year[i + 1]
        pop_start <- completed_years$`pop_density(person/km2)`[i]
        pop_end <- completed_years$`pop_density(person/km2)`[i + 1]
        
        # Calcular a variação de densidade por ano
        diff_per_year <- (pop_end - pop_start) / (year_end - year_start)
        
        # Criar os anos intermediários
        for(year in (year_start + 1):(year_end - 1)) {
          pop_value <- pop_start + (year - year_start) * diff_per_year
          filled_rows <- append(filled_rows, list(tibble(UF = completed_years$UF[i], 
                                                          Year = year, 
                                                          `pop_density(person/km2)` = pop_value)))
        }
      }
      
      # Retornar as linhas preenchidas, além das linhas originais
      bind_rows(df_group, bind_rows(filled_rows)) %>%
        arrange(Year)  # Ordenar as linhas por ano
    })
  
  # Retornar o dataframe preenchido
  return(df_filled)
}

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
  select(adm_unit, sigla_admu) %>%
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

pib_estadual_dollar <- pib_estadual %>%
  left_join(taxa_cambio, by = "year") %>%
  mutate(pib_dollar_current = value/taxa_cambio) %>%
  select(state, year, pib_dollar_current) %>%
  rename(iso3c = "state")

population_brazil <- brazil_population %>%
  mutate(across(`1872`:`1991`, ~ as.numeric(gsub(",", ".", .)))) %>%
  pivot_longer(cols = -UF, names_to = "Year", values_to = "pop_density(person/km2)") %>%
  mutate(Year = as.numeric(Year))

# Para preencher os dados de densidade populacional entre décadas,
# foi adotada uma abordagem que distribui proporcionalmente a variação
# observada entre os anos de início e fim de cada década, considerando
# o aumento ou declínio da população. A diferença de densidade 
# populacional entre os anos conhecidos foi calculada e dividida 
# igualmente entre os anos intermediários, refletindo as tendências 
# de crescimento ou diminuição populacional ao longo do tempo. 
# Essa metodologia permite capturar as variações populacionais de forma
# contínua, assumindo uma tendência linear entre as décadas.
# Aplicando a função para preencher os dados
population_brazil_filled <- fill_population_data(population_brazil) %>%
  filter(Year >= 1900)

# State Brazil
uf_table <- tibble::tibble(
  UF = c("Acre", "Alagoas", "Amapa", "Amazonas", "Bahia", "Ceara", 
             "Distrito Federal", "Espirito Santo", "Goias", "Maranhao", 
             "Mato Grosso", "Mato Grosso do Sul", "Minas Gerais", "Para", 
             "Paraiba", "Parana", "Pernambuco", "Piaui", "Rio de Janeiro", 
             "Rio Grande do Norte", "Rio Grande do Sul", "Rondonia", 
             "Roraima", "Santa Catarina", "Sao Paulo", "Sergipe", "Tocantins"),
  iso3c = c("AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES", "GO", "MA", 
            "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI", "RJ", "RN", 
            "RS", "RO", "RR", "SC", "SP", "SE", "TO")
)

# Converte para siglas usando dplyr
population_brazil_filled_admunit <- left_join(
  population_brazil_filled, 
  uf_table, by = "UF"
) %>%
  ungroup() %>%
  select(-UF) %>%
  rename(year = Year) %>%
  remove_missing()

# join pop_density and pib_dollar_current
brazil_data <- full_join(population_brazil_filled_admunit,
  pib_estadual_dollar,
  by = c("iso3c", "year"))

# Agora você pode fazer o merge com o dataset socieconômico
dados_socieconomic_merge <- bind_rows(
  data_socieconomic,
  brazil_data) %>%
  select("iso3c","year", "pop_density(person/km2)", "pib_dollar_current") %>%
  rename(sigla_admu = "iso3c")

# unidades administrativas nao encontradas
setdiff(adm_unit_list$sigla_admu,
   dados_socieconomic_merge$sigla_admu) 

# Save in 01_data_cleaned
save(
  dados_socieconomic_merge,
  file = file.path(
    "01_data_cleaned",
    "data_socieconomic_cleaned.RData")
)
