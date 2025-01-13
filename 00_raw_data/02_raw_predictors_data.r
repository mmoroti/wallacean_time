library(WDI) # API world bank
library(ipeadatar) # API ipea (Brazilian states data)
library(tidyverse) # data handling

options(scipen = 999, digits = 15) # evitar anotacoes cientificas 

head(data.frame(WDIsearch("GDP")))
WDIsearch("education")  # Busca por indicadores relacionados à educação
WDIsearch("health")     # Busca por indicadores relacionados à saúde

indicadores_economicos <- c(
  "NY.GDP.MKTP.KD.ZG",  # Crescimento do PIB (% anual)
  "NY.GDP.PCAP.KD.ZG",  # Crescimento do PIB per capita (% anual)
  "EN.POP.DNST",        # Densidade populacional (pessoas por km² de terra)
  "NY.GDP.MKTP.CD",     # PIB em dólares correntes
  "SP.POP.TOTL",        # População total
  "GB.XPD.RSDV.GD.ZS"   # Despesa em P&D como % do PIB (indicador similar para "Human resources in R&D")
)

# Baixar dados de PIB e população
data_worldbank <- WDI(
  indicator = indicadores_economicos,
    start = 1960, end = 2024
)

# rename_cols
data_socieconomic <- data_worldbank %>%
  rename(
    "pib_growth_anual(%)" = "NY.GDP.MKTP.KD.ZG",
    "pib_growth_percapta_anual(%)" = "NY.GDP.PCAP.KD.ZG",
    "pop_density(person/km2)" = "EN.POP.DNST",
    "pib_dollar_current" = "NY.GDP.MKTP.CD",
    "total_population" = "SP.POP.TOTL",
    "invest_pd_pib" = "GB.XPD.RSDV.GD.ZS",
  ) %>%
  select(-iso2c)
#names(data_socieconomic)

# IPEA data (States of Brazilian) ----
# podem ser acessadas em
# http://www.ipeadata.gov.br/ExibeSerieR.aspx?stub=1&serid=1540855420&MINDATA=2014&MAXDATA=2030&TNIVID=2&TPAID=1&module=R
search_series("Taxa de câmbio")
View(available_series())

# PIB estadual
#print(metadata("PIBPMCE")$comment
# O produto interno bruto (PIB) é o total dos bens e serviços
# produzidos pelas unidades produtoras residentes destinados
# ao consumo final, sendo equivalente à soma dos valores adicionados
# pelas diversas atividades econômicas acrescida dos impostos,
# líquidos de subsídios, sobre produtos. Equivale também à soma
# dos consumos finais de bens e serviços valorados a preço de mercado,
# sendo também equivalente à soma das rendas primárias.
# Neste caso é a medição do PIB de um determinado
# Exemplo de extração dos dados
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
  select(year, value)

# consultar worldbank
# Para maior precisão, busque ajustar os dados para PPP se estiver disponível.

# Save socieconomic dataset
save(
  data_socieconomic,
  file = file.path(
    "00_raw_data",
    "data_socieconomic.RData"
  )
)
