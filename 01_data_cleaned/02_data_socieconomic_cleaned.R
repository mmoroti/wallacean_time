library(sf)
library(tidyverse)
library(units)
library(countrycode)

rm(list=ls()); gc() # clean local enviroment
# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time")

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

# LOAD DATA ----
# Geographic data
load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData")
  )

# Socieconomic data
load(file.path(
  local_directory,
  "00_raw_data",
  "data_socieconomic.RData")
)

# Occurences cleaned
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "dataset_occurences.RData")
)

rm(list = setdiff(ls(), c("local_directory",
                          "dose_data",
                          "maddinson_data",
                          "politic_data",
                          "qog_data",
                          "institutions",
                          "geographic_shape_data",
                          "list_per_admunit",
                          "data_wallacean_unnested"))); gc()

# HANDLING DATA ----
admunit <- list_per_admunit %>%
  select(name_en, adm0_a3) %>%
  distinct() %>%
  arrange(adm0_a3)

politic_data_summary <- politic_data %>%
  group_by(country_text_id) %>%
  # se quiser pegar os 10 ultimos anos
  slice_tail(n = 10) %>%
  summarise(
    min_year_politic = min(year, na.rm = TRUE),
    max_year_politic = max(year, na.rm = TRUE),
    democracy_index_mean = mean(v2x_polyarchy, na.rm = TRUE),
    democracy_index_sd = sd(v2x_polyarchy, na.rm = TRUE),
    academic_freedom_mean = mean(v2clacfree, na.rm = TRUE),
    academic_freedom_sd = sd(v2clacfree, na.rm = TRUE),
    #educational_inequality = mean(e_peedgini, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    country_text_id = recode(
      country_text_id,
      "XKX" = "KOS",
      "PSB" = "PSX",
      "SSD" = "SDS",
      "SML" = "SOL"
    )
  )

admunit_politic <- left_join(
  admunit,
  politic_data_summary,
  by = c("adm0_a3" = "country_text_id")
) 

# Primeiro vamos tentar unir uma media para as unidades admnistrativas
# Dentro do DOSE, temos dados de Unidades administrativas para diversos paises
# mas so quebramos as unidades administrativas para 6 paises. Entao vamos
# separar em dois datasets, um com as regioes e outro so com os paises inteiros

# ajuste dos nomes entre as bases
map_names <- c(
  "northwest territories including nunavut" = "nunavut",
  "adygea" = "republic of adygea",
  "buryatia" = "republic of buryatia",
  "chechen" = "chechen republic",
  "chukotka" = "chukotka autonomous okrug",
  "chuvash" = "chuvash republic",
  "dagestan" = "republic of dagestan",
  "ingushetia" = "republic of ingushetia",
  "jewish autonomous region" = "jewish",
  "kabardino-balkar" = "kabardino-balkaria",
  "kalmykia" = "republic of kalmykia",
  "karachay-cherkess" = "karachay-cherkess republic",
  "khabarovsk" = "khabarovsk krai",
  "khakassia" = "republic of khakassia",
  "komi" = "komi republic",
  "krasnodar" = "krasnodar krai",
  "krasnoyarsk" = "krasnoyarsk krai",
  "mari el" = "mari el republic",
  "mordovia" = "republic of mordovia",
  "moscow city & moscow region" = "moscow",
  "north ossetia-alania" = "republic of north ossetia-alania",
  "primorsky" = "primorsky krai",
  "sakha" = "sakha republic",
  "stavropol" = "stavropol krai",
  "tatarstan" = "republic of tatarstan",
  "udmurt" = "udmurt republic",
  "kamchatka" = "kamchatka krai",
  "altai krai" = "altai republic",
  "distrito federal" = "federal"
)

# Essas duas unidades administrativas sao tratadas como uma só em nossos dados
# regiao de Moscow, por isso vou unir region e city em uma só. Para isso o correto
# é somar 'moscow region' + 'moscow city' = moscow
df_moscow <- dose_data %>%
  filter(region %in% c("Moscow city", "Moscow region")) %>%
  group_by(year) %>%
  summarise(
    region = "Moscow",
    gdp_2015 = sum(gdp_2015, na.rm = TRUE),
    pop_abs = sum(pop_abs, na.rm = TRUE),
    country = "Russia"
  )

# Regional data
region_data <- dose_data %>%
  # removendo dados de moscow desagregados
  filter(!region %in% c("Moscow city", "Moscow region")) %>%
  # voltando dados de moscow agregados
  bind_rows(df_moscow) %>%
  filter(country %in% c("Brazil", "Russia", "China", 
                        "Australia", "USA", "Canada")) %>%
  distinct(region, year, .keep_all = TRUE) %>% 
  mutate(region_tolower = tolower(region)) %>%
  mutate(region_tolower = str_replace_all(region_tolower, map_names)) %>%
  group_by(region_tolower) %>%
  slice_tail(n = 10) %>%
  summarise(
    min_year_dose = min(year, na.rm = TRUE),
    max_year_dose = max(year, na.rm = TRUE),
    gdp_2015_mean = mean(gdp_2015, na.rm = TRUE),
    gdp_2015_sd = sd(gdp_2015, na.rm = TRUE),
    pop_abs_mean = mean(pop_abs, na.rm = TRUE),
    pop_abs_sd = sd(pop_abs, na.rm = TRUE),
  ) %>%
  ungroup() 

admunit_sociopolitic <- admunit_politic %>%
  filter(adm0_a3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  mutate(region_tolower = tolower(name_en)) %>%
  left_join(region_data, by = "region_tolower") #%>%
  # verificar as unidades faltantes
  #select(name_en, adm0_a3, region_tolower, max_year_dose) %>%
  #filter(is.na(max_year_dose)) %>%
  #View()
  
# Country data
#country_data <- dose_data %>%
#  filter(!country %in% c("Brazil", "Russia", "China", 
#                        "Australia", "USA", "Canada")) %>%
#  group_by(country, GID_0, year) %>%
#  summarise(
#    pop_abs = sum(pop_abs, na.rm = TRUE),
#    gdp_2015 = sum(gdp_2015, na.rm = TRUE),
#    .groups = "drop"
#  ) %>%
#  group_by(GID_0) %>%
#  arrange(year, .by_group = TRUE) %>%
#  slice_tail(n = 10) %>%
#  summarise(
#    min_year_dose = min(year, na.rm = TRUE),
#    max_year_dose = max(year, na.rm = TRUE),
#    gdp_2015_mean = mean(gdp_2015, na.rm = TRUE),
#    gdp_2015_sd = sd(gdp_2015, na.rm = TRUE),
#    pop_abs_mean = mean(pop_abs, na.rm = TRUE),
#    pop_abs_sd = sd(pop_abs, na.rm = TRUE),
#  ) %>%
#  ungroup()
#
#country_sociopolitic <- admunit_politic %>%
#  filter(!adm0_a3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
#  left_join(country_data, by = c("adm0_a3" = "GID_0")) #%>%
#  #filter(is.na(max_year_dose)) %>%
#  #View()

# o DOSE tem cobertura esparsa para países da África e do Oriente Médio, 
# onde o Maddison (que usa dados do Banco Mundial e ONU para o período recente) 
# será essencial para completar a sua média mundial
maddinson_data_edit <- maddinson_data %>%
  filter(!countrycode %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  mutate(pop = pop * 1000) %>%
  group_by(countrycode) %>%
  # se quiser pegar os 10 ultimos anos
  slice_tail(n = 10) %>%
  summarise(
    min_year_maddison = min(year, na.rm = TRUE),
    max_year_maddison = max(year, na.rm = TRUE),
    gdp_2011_mean = mean(gdp_2011, na.rm = TRUE),
    gdp_2011_sd = sd(gdp_2011, na.rm = TRUE),
    pop_mean = mean(pop, na.rm = TRUE),
    pop_sd = sd(pop, na.rm = TRUE),
  ) 

country_sociopolitic <- admunit_politic %>%
  filter(!adm0_a3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  left_join(maddinson_data_edit, by = c("adm0_a3" = "countrycode")) 

# Obtain colonial origin
qog_data_edit <- qog_data %>%
  remove_missing() %>%
  distinct(ccodealp, .keep_all = TRUE) %>%
  mutate(colonial_origin = if_else(ht_colonial == 0, 0L, 1L)) %>%
  select(ccodealp, colonial_origin)

# Obtain area 
geographic_shape_data <- geographic_shape_data %>%
  st_transform(crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  mutate(
    area_m2 = st_area(geometry),  # Área em metros quadrados
    area_km2 = set_units(area_m2, km^2),  # Converter para km²
    area_km2 = as.numeric(area_km2)  
  )

# Obtain institutions and institutions of each adm. unit
# Aqui vou pegar os paises grandes e aqueles sem classificacao de country para
# determinar onde eles caem no nosso poligono, com isso consigo identificar
# aonde nas unidades administrativas eles caem e classificar os nao classificados
institutions_sf <- institutions %>%
  filter(is.na(country) | country %in% c("BRA", "CHN", "RUS", "USA", "CAN", "AUS")) %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Contagem das instituicoes dentro dos name_en
institutions_with_location <- institutions_sf %>%
  st_join(geographic_shape_data %>% select(name_en)) %>%
  filter(!is.na(name_en)) %>%
  st_drop_geometry() %>%
  group_by(name_en) %>%
  summarise(n_institutions = n())

# Agora posso agrupar de acordo com os paises, mas antes, quero identificar
# pelo codigo ISO3 o nome exato em nossos dados (name_en)
nameen_list <- geographic_shape_data %>%
  select(name_en, ISO3) %>%
  filter(!ISO3 %in% c("BRA", "CHN", "RUS", "USA", "CAN", "AUS")) %>%
  st_drop_geometry() %>%
  mutate(
    ISO3 = countrycode(
      name_en,
      origin = "country.name",
      destination = "iso3c",
      nomatch = NA  # Mantém NA para não-encontrados
    )
  )

df_institutions <- institutions %>%
  filter(!country %in% c("BRA", "CHN", "RUS", "USA", "CAN", "AUS")  & !is.na(country)) %>%
  group_by(country) %>%
  summarise(n_institutions = n()) %>%
  left_join(nameen_list, by = c("country" = "ISO3")) %>%
  bind_rows(institutions_with_location) %>%
  filter(country != "GEO" | is.na(country)) %>%  # Remove linhas com country == "GEO"
  group_by(name_en) %>%                          # Agrupa por name_en
  summarise(
    n_institutions = sum(n_institutions, na.rm = TRUE),  # Soma os valores
    country = first(na.omit(country)),          # Pega o primeiro país não-NA
    GEO = any(country == "GEO", na.rm = TRUE)   # Marca se tinha GEO antes
  ) %>%
  ungroup() %>%
  select(name_en, n_institutions)

# Obtain 
data_sociopolitic <- bind_rows(
  country_sociopolitic,
  admunit_sociopolitic
) %>%
  select(-region_tolower) %>%
  # Georgia país aparece duplicado mas nao temos nenhum grupo com occ la
  filter(adm0_a3 != "GEO") %>%
  left_join(geographic_shape_data %>% 
              # se nao colocar is.na(ISO3) ele remove os NA's tambem
              filter(is.na(ISO3) | ISO3 != "GEO") %>%
              select(name_en, area_km2, GlobalNorth) %>%
              st_drop_geometry(),
            by = "name_en") %>%
  distinct() %>%
  mutate(
    # GDP consolidado
    gdp_mean = coalesce(gdp_2011_mean, gdp_2015_mean),
    # População consolidada
    pop = coalesce(pop_mean, pop_abs_mean),
    Origin = case_when(
      !is.na(gdp_2011_mean) ~ "Maddison",
      !is.na(gdp_2015_mean) ~ "DOSE"
    )) %>%
  left_join(qog_data_edit,
            by = c("adm0_a3" = "ccodealp")) %>%
  left_join(df_institutions, by = "name_en")  %>%
  mutate(
    n_institutions = ifelse(is.na(n_institutions), 0, n_institutions)
  )

save(
  data_sociopolitic,
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_agregatte.RData")
)

# Agregando temporalmente os dados ----
# Obtain centroid each adm unit
centroid_undadm <- geographic_shape_data %>% 
  mutate(geometry = st_centroid(geometry)) %>%  # substitui geometria pelo centróide
  mutate(
    longitude = st_coordinates(geometry)[,1],
    latitude  = st_coordinates(geometry)[,2]
  ) %>%
  mutate(
    name_en = stri_trans_general(tolower(name_en), "Latin-ASCII")) %>%
  st_drop_geometry() %>%
  select(-adm0_a3) %>%
  distinct(name_en, .keep_all = TRUE)
