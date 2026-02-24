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
adjust_country_years <- function(df_country){
  
  first_year <- min(df_country$year, na.rm = TRUE)
  
  # ---------- PARTE 1: ajustar anos < 1820 ----------
  if(first_year < 1820 & 1820 %in% df_country$year){
    
    ref_1820 <- df_country %>% 
      filter(year == 1820)
    
    df_country <- df_country %>%
      mutate(
        log_gdp_final = if_else(
          year < 1820,
          ref_1820$log_gdp_final,
          log_gdp_final
        ),
        gdp_pc_final = if_else(
          year < 1820,
          ref_1820$gdp_pc_final,
          gdp_pc_final
        )
      )
  }
  
  # ---------- PARTE 2: criar 2023–2025 ----------
  if(2022 %in% df_country$year){
    
    ref_2022 <- df_country %>% 
      filter(year == 2022)
    
    future_years <- setdiff(2023:2025, df_country$year)
    
    if(length(future_years) > 0){
      
      future_rows <- tibble(
        year = future_years,
        country = ref_2022$country,
        countrycode = ref_2022$countrycode,
        name_en = ref_2022$name_en,
        macroregion = ref_2022$macroregion,
        log_gdp_final = ref_2022$log_gdp_final,
        gdp_pc_final = ref_2022$gdp_pc_final,
        gdp_source = ref_2022$gdp_source,
        gdp_uncertainty = NA_character_
      )
      
      df_country <- bind_rows(df_country, future_rows)
    }
  }
  
  # ---------- PARTE 3: marcar edge_imputed ----------
  df_country <- df_country %>%
    mutate(
      gdp_uncertainty = case_when(
        year < 1820 ~ "edge_imputed",
        year > 2022 ~ "edge_imputed",
        TRUE ~ gdp_uncertainty
      )
    ) %>%
    arrange(year)
  
  return(df_country)
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
                          "maddison_macroregion_data", 
                          "politic_data",
                          "qog_data",
                          "institutions",
                          "geographic_shape_data",
                          "list_per_admunit",
                          "data_wallacean_unnested"))); gc()

# HANDLING DATA ----
admunit <- list_per_admunit %>%
  select(name_en, ISO3) %>%
  distinct() %>%
  arrange(ISO3)

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
  by = c("ISO3" = "country_text_id")
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
  filter(ISO3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
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
  filter(!ISO3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  left_join(maddinson_data_edit, by = c("ISO3" = "countrycode")) 

# Obtain colonial origin
qog_data_edit <- qog_data %>%
  remove_missing() %>%
  distinct(ccodealp, .keep_all = TRUE) %>%
  mutate(colonial_origin = if_else(ht_colonial == 0, 0L, 1L)) %>%
  select(ccodealp, colonial_origin)

# Obtain area 
geographic_shape_data <- geographic_shape_data %>%
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
  st_transform(crs = st_crs(geographic_shape_data))

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
  #filter(ISO3 != "GEO") %>%
  left_join(geographic_shape_data %>% 
              # se nao colocar is.na(ISO3) ele remove os NA's tambem
              #filter(is.na(ISO3) | ISO3 != "GEO") %>%
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
            by = c("ISO3" = "ccodealp")) %>%
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

# AGGREGATING THE DATA OVER TIME ----
occ_year <- data_wallacean_unnested %>%
  filter(year > 1500) %>% 
  # spp so pode ser encontrada depois de ja ter sido descrita
  mutate(year_modified = if_else(year < YearOfDescription, YearOfDescription, year)) %>%
  ungroup() %>%
  distinct(name_en, year_modified) %>%
  arrange(name_en, year_modified) %>%
  rename(year = year_modified) %>%
  left_join(geographic_shape_data %>%
              select(ISO3, name_en) %>%
              st_drop_geometry(), by = "name_en")

bigsix_occ <- occ_year %>%
  filter(ISO3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  select(-name_en) %>%
  distinct() %>%
  mutate(
    name_en = recode(
      ISO3,
      "BRA" = "Brazil",
      "RUS" = "Russian Federation",
      "CHN" = "China",
      "AUS" = "Australia",
      "CAN" = "Canada",
      "USA" = "United States"
    )
  ) %>%
  select(-ISO3)

first_year <- occ_year %>%
  filter(!ISO3 %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  bind_rows(bigsix_occ) %>%
  group_by(name_en) %>%
  summarise(
    ISO3       = first(ISO3),
    first_year = min(year),
    last_year  = max(year),
    n_points   = n(),
  ) %>%
  mutate(
    ISO3 = case_when(
      name_en == "Australia"           ~ "AUS",
      name_en == "Brazil"              ~ "BRA",
      name_en == "Canada"              ~ "CAN",
      name_en == "China"               ~ "CHN",
      name_en == "Russian Federation"  ~ "RUS",
      name_en == "United States"       ~ "USA",
      TRUE ~ ISO3   # mantem o valor original se não estiver na lista
    )
  )

maddison_filtered <- maddinson_data %>%
  mutate(pop = pop * 1000) %>%
  left_join(first_year, by = c("countrycode" = "ISO3")) %>%
  filter(year >= first_year) %>%
  mutate(log_gdp = log(gdp_2011)) %>%
  select(-pop, -first_year, -last_year, -n_points)

df_interp <- maddison_filtered %>%
  # Agrupa por país para interpolar cada série separadamente
  group_by(country) %>%
  group_modify(~{
    
    # Cria uma sequência completa de anos entre o mínimo e o máximo disponível
    years_full <- tibble(
      year = seq(min(.x$year), max(.x$year), by = 1)
    )
    
    # Junta a sequência completa com os dados originais do país
    # Isso garante que anos faltantes apareçam como NA
    df_full <- years_full %>%
      left_join(.x, by = "year") %>%
      arrange(year)
    
    # Realiza interpolação log-linear:
    # approx() usa apenas os anos onde log_gdp não é NA
    interp_vals <- approx(
      x = df_full$year[!is.na(df_full$log_gdp)],  # anos observados
      y = df_full$log_gdp[!is.na(df_full$log_gdp)],  # valores observados
      xout = df_full$year  # anos completos
    )$y
    
    # Adiciona colunas interpoladas e indicador de interpolação
    df_full %>%
      mutate(
        # log do PIB interpolado
        log_gdp_interp = interp_vals,
        
        # PIB per capita interpolado (voltando do log)
        gdp_pc_interp = exp(log_gdp_interp),
        
        # Indicador: 1 = interpolado, 0 = observado
        # Baseado no valor original (antes da interpolação)
        is_interpolated = if_else(is.na(log_gdp), 1, 0)
      )
  }) %>%
  # Remove agrupamento
  ungroup() #%>%

# Plot
ggplot(df_interp %>% filter(country == "Russian Federation"), aes(x = year)) +
  geom_point(aes(y = log_gdp), alpha = 0.6, color = "blue") + # pontos observados
  # Intervalo da Primeira Guerra
  geom_rect(aes(xmin = 1914, xmax = 1918, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = 0.1, inherit.aes = FALSE) +
  # Intervalo da Segunda Guerra
  geom_rect(aes(xmin = 1939, xmax = 1945, ymin = -Inf, ymax = Inf),
            fill = "grey50", alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = log_gdp_interp), color = "red", linewidth = 0.8) + # interpolação
  theme_minimal() +
  labs(
    title = "Log-linear interpolation of GDP per capta",
    x = "Year",
    y = "log(GDP per capta)",
    color = "Legenda"
  ) 

# Macroregioes
head(maddison_macroregion_data)

macro_long <- maddison_macroregion_data |>
  pivot_longer(-year,
               names_to = "macroregion",
               values_to = "gdp_region")

macro_interp <- macro_long %>% # Agrupa por país para interpolar cada série separadamente
  mutate(log_gdp = log(gdp_region)) %>%
  group_by(macroregion) %>%
  group_modify(~{
    # Cria uma sequência completa de anos entre o mínimo e o máximo disponível
    years_full <- tibble(
      year = seq(min(.x$year), max(.x$year), by = 1)
    )
    # Junta a sequência completa com os dados originais do país
    # Isso garante que anos faltantes apareçam como NA
    df_full <- years_full %>%
      left_join(.x, by = "year") %>%
      arrange(year)
    # Realiza interpolação log-linear:
    # approx() usa apenas os anos onde log_gdp não é NA
    interp_vals <- approx(
      x = df_full$year[!is.na(df_full$log_gdp)],  # anos observados
      y = df_full$log_gdp[!is.na(df_full$log_gdp)],  # valores observados
      xout = df_full$year  # anos completos
    )$y
    # Adiciona colunas interpoladas e indicador de interpolação
    df_full %>%
      mutate(
        # log do PIB interpolado
        log_gdp_interp_macro = interp_vals,
        
        # PIB per capita interpolado (voltando do log)
        gdp_pc_interp_macro = exp(log_gdp_interp_macro),
        
        # Indicador: 1 = interpolado, 0 = observado
        # Baseado no valor original (antes da interpolação)
        is_interpolated_macro = if_else(is.na(log_gdp), 1, 0),
        is_macroregion = 1,
      )
  }) %>%
  # Remove agrupamento
  ungroup() %>%
  select(-gdp_region, -log_gdp)

# Unindo os dados em ambas as escalas
df_economy <- left_join(df_interp, macro_interp,
                        by = c("macroregion", "year")) 

# Paises da uniao sovietica que faziam parte do que era a Russia entre
# a data de 1950 ate eles comecarem a ter os proprios dados. Entao vamos
# lidar com eles separadamente
ex_ussr <- c(
  "Armenia", "Azerbaijan", "Belarus", "Georgia", 
  "Kazakhstan", "Kyrgyzstan", "Latvia", "Lithuania", 
  "Republic of Moldova", "Tajikistan", "Turkmenistan", 
  "Ukraine", "Uzbekistan"
)

# Vamos colocar os dados macroeconomicos das regioes ate 1950 para
# tornar as regioes todas comparadas. A maioria dos dados de PIB começam
# em 1950, como vamos construir ranks anuais, precisamos que as medidas
# de cada ano sejam comparaveis. Entao a comparacao justa seria
# regiao vs regiao ate 1950
# pais vs pais a partir de 1950
df_rest <- df_economy %>% filter(!country %in% ex_ussr)

df_rest_economy <- df_rest %>%
  mutate(
    # definir regime temporal primeiro
    economic_regime = if_else(year <= 1949, "macro_period", "national_period"),
    
    # definir fonte do dado respeitando o regime
    gdp_source = case_when(
      economic_regime == "macro_period"   ~ "macroregion",
      economic_regime == "national_period" & !is.na(log_gdp_interp) ~ "country",
      economic_regime == "national_period" & is.na(log_gdp_interp)  ~ "macroregion_fallback",
      TRUE ~ NA_character_
    ),
    
    # série final respeitando o regime
    log_gdp_final = case_when(
      economic_regime == "macro_period"   ~ log_gdp_interp_macro,
      economic_regime == "national_period" & !is.na(log_gdp_interp) ~ log_gdp_interp,
      TRUE ~ log_gdp_interp_macro
    ),
    
    gdp_pc_final = case_when(
      economic_regime == "macro_period"   ~ gdp_pc_interp_macro,
      economic_regime == "national_period" & !is.na(gdp_pc_interp) ~ gdp_pc_interp,
      TRUE ~ gdp_pc_interp_macro
    ),
    
    # incerteza agora fica interpretável historicamente
    gdp_uncertainty = case_when(
      economic_regime == "macro_period" ~ "regional_imputed",
      gdp_source == "country" & is_interpolated == FALSE ~ "observed",
      gdp_source == "country" & is_interpolated == TRUE  ~ "interpolated",
      gdp_source == "macroregion_fallback"               ~ "missing_country_estimated",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-gdp_2011, -log_gdp, -log_gdp_interp,
         -is_interpolated_macro, -gdp_pc_interp, -is_interpolated,
         -log_gdp_interp_macro, -gdp_pc_interp_macro, 
         -is_macroregion, -economic_regime)

# Aqui vamos preencher para paises que pertenciam a ex-uniao sovietica
# dados da federacao russa para eles ate terem seus proprios indicadores
first_year_df <- maddinson_data_edit %>%
  filter(!is.na(gdp_2011)) |>
  group_by(country) |>
  summarise(
    first_year = min(year),
    last_year  = max(year),
    n_points   = n()
  ) %>% 
  filter(first_year > 1950) %>%
  select(country, first_year) %>%
  filter(country %in% ex_ussr)

# Pega só os anos e PIB da Rússia
russia_gdp <- df_economy %>%
  filter(country == "Russian Federation") %>%
  select(year, log_gdp_interp, gdp_pc_interp)
# Para os países ex-URSS
df_ussr <- df_economy %>% 
  filter(country %in% ex_ussr) %>%
  left_join(first_year_df, by = "country") %>%
  
  # junta os dados da Rússia
  left_join(russia_gdp, by = "year", suffix = c("", "_russia")) %>%
  
  mutate(
    sovereign_economy = case_when(
      year <= 1949 ~ "macroregion",
      year >= 1950 & year < first_year ~ "Russia",
      TRUE ~ country
    ),
    
    log_gdp_final = case_when(
      sovereign_economy == "macroregion" ~ log_gdp_interp_macro,
      sovereign_economy == "Russia" ~ log_gdp_interp_russia,
      TRUE ~ log_gdp_interp
    ),
    
    gdp_pc_final = case_when(
      sovereign_economy == "macroregion" ~ gdp_pc_interp_macro,
      sovereign_economy == "Russia" ~ gdp_pc_interp_russia,
      TRUE ~ gdp_pc_interp
    ),
    
    gdp_source = sovereign_economy,
    
    gdp_uncertainty = case_when(
      sovereign_economy == "macroregion" ~ "regional_imputed",
      sovereign_economy == "Russia" ~ "sovereign_imputed",
      sovereign_economy == country & is_interpolated == 1 ~ "interpolated",  # respeita interpolação original
      sovereign_economy == country & is_interpolated == 0 ~ "observed"
    )
  ) %>%
  select(-first_year, -log_gdp_interp_russia, -gdp_pc_interp_russia,
         -gdp_2011, -log_gdp, -log_gdp_interp,
         -gdp_pc_interp, -log_gdp_interp_macro, -gdp_pc_interp_macro,
         -sovereign_economy, -is_interpolated, -is_interpolated_macro, 
         -is_macroregion)

names(df_ussr)
names(df_rest_economy)

df_economy_all <- bind_rows(
  df_rest_economy,
  df_ussr
) %>%
  group_by(country) %>%
  group_modify(~ adjust_country_years(.x)) %>%
  ungroup() %>%
  arrange(country, year)

View(df_economy_all)

# edge_imputed: anos que nao tem dados e usamos os ultimos dados presentes 1820 < & 2022 >
# interpolated sao dados que foram interpolados entre anos de dados observados ausentes
# missing_country_estimated occ aparecia, mas o pais ainda nao tinha dados observados na epoca, usamos dados da regiao
# observed dados observados no maddison
# regional_imputed assumimos os dados regionais entre 1820 ate 1950
# sovereign_imputed pais pertencia a URSS, entao usamos dados da Russia para eles

# Vamos separar o que eh unidade administrativa do que eh paises, pois existem
# tratamentos diferentes entre eles
df_economy_bigsix <- df_economy_all %>%
  filter(countrycode %in% c("BRA", "USA", "RUS", "CHN", "AUS", "CAN"))

df_economy_countries <- df_economy_all %>%
  filter(!countrycode %in% c("BRA", "USA", "RUS", "CHN", "AUS", "CAN"))

# Region data
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
  mutate(
    GID_0 = case_when(
      region == "Moscow" ~ "RUS",
      TRUE ~ GID_0
    )
  )

# Ajustes para as und. adm nos paises. 

# Brazil
# Os estados brasileiros passam a ter dados a partir de 1970
# usar o valor de Goiás de 1970 ate 1984 para Tocantins
# usar o valor de Mato Grosso de 1970 ate 1974 para Mato Grosso do Sul
brazil_data <- region_data %>% 
  filter(GID_0 == "BRA" & year >= 1970)
# 1. Valores de referência
goias_ref <- brazil_data %>%
  filter(region == "Goiás", year >= 1970, year <= 1988) 
mt_ref <- brazil_data %>%
  filter(region == "Mato Grosso", year >= 1970, year <= 1974) 
# 2. Criar linhas artificiais para Tocantins e MS
tocantins_fill <- goias_ref %>%
  mutate(region = "Tocantins")
ms_fill <- mt_ref %>%
  mutate(region = "Mato Grosso Do Sul")
# 3. Unir tudo ao dataset original
brazil_fixed <- brazil_data %>%
  bind_rows(tocantins_fill, ms_fill) %>%
  arrange(region, year)

# Canada
# Antes de 1º de abril de 1999, o que hoje conhecemos como Nunavut não existia 
# como uma entidade administrativa separada. Ele fazia parte dos 
# Northwest Territories (Territórios do Noroeste). No DOSE chama 
# "Northwest Territories Including Nunavut" entao vamos usar os valores dessas
# series para os estados ate eles terem dados proprios
canada_edit_data <- region_data %>% 
  filter(region == "Northwest Territories Including Nunavut" & year > 1983)
nunavut_ref <- canada_edit_data %>%
  mutate(region = "Nunavut") 
northwest_ref <- canada_edit_data %>%
  mutate(region = "Northwest Territories") 
# Unir ao dataset original tirando o territorio 'Northwest Territories Including Nunavut'
canada_fixed <- region_data %>%
  filter(GID_0 == "CAN" & region != "Northwest Territories Including Nunavut") %>%
  bind_rows(nunavut_ref, northwest_ref) %>%
  arrange(region, year)

# China
# 1988 Hainan foi oficialmente separada de Guangdong e elevada à categoria 
# de província independente e, simultaneamente
china_edit_data <- region_data %>% 
  filter(region == "Guangdong" & year < 1978)
hainan_ref <- china_edit_data %>%
  mutate(region = "Hainan") 
china_fixed <- region_data %>%
  filter(GID_0 == "CHN") %>%
  bind_rows(hainan_ref) %>%
  arrange(region, year)

# Russia
# a Ingushetia só existe como unidade separada a partir de 1994, e a Chechênia 
# a partir de 2005, antes disso o PIB das duas era computado junto em 
# uma única unidade administrativa.
russia_edit_data <- region_data %>% 
  filter(region == "Ingushetia" & year < 2005)
chechen_ref <- russia_edit_data %>%
  mutate(region = "Chechen") 
russia_fixed <- region_data %>%
  filter(GID_0 == "RUS") %>%
  bind_rows(chechen_ref) %>%
  arrange(region, year)

# Australia comeca a ter dados a partir de 1970 para PIB
australia_fixed <- region_data %>%
  filter(GID_0 == "AUS" & year > 1953)

# Unindo os dados para as und adm.
# AUS e EUA estao todos partindo do mesmo ano, nao precisam de ajustes
region_data_edit <- region_data %>%
  filter(!GID_0 %in% c("BRA", "CAN", "CHN", "RUS", "AUS")) %>%
  bind_rows(brazil_fixed, canada_fixed, china_fixed, russia_fixed, australia_fixed)

region_data_edit %>%
  group_by(region) %>%
  summarise(
    region       = first(region),
    GID_0 = first(GID_0),
    first_year = min(year),
    last_year  = max(year),
    n_points   = n(),
  ) %>%
  arrange(GID_0, region) %>% head()

region_data_all <- region_data_edit %>%
  mutate(log_gdp = log(gdp_2015)) %>%
  group_by(region) %>%
  group_modify(~{
    
    years_full <- tibble(year = seq(min(.x$year), max(.x$year), 1))
    
    df_full <- years_full %>%
      left_join(.x, by = "year") %>%
      arrange(year)
    
    interp_vals <- approx(
      x = df_full$year[!is.na(df_full$log_gdp)],
      y = df_full$log_gdp[!is.na(df_full$log_gdp)],
      xout = df_full$year
    )$y
    
    df_full %>%
      mutate(
        log_gdp_final = interp_vals,
        gdp_pc_final = exp(log_gdp_final),
        is_interpolated = if_else(is.na(log_gdp), 1, 0),
        # Copia metadados do grupo
        countrycode    = .x$countrycode[1],
        GID_0          = .x$GID_0[1],
        #region_tolower = .x$region_tolower[1]
      )
  }) %>%
  ungroup() %>% 
  select(-country, -GID_1, -pop_abs,
         -gdp_2015, -region_tolower, -log_gdp)

# Unir os dois agora
states <- region_data_all %>%
  select(GID_0, region) %>%
  distinct()

min_year <- region_data_all %>%
  group_by(GID_0) %>%
  summarise(
    min_year = first(year)
  )

df_economy_bigsix_all <- left_join(
  states,
  df_economy_bigsix,
  by = c("GID_0" = "countrycode")
) %>%
  left_join(min_year, by = "GID_0") %>%
  filter(year < min_year) %>% 
  select(-min_year)

# preechendo a macroregion
ref_meta <- df_economy_bigsix_all %>%
  select(GID_0, macroregion) %>% 
  distinct()

# Criando a serie do ultimo ano de cada und adm ate 2025
last_rows_extended <- region_data_all %>%
  group_by(GID_0) %>%
  filter(year == max(year)) %>%      # pega a última linha de cada país
  ungroup() %>%
  rowwise() %>%
  do({
    last_row <- .
    tibble(
      GID_0 = last_row$GID_0,
      region = last_row$region,
      year = seq(last_row$year + 1, 2025, 1),
      gdp_pc_final = last_row$gdp_pc_final,
      log_gdp_final = last_row$log_gdp_final,
      gdp_uncertainty = "edge_imputed",
      gdp_source = "adm_unit"
    )
  }) %>%
  ungroup() %>%
  left_join(ref_meta, by = "GID_0")

df_economy_admunit <- bind_rows(
  last_rows_extended,
  df_economy_bigsix_all,
  region_data_all %>%
    left_join(ref_meta, by = "GID_0") %>% 
    mutate(
      gdp_source = "adm_unit",
      gdp_uncertainty = case_when(
        is_interpolated == 1 ~ "interpolated",
        is_interpolated == 0 ~ "observed",
        TRUE ~ NA_character_)
)) %>%
  arrange(region, year) %>%
  select(-is_interpolated, -country, -name_en) %>%
  rename(
    countrycode = GID_0,
    name_en = region
  )

names(df_economy_admunit)
names(df_economy_countries)

data_socieconomic_temporal <- bind_rows(
  df_economy_admunit,
  df_economy_countries
) %>% 
  mutate(
    country = case_when(
      countrycode == "AUS" ~ "Australia", 
      countrycode == "BRA" ~ "Brazil",
      countrycode == "CAN" ~ "Canada",
      countrycode == "CHN" ~ "China",
      countrycode == "RUS" ~ "Russian Federation",
      countrycode == "USA" ~ "United States",
      TRUE ~ country))

data_socieconomic_temporal %>%
  group_by(name_en) %>%
  summarise(
    countrycode       = first(countrycode),
    first_year = min(year),
    last_year  = max(year),
    n_points   = n(),
  ) %>% View()

# Precisamos harmonizar os nomes das unidades administrativas nas bases
# de dados socieconomicos e na nossa base de ocorrencia
setdiff(occ_year$name_en,
        data_socieconomic_temporal$name_en)

setdiff(data_socieconomic_temporal$name_en,
        occ_year$name_en)

harmonize_names <- c(
  # Rússia - abreviações / traduções
  "Adygea" = "Republic of Adygea",
  "Arhangelsk" = "Arkhangelsk",
  "Buryatia" = "Republic of Buryatia",
  "Chechen" = "Chechen Republic",
  "Chukotka" = "Chukotka Autonomous Okrug",
  "Chuvash" = "Chuvash Republic",
  "Dagestan" = "Republic of Dagestan",
  "Ingushetia" = "Republic of Ingushetia",
  "Jewish Autonomous Region" = "Jewish",
  "Kabardino-Balkar" = "Kabardino-Balkaria",
  "Kalmykia" = "Republic of Kalmykia",
  "Kamchatka" = "Kamchatka Krai",
  "Karachay-Cherkess" = "Karachay-Cherkess Republic",
  "Khabarovsk" = "Khabarovsk Krai",
  "Khakassia" = "Republic of Khakassia",
  "Komi" = "Komi Republic",
  "Krasnodar" = "Krasnodar Krai",
  "Krasnoyarsk" = "Krasnoyarsk Krai",
  "Mari El" = "Mari El Republic",
  "Mordovia" = "Republic of Mordovia",
  "North Ossetia-Alania" = "Republic of North Ossetia-Alania",
  "Primorsky" = "Primorsky Krai",
  "Sakha" = "Sakha Republic",
  "Stavropol" = "Stavropol Krai",
  "Tatarstan" = "Republic of Tatarstan",
  "Transbaikal" = "Zabaykalsky Krai",
  "Tyva" = "Tuva Republic",
  "Udmurt" = "Udmurt Republic",
  
  # Brasil - capitalização
  "Mato Grosso Do Sul" = "Mato Grosso do Sul",
  "Rio De Janeiro" = "Rio de Janeiro",
  "Rio Grande Do Norte" = "Rio Grande do Norte",
  "Rio Grande Do Sul" = "Rio Grande do Sul",
  "Distrito Federal" = "Federal",
  
  # Canadá - capitalização
  "Newfoundland And Labrador" = "Newfoundland and Labrador",
  
  # EUA - nome alternativo
  "District of Columbia" = "Washington D.C."
)

data_socieconomic_temporal <- data_socieconomic_temporal %>%
  mutate(
    name_en = recode(name_en, !!!harmonize_names)
  ) %>%
  group_by(year) %>%
  mutate(gdp_rank = rank(log_gdp_final) / n(),
         gdp_percentile = percent_rank(log_gdp_final))

setdiff(occ_year$name_en,
        data_socieconomic_temporal$name_en)
setdiff(teste$name_en,
        data_socieconomic_temporal$name_en)

data <- left_join(
  occ_year,
  data_socieconomic_temporal,
  by = c("name_en", "year")) 

countries_without_gdp <- data %>%
  filter(is.na(gdp_pc_final)) %>% 
  distinct(name_en) 

# Save temporal series
save(
  data_socieconomic_temporal,
  countries_without_gdp,
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_temporal.RData")
)
