library(sf)
library(tidyverse)
library(stringi)

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

# explorar dados socieconomicos
# https://www.nature.com/articles/s41597-023-02323-8
load(file.path("00_raw_data",
               "geographic_shape_data.RData"))

data_dose <- read.csv2("00_raw_data/raw_data/DOSE_V2.11.csv", sep = ',') %>%
  select(country, GID_0, region, GID_1, year, pop, grp_pc_usd_2015) 

# Conta a frequência de cada adm0_a3
counts <- geographic_shape_data %>%
  count(adm0_a3)

# Dataset 1: apenas os adm0_a3 que aparecem uma vez
shp_country <- geographic_shape_data %>%
  filter(adm0_a3 %in% counts$adm0_a3[counts$n == 1]) %>%
  st_transform(crs = 6933) %>%   # reprojeta para Equal Area
  mutate(area_km2 = round(as.numeric(st_area(geometry)) / 10^6),2)  # converte m² → km²

country_economic <- data_dose %>%
  select(GID_0, pop, grp_pc_usd_2015) %>%
  group_by(GID_0) %>%
  mutate(
    pop = as.numeric(pop),
    grp_pc_usd_2015 = as.numeric(grp_pc_usd_2015)
  ) %>%
  slice_tail(n = 10) %>%
  summarise(
    pop = mean(pop, na.rm = TRUE),
    grp_pc_usd_2015 = mean(grp_pc_usd_2015, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(everything(), ~replace(., is.nan(.), NA))) 

#  summarise(
#    mean_pop = mean(as.numeric(pop), na.rm = TRUE),
#    mean_grp = mean(as.numeric(grp_pc_usd_2015), na.rm = TRUE)
#  ) %>%
#  mutate(across(everything(), ~replace(., is.nan(.), NA))) 

socieconomic_country <- left_join(shp_country,
                                  country_economic,
                                  by = c("adm0_a3" = "GID_0")) %>%
  st_drop_geometry()
visdat::vis_miss(socieconomic_country)
# Dataset 2: adm0_a3 que aparecem duas ou mais vezes
shp_adm_unit <- geographic_shape_data %>%
  filter(adm0_a3 %in% counts$adm0_a3[counts$n >= 2]) %>%
  mutate(
    name_en = stri_trans_general(tolower(name_en), "Latin-ASCII"))
adm_unit <- unique(shp_adm_unit$adm0_a3)

region_economic <- data_dose %>%
  filter(GID_0 %in% adm_unit) %>%
  select(region, pop, grp_pc_usd_2015) %>%
  mutate(
    region = stri_trans_general(tolower(region), "Latin-ASCII"),
    pop = as.numeric(pop),
    grp_pc_usd_2015 = as.numeric(grp_pc_usd_2015)
  ) %>%
  group_by(region) %>%
  slice_tail(n = 10) %>%
  summarise(
    pop = mean(pop, na.rm = TRUE),
    grp_pc_usd_2015 = mean(grp_pc_usd_2015, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(everything(), ~replace(., is.nan(.), NA)))


socieconomic_region <- left_join(shp_adm_unit,
                                 region_economic,
                                 by = c("name_en" = "region")) %>%
  st_transform(crs = 6933) %>%   # reprojeta para Equal Area
  mutate(area_km2 = round(as.numeric(st_area(geometry)) / 10^6),2) %>%  # converte m² → km²
  st_drop_geometry()

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

visdat::vis_miss(centroid_undadm)

socieconomic_data_raw <- bind_rows(socieconomic_region,
                               socieconomic_country) %>%
  mutate(
    name_en = stri_trans_general(tolower(name_en), "Latin-ASCII"))

socieconomic_data <- left_join(socieconomic_data_raw,
          centroid_undadm,
          by = "name_en") 

visdat::vis_miss(socieconomic_data)

save(
  socieconomic_data,
  file = file.path(
    "01_data_cleaned",
    "data_socieconomic_cleaned.RData")
)
