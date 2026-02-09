library(sf)
library(tidyverse)
library(rgbif)
library(rnaturalearth)

## set your directory in this block
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# CHANGE GRID-CELLS TO POLYGONS---- 
load(file.path(
      local_directory,
       "00_raw_data",
       "tetrapodstraits_data.RData"))

vert_assemblages <- data.table::fread(
  file.path(
    local_directory,
    "Shapefiles",
    "TetrapodTraits_Community",
    "Tetrapod_360.csv"),
  stringsAsFactors=T
)

grid_cells_sf <- sf::st_read(
  dsn = file.path(
  local_directory,
  "Shapefiles",
  "TetrapodTraits_Community"),
    layer='gridcells_110km'
)

vert_assemblages_sf <- vert_assemblages %>%
  inner_join(grid_cells_sf, by = c("Cell_Id110" = "Cl_I110"))

poligonos_especies <- vert_assemblages_sf %>%
  group_by(Scientific.Name) %>%     # Agrupa por espécie
  summarise(geometry = st_union(geometry)) %>% # Une as células
  st_as_sf() %>%
  st_transform(st_crs(grid_cells_sf)) %>%
  mutate(TetrapodTraits.Name = Scientific.Name)

tetrapods_polygons_key <- species_list_tetrapods_filter %>%
  select(verbatim_name, speciesKey) %>%
  left_join(poligonos_especies, by = c('verbatim_name' = 'Scientific.Name')) %>%
  st_as_sf() %>%
  # algumas espécies ficaram com NA's pq nao tem dados espaciais
  remove_missing()

# GEOGRAPHIC SHAPEFILES ----
big_six <-  c("Brazil",
              "United States of America",
              "China",
              "Australia",
              "Russia", 
              "Canada")

world <- ne_countries(scale = "medium",
                      returnclass = "sf")
admin_divisions <- ne_states(
  country = big_six,
  returnclass = "sf",
)

# check projection
st_crs(world)
# colunas compartilhadas
intersect(names(admin_divisions), names(world))

admin_divisions_filter <- admin_divisions %>%
  select("admin", "adm0_a3", "name_en", "geometry")

geographic_shape_data <- world %>%
  select("admin", "adm0_a3", "name_en", "geometry") %>%
  filter(!(adm0_a3 %in% unique(admin_divisions_filter$adm0_a3))) %>%
  bind_rows(admin_divisions_filter) 

# Identificar feições inválidas
invalid_geometries <- geographic_shape_data[!st_is_valid(geographic_shape_data), ]
invalid_geometries <-  st_make_valid(invalid_geometries) 

geographic_shape_data <- geographic_shape_data %>%
  filter(st_is_valid(geographic_shape_data)) %>% # 454 unidades administrativas
  bind_rows(invalid_geometries) # 456 unidades administrativas

# ADD GLOBAL NORTH CLASSIFICATION ----
world_continents <- sf::st_read(dsn = file.path(
  "Shapefiles", "north_south_global"),
  layer='world-administrative-boundaries'
)

# Define a new grouping for major regions within continents (https://data.ipu.org/content/regional-groupings):
world_continents<-world_continents[,c("name", "continent", "region", "iso3", "geometry")]
names(world_continents)<-c("CountryName", "Continent", "Subregion", "ISO3", "geometry")
world_continents$Georegion<-factor(world_continents$Subregion,
                                   levels=c("Australia and New Zealand", "Caribbean", "Central America", 
                                            "Central Asia", "Eastern Africa", "Eastern Asia", "Eastern Europe",
                                            "Melanesia", "Micronesia", "Middle Africa", "Northern Africa",
                                            "Northern America", "Northern Europe", "Polynesia", "South-Eastern Asia",
                                            "South America", "Southern Africa", "Southern Asia", "Southern Europe",
                                            "Western Africa", "Western Asia", "Western Europe"),
                                   labels=c("Australia and New Zealand", "Latin America and Caribbean", "Latin America and Caribbean",
                                            "Central, East and South Asia", "Sub-Saharan Africa", "Central, East and South Asia", "Europe", 
                                            "Southeast Asia and Pacific Islands", "Southeast Asia and Pacific Islands", "Sub-Saharan Africa", "Near East and Northern Africa", 
                                            "Northern America", "Europe", "Southeast Asia and Pacific Islands", "Southeast Asia and Pacific Islands",
                                            "Latin America and Caribbean", "Sub-Saharan Africa",  "Central, East and South Asia", "Europe", 
                                            "Sub-Saharan Africa", "Near East and Northern Africa", "Europe"))
world_continents$GlobalNorth<-factor(world_continents$Subregion,
                                     levels=c("Australia and New Zealand", "Caribbean", "Central America", 
                                              "Central Asia", "Eastern Africa", "Eastern Asia", "Eastern Europe",
                                              "Melanesia", "Micronesia", "Middle Africa", "Northern Africa",
                                              "Northern America", "Northern Europe", "Polynesia", "South-Eastern Asia",
                                              "South America", "Southern Africa", "Southern Asia", "Southern Europe",
                                              "Western Africa", "Western Asia", "Western Europe"),
                                     labels=c("1", "0", "0",
                                              "0", "0", "0", "1", 
                                              "0", "0", "0", "0", 
                                              "1", "1", "0", "0",
                                              "0", "0",  "0", "1", 
                                              "0", "0", "1"))

# Adjust South Korea, Japan, and French Guiana as Global North:
world_continents$GlobalNorth<-as.character(world_continents$GlobalNorth)
world_continents[which(world_continents$CountryName=="Japan"),]$GlobalNorth<-"1"
world_continents[which(world_continents$CountryName=="Republic of Korea"),]$GlobalNorth<-"1"
world_continents[which(world_continents$CountryName=="French Guiana"),]$GlobalNorth<-"1"
world_continents[which(world_continents$CountryName=="Israel"),]$GlobalNorth<-"1"

# Make Russian Federation a separated region:
world_continents$Georegion<-as.character(world_continents$Georegion)
world_continents[world_continents$CountryName=="Russian Federation",]$Georegion<-"Russian Federation"

# Garantindo a mesma projecao
geographic_shape_data <- st_transform(geographic_shape_data,
                                      st_crs(world_continents))


world_globalnorth <- world_continents %>% 
  select(CountryName, GlobalNorth, Georegion, ISO3) %>%
  mutate(
    CountryName = if_else(CountryName == "Russian Federation",
                          "Russia",
                          CountryName)
  ) %>%
  st_drop_geometry()

# fazendo operacoes para recuperar o que eh GlobalNorth e o que nao eh
global_to_fill <- left_join(
  geographic_shape_data,
  world_globalnorth,
  by = c("admin" = "CountryName")
) 

# Filtrando os NA's
to_fill <- global_to_fill %>%
  filter(is.na(GlobalNorth)) %>%
  select(-GlobalNorth, -Georegion, -ISO3)

# Recuperando via ISO3
to_fill_2 <- left_join(
  to_fill,
  world_globalnorth,
  by = c("adm0_a3" = "ISO3")
) 

# Ainda nao conseguimos
to_fill_3 <- to_fill_2 %>%
  filter(is.na(GlobalNorth)) %>%
  select(-GlobalNorth, -Georegion)

to_fill_4 <- to_fill_3 %>%
  mutate(
    GlobalNorth = case_when(
      admin %in% c(
        "Saint Helena",
        "Saint Pierre and Miquelon",
        "Wallis and Futuna",
        "Saint Martin",
        "Saint Barthelemy",
        "French Southern and Antarctic Lands",
        "Aland",
        "Sint Maarten",
        "Curaçao",
        "Kosovo",
        "South Georgia and the Islands"
      ) ~ 1,
      
      admin %in% c(
        "Palestine",
        "Somaliland",
        "Northern Cyprus"
      ) ~ 0,
      
      admin %in% c(
        "Indian Ocean Territories",
        "Heard Island and McDonald Islands",
        "Ashmore and Cartier Islands",
        "Siachen Glacier"
      ) ~ 0,
      
      TRUE ~ NA_real_
    ),
    
    # Classificação das regiões
    Georegion = case_when(
      admin %in% c("Saint Helena", "Saint Pierre and Miquelon", 
                   "Wallis and Futuna", "Saint Martin", 
                   "Saint Barthelemy", "Aland", 
                   "Sint Maarten", "Curaçao") ~ "Europe",
      
      admin == "French Southern and Antarctic Lands" ~ "Australia and New Zealand",
      
      admin == "Kosovo" ~ "Europe",
      
      admin == "South Georgia and the Islands" ~ "Latin America and Caribbean",
      
      admin %in% c("Palestine", "Somaliland", "Northern Cyprus") ~ 
        case_when(
          admin == "Palestine" ~ "Near East and Northern Africa",
          admin == "Somaliland" ~ "Sub-Saharan Africa",
          admin == "Northern Cyprus" ~ "Near East and Northern Africa"
        ),
      
      admin %in% c("Indian Ocean Territories", 
                   "Heard Island and McDonald Islands",
                   "Ashmore and Cartier Islands") ~ "Australia and New Zealand",
      
      admin == "Siachen Glacier" ~ "Central, East and South Asia",
      
      TRUE ~ NA_character_
    )
  )

# remover NA's para juntar tudo em um dataframe soh
global_to_fill <- global_to_fill %>% 
  filter(!is.na(GlobalNorth)) %>%
  mutate(GlobalNorth = as.integer(GlobalNorth))
to_fill_2 <- to_fill_2 %>% 
  filter(!is.na(GlobalNorth)) %>%
  mutate(GlobalNorth = as.integer(GlobalNorth))
to_fill_4 <- to_fill_4 %>% 
  filter(!is.na(GlobalNorth)) %>%
  mutate(GlobalNorth = as.integer(GlobalNorth))


global_north_exceptions <- c(
  "United States Virgin Islands"
)

geographic_shape_data <- bind_rows(
  global_to_fill,
  to_fill_2,
  to_fill_4
) %>%
  mutate(
    GlobalNorth = if_else(
      admin %in% global_north_exceptions,
      1,
      GlobalNorth
    )
  ) %>%
  select(-"ISO3", -"CountryName") %>%
  rename("ISO3" = "adm0_a3") %>%
  filter(!is.na(name_en))

# Washington e Moscow aparecem duplicados
# Georgia aparece duas vezes no name_en pq tem um pais e um estado com 
# mesmo nome
# Identificar quais name_en são duplicados
# Unir polígonos duplicados

unified_data <- geographic_shape_data %>%
  filter(name_en %in% c("Altai Republic", "Washington", "Moscow", "Georgia")) %>%
  # renomear antes de agrupar
  mutate(
    name_en = case_when(
      # Washington mais ao leste → DC
      name_en == "Washington" & st_coordinates(st_centroid(geometry))[,1] > -90 ~ "Washington D.C.",
      # Washington mais a oeste → estado
      name_en == "Washington" & st_coordinates(st_centroid(geometry))[,1] <= -90 ~ "Washington",
      # Altai Republic mais ao norte → Altai Krai
      name_en == "Altai Republic" & st_coordinates(st_centroid(geometry))[,1] >= 85 ~ "Altai Republic",
      # Altai Republic mais ao sul → mantém Altai Republic
      name_en == "Altai Republic" & st_coordinates(st_centroid(geometry))[,1] < 85 ~ "Altai Krai",
      # Georgia no país → Georgia Country
      name_en == "Georgia" & ISO3 == "GEO" ~ "Georgia Country",
      # Moscow mantém mesmo nome
      TRUE ~ name_en
    )
  ) %>%
  # agrupar pelo novo nome
  group_by(name_en) %>%
  summarize(
    admin = first(admin),
    ISO3 = first(ISO3),
    GlobalNorth = first(GlobalNorth),
    Georegion = first(Georegion),
    geometry = st_union(geometry),
    .groups = "drop"
  )

geographic_shape_data <- geographic_shape_data %>%
  filter(!name_en %in% c("Altai Republic", "Washington", "Moscow", "Georgia")) %>%
  bind_rows(unified_data)

geographic_shape_data %>%
  group_by(name_en) %>%
  filter(n() > 1) %>%
  ungroup()
  

# SENSITIVITY ANALYSIS ----
## Expanded polygons ----
expand_species_polygons <- function(sf_polygons,
                                    buffer_distance_m,
                                    crs_metric = st_crs(grid_cells_sf)) {
  sf_polygons %>%
    st_transform(crs_metric) %>%
    mutate(
      geometry = if (buffer_distance_m > 0) {
        st_buffer(geometry, dist = buffer_distance_m)
      } else {
        geometry
      }
    ) %>%
    st_make_valid() 
}

scenarios <- tibble::tibble(
  scenario = c("base", "buffer_110km", "buffer_220km"),
  buffer_m = c(0, 110000, 220000)
)

tetrapods_polygons_expanded <- scenarios %>%
  mutate(
    polygons = map(
      buffer_m,
      ~ expand_species_polygons(
        sf_polygons = tetrapods_polygons_key,
        buffer_distance_m = .x
      )
    )
  ) %>%
  select(scenario, polygons) %>%
  unnest(polygons) %>%
  st_as_sf() %>%
  st_make_valid()

table(tetrapods_polygons_expanded$scenario)
# 32.002 poligonos

# only_land
only_land <- geographic_shape_data %>%
  st_union() %>%      # dissolve todos os estados
  st_transform(st_crs(grid_cells_sf)) %>%
  st_make_valid()

st_crs(only_land)
st_crs(tetrapods_polygons_expanded)

tetrapods_polygons_expanded_land <- st_intersection(
  tetrapods_polygons_expanded,
  only_land
)

table(tetrapods_polygons_expanded_land$scenario)
# 32002-31882 = 120 da base sairam
# 32002-31943 = 59 dos 110km sairam
# 32002-31966 = 36 dos 220km sairam

# conferindo quem saiu
completo <- tetrapods_polygons_expanded %>%
  filter(scenario == "buffer_110km") 
  
perdidos <- tetrapods_polygons_expanded_land %>%
  filter(scenario == "buffer_110km") 

sp_perdidas <- setdiff(completo$verbatim_name, perdidos$verbatim_name)

# Basicamente das 120 spp que sairam da base, são de ilhas!
mapview::mapview(tetrapods_polygons_expanded %>%
                   filter(verbatim_name %in% sp_perdidas & scenario == "buffer_110km"))

mapview::mapview(tetrapods_polygons_expanded %>%
                   filter(verbatim_name %in% "Cyanoramphus unicolor"))

tetrapods_polygons_sensitivity <- tetrapods_polygons_expanded_land %>%
  filter(scenario != "base") %>%
  arrange(verbatim_name)

# SAVE ----
# Base
save(tetrapods_polygons_key,
     file = file.path(
       local_directory,
       "00_raw_data",
       "tetrapods_polygons_key.RData"
     )
)

# Sensitivity analysis expanded 110km and 220km
save(tetrapods_polygons_key,
     file = file.path(
       local_directory,
       "00_raw_data",
       "tetrapods_polygons_sensitivity.RData"
     )
)

# Geographic file
save(geographic_shape_data,
     file = file.path(
       local_directory,
       "00_raw_data",
       "geographic_shape_data.RData"
     )
)
