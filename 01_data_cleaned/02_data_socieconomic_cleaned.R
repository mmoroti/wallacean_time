# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    #renv::restore()
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

# preencher os anos intermediarios com base na variacao entre decadas
adjust_country_years <- function(df_country){
  
  first_year <- min(df_country$year, na.rm = TRUE)
  
  # PARTE 1: ajustar anos < 1820
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
  
  # PARTE 2: criar 2023–2025
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
  
  # PARTE 3: marcar edge_imputed
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

# names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "sf",
  "tidyverse",
  "units",
  "countrycode"
)

# executing install & load for each package
sapply(package_vec, install.load.package)

# Set directory 
local_directory <- file.path("C:",
                             "Users", "Mathias-190925-2", "OneDrive",
                             "Pesquisa", "wallacean_time") 

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
                          "adjust_country_years",
                          "dose_data",
                          "maddinson_data",
                          "maddison_macroregion_data", 
                          "maddison_ppp_data",
                          "politic_data",
                          "qog_data",
                          "institutions",
                          "geographic_shape_data",
                          "list_per_admunit",
                          "onu_alignment_data",
                          "data_wallacean_unnested"))); gc()

# HANDLING DATA ----
admunit <- geographic_shape_data %>%
  select(name_en, ISO3) %>%
  distinct() %>%
  arrange(ISO3) %>%
  sf::st_drop_geometry()

politic_data_summary <- politic_data %>%
  group_by(country_text_id) %>%
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
    country_text_id = recode_factor(
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

# First, we will attempt to calculate a mean value for each administrative unit.
# DOSE provides administrative-unit-level data for several countries, but we
# only subdivided six countries into administrative units. Therefore, we will
# split the data into two datasets: one containing the administrative regions
# and another containing countries as whole units.

# Harmonization of names across datasets.
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

# These two administrative units are treated as a single unit in our dataset
# (Moscow region). Therefore, I will merge region and city into one unit by
# combining 'moscow region' and 'moscow city' into 'moscow'.
df_moscow <- dose_data %>%
  filter(region %in% c("Moscow city", "Moscow region")) %>%
  group_by(year) %>%
  summarise(
    region = "Moscow",
    gdp_2015 = sum(gdp_2015, na.rm = TRUE),
    pop_abs = sum(pop_abs, na.rm = TRUE),
    country = "Russia",
    GID_0 = "RUS"
  ) 

# Regional data
region_data <- dose_data %>%
  filter(!region %in% c("Moscow city", "Moscow region")) %>%
  bind_rows(df_moscow) %>%
  filter(country %in% c("Brazil", "Russia", "China", 
                        "Australia", "USA", "Canada")) %>%
  distinct(region, year, .keep_all = TRUE) %>% 
  # Conversion of GDP per capita to PPP-adjusted US dollars.
  left_join(maddison_ppp_data, by = c("GID_0" = "countrycode", "year")) %>%
  mutate(gdp_ppp_2015 = gdp_2015 / pl_gdpo) %>%
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
  left_join(region_data, by = "region_tolower") 

# DOSE has sparse coverage for countries in Africa and the Middle East,
# where Maddison (which uses World Bank and UN data for the recent period)
# will be essential to complement its global average.
maddinson_data_edit <- maddinson_data %>%
  filter(!countrycode %in% c("BRA", "RUS", "CHN", "AUS", "CAN", "USA")) %>%
  mutate(pop = pop * 1000) %>%
  group_by(countrycode) %>%
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

# Obtain institutions and institutions of each administrative unit.
# Here, I will retrieve large countries and those without a country classification
# to determine where they fall within our polygons. This will allow me to identify
# the administrative units in which they are located and classify previously
# unclassified institutions.
institutions_sf <- institutions %>%
  filter(is.na(country) | country %in% c("BRA", "CHN", "RUS", "USA", "CAN", "AUS")) %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = st_crs(geographic_shape_data))

# Count institutions within each name_en.
institutions_with_location <- institutions_sf %>%
  st_join(geographic_shape_data %>% select(name_en)) %>%
  filter(!is.na(name_en)) %>%
  st_drop_geometry() %>%
  group_by(name_en) %>%
  summarise(n_institutions = n())

# Now I can group the data by country, but first I want to identify the exact
# name used in our dataset (name_en) for each ISO3 code.
nameen_list <- geographic_shape_data %>%
  select(name_en, ISO3) %>%
  filter(!ISO3 %in% c("BRA", "CHN", "RUS", "USA", "CAN", "AUS")) %>%
  st_drop_geometry() %>%
  mutate(
    ISO3 = countrycode(
      name_en,
      origin = "country.name",
      destination = "iso3c",
      nomatch = NA  # NA for not found
    )
  )

df_institutions <- institutions %>%
  filter(!country %in% c("BRA", "CHN", "RUS", "USA", "CAN", "AUS")  & !is.na(country)) %>%
  group_by(country) %>%
  summarise(n_institutions = n()) %>%
  left_join(nameen_list, by = c("country" = "ISO3")) %>%
  bind_rows(institutions_with_location) %>%
  filter(country != "GEO" | is.na(country)) %>%  
  group_by(name_en) %>%                          
  summarise(
    n_institutions = sum(n_institutions, na.rm = TRUE),  
    country = first(na.omit(country)),          
    GEO = any(country == "GEO", na.rm = TRUE)   
  ) %>%
  ungroup() %>%
  select(name_en, n_institutions)

# Add ONU Alignment
onu_alignment <- onu_alignment_data %>%
  group_by(iso3c) %>%
  # se quiser pegar os 10 ultimos anos
  slice_tail(n = 10) %>%
  summarise(
    #min_year_politic = min(year, na.rm = TRUE),
    #max_year_politic = max(year, na.rm = TRUE),
    IdealPointFP_mean = mean(IdealPointFP, na.rm = TRUE)
  ) %>%
  ungroup()

# Obtain 
data_sociopolitic <- bind_rows(
  country_sociopolitic,
  admunit_sociopolitic
) %>%
  select(-region_tolower) %>%
  #filter(ISO3 != "GEO") %>%
  left_join(geographic_shape_data %>% 
              select(name_en, area_km2, GlobalNorth) %>%
              st_drop_geometry(),
            by = "name_en") %>%
  distinct() %>%
  mutate(
    gdp_mean = coalesce(gdp_2011_mean, gdp_2015_mean),
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
  ) %>%
  left_join(onu_alignment, by = c("ISO3" = "iso3c")) 

save(
  data_sociopolitic,
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_agregatte.RData")
)

# AGGREGATING THE DATA OVER TIME ----
## GDP over time ----
occ_year <- data_wallacean_unnested %>%
  filter(year > 1500) %>% 
  # A species can only be found after it has been formally described.
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
    name_en = recode_factor(
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
      TRUE ~ ISO3   # keeps the original value if it is not in the list
    )
  )

maddison_filtered <- maddinson_data %>%
  mutate(pop = pop * 1000) %>%
  left_join(first_year, by = c("countrycode" = "ISO3")) %>%
  filter(year >= first_year) %>%
  mutate(log_gdp = log(gdp_2011)) %>%
  select(-pop, -first_year, -last_year, -n_points)

df_interp <- maddison_filtered %>%
  # Group by country to interpolate each series separately
  group_by(country) %>%
  group_modify(~{
    # Creates a complete sequence of years between the minimum and maximum available
    years_full <- tibble(
      year = seq(min(.x$year), max(.x$year), by = 1)
    )
    # Combine the complete sequence with the original country data
    # This ensures that missing years appear as NA
    df_full <- years_full %>%
      left_join(.x, by = "year") %>%
      arrange(year)
    # Performs log-linear interpolation:
    # approx() uses only the years where log_gdp is not NA
    interp_vals <- approx(
      x = df_full$year[!is.na(df_full$log_gdp)],  # years observed
      y = df_full$log_gdp[!is.na(df_full$log_gdp)],  # values observed
      xout = df_full$year  # anos completos
    )$y
    # Adds interpolated columns and an interpolation indicator
    df_full %>%
      mutate(
        # log of interpolated GDP
        log_gdp_interp = interp_vals,
        # Interpolated GDP per capita (coming back from log)
        gdp_pc_interp = exp(log_gdp_interp),
        # Indicator: 1 = interpolated, 0 = observed
        # Based on the original value (before interpolation)
        is_interpolated = if_else(is.na(log_gdp), 1, 0)
      )
  }) %>%
  # Remove grouping
  ungroup() #%>%

# Plot
#ggplot(df_interp %>% filter(country == "Libya"), aes(x = year)) +
#  geom_point(aes(y = log_gdp), alpha = 0.6, color = "blue") + 
#  # Intervalo da Primeira Guerra
#  geom_rect(aes(xmin = 1914, xmax = 1918, ymin = -Inf, ymax = Inf),
#            fill = "grey70", alpha = 0.1, inherit.aes = FALSE) +
#  # Intervalo da Segunda Guerra
#  geom_rect(aes(xmin = 1939, xmax = 1945, ymin = -Inf, ymax = Inf),
#            fill = "grey50", alpha = 0.1, inherit.aes = FALSE) +
#  geom_line(aes(y = log_gdp_interp), color = "red", linewidth = 0.8) + 
#  theme_minimal() +
#  labs(
#    title = "Log-linear interpolation of GDP per capta",
#    x = "Year",
#    y = "log(GDP per capta)",
#    color = "Legenda"
#  ) 

# Macroregions
head(maddison_macroregion_data)

macro_long <- maddison_macroregion_data |>
  pivot_longer(-year,
               names_to = "macroregion",
               values_to = "gdp_region")

macro_interp <- macro_long %>%
  mutate(log_gdp = log(gdp_region)) %>%
  group_by(macroregion) %>%
  group_modify(~{
    years_full <- tibble(
      year = seq(min(.x$year), max(.x$year), by = 1)
    )
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
        log_gdp_interp_macro = interp_vals,
        gdp_pc_interp_macro = exp(log_gdp_interp_macro),
        is_interpolated_macro = if_else(is.na(log_gdp), 1, 0),
        is_macroregion = 1,
      )
  }) %>%
  ungroup() %>%
  select(-gdp_region, -log_gdp)

# Combining data across both scales
df_economy <- left_join(df_interp, macro_interp,
                        by = c("macroregion", "year")) 

# Former Soviet republics that were part of the Russian territory in the
# historical data from 1950 until they began reporting their own data.
# Therefore, we will handle them separately.
ex_ussr <- c(
  "Armenia", "Azerbaijan", "Belarus", "Georgia", 
  "Kazakhstan", "Kyrgyzstan", "Latvia", "Lithuania", 
  "Republic of Moldova", "Tajikistan", "Turkmenistan", 
  "Ukraine", "Uzbekistan"
)

# Include regional macroeconomic data up to 1950 to ensure that all regions
# are compared on a consistent basis. Most GDP data start in 1950, and since
# we will construct annual rankings, measurements must be comparable across
# units within each year. Therefore, the fairest approach is to compare
# regions with regions before 1950 and countries with countries from 1950 onward.
df_rest <- df_economy %>% filter(!country %in% ex_ussr)

df_rest_economy <- df_rest %>%
  mutate(
    economic_regime = if_else(year <= 1949, "macro_period", "national_period"),
    gdp_source = case_when(
      economic_regime == "macro_period"   ~ "macroregion",
      economic_regime == "national_period" & !is.na(log_gdp_interp) ~ "country",
      economic_regime == "national_period" & is.na(log_gdp_interp)  ~ "macroregion_fallback",
      TRUE ~ NA_character_
    ),
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

# Here, we will use Russian Federation data to fill in the missing values
# for countries that were part of the former Soviet Union until their own
# indicators become available.
first_year_df <- maddinson_data %>%
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

# Keep only the years and GDP values for Russia.
russia_gdp <- df_economy %>%
  filter(country == "Russian Federation") %>%
  select(year, log_gdp_interp, gdp_pc_interp)
# For former Soviet Union countries.
df_ussr <- df_economy %>% 
  filter(country %in% ex_ussr) %>%
  left_join(first_year_df, by = "country") %>%
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

# edge_imputed: years with no available data, for which we carried forward the
# most recent available values (before 1820 and after 2022).
# interpolated: missing values interpolated between years with observed data.
# missing_country_estimated: occurrence records were available, but the country
# had no observed data for that period; regional data were used instead.
# observed: values directly observed in Maddison.
# regional_imputed: regional data were used for all units from 1820 to 1950.
# sovereign_imputed: the country belonged to the former USSR, so Russian data
# were used as a proxy until country-specific data became available.

# Separate administrative units from countries, as they require different
# treatments.
df_economy_bigsix <- df_economy_all %>%
  filter(countrycode %in% c("BRA", "USA", "RUS", "CHN", "AUS", "CAN"))

df_economy_countries <- df_economy_all %>%
  filter(!countrycode %in% c("BRA", "USA", "RUS", "CHN", "AUS", "CAN"))

# Region data
region_data <- dose_data %>%
  filter(!region %in% c("Moscow city", "Moscow region")) %>%
  bind_rows(df_moscow) %>%
  filter(country %in% c("Brazil", "Russia", "China", 
                        "Australia", "USA", "Canada")) %>%
  distinct(region, year, .keep_all = TRUE) %>% 
  left_join(maddison_ppp_data, by = c("GID_0" = "countrycode", "year")) %>%
  mutate(gdp_ppp_2015 = gdp_2015 / pl_gdpo) %>%
  mutate(region_tolower = tolower(region)) %>%
  mutate(region_tolower = str_replace_all(region_tolower, map_names)) %>%
  mutate(
    GID_0 = case_when(
      region == "Moscow" ~ "RUS",
      TRUE ~ GID_0
    )
  )

# Brazil
# Brazilian states have data available from 1970 onward.
# Use Goiás values from 1970 to 1984 for Tocantins.
# Use Mato Grosso values from 1970 to 1974 for Mato Grosso do Sul.
brazil_data <- region_data %>% 
  filter(GID_0 == "BRA" & year >= 1970)
# 1. Reference values
goias_ref <- brazil_data %>%
  filter(region == "Goiás", year >= 1970, year <= 1988) 
mt_ref <- brazil_data %>%
  filter(region == "Mato Grosso", year >= 1970, year <= 1974) 
# 2. Create artificial lines for Tocantins and MS
tocantins_fill <- goias_ref %>%
  mutate(region = "Tocantins")
ms_fill <- mt_ref %>%
  mutate(region = "Mato Grosso Do Sul")
# 3. Join everything to the original dataset
brazil_fixed <- brazil_data %>%
  bind_rows(tocantins_fill, ms_fill) %>%
  arrange(region, year)

# Canada
# Before April 1, 1999, what is now known as Nunavut did not exist as a
# separate administrative entity. It was part of the Northwest Territories.
# In DOSE, this unit is referred to as "Northwest Territories Including Nunavut".
# Therefore, we will use values from this series for the territories until
# they have their own data.
canada_edit_data <- region_data %>% 
  filter(region == "Northwest Territories Including Nunavut" & year > 1983)
nunavut_ref <- canada_edit_data %>%
  mutate(region = "Nunavut") 
northwest_ref <- canada_edit_data %>%
  mutate(region = "Northwest Territories") 
# Join to the original dataset, excluding the territory 'Northwest Territories
# Including Nunavut'
canada_fixed <- region_data %>%
  filter(GID_0 == "CAN" & region != "Northwest Territories Including Nunavut") %>%
  bind_rows(nunavut_ref, northwest_ref) %>%
  arrange(region, year)

# China
# In 1988, Hainan was officially separated from Guangdong and established
# as an independent province, becoming a separate administrative unit.
china_edit_data <- region_data %>% 
  filter(region == "Guangdong" & year < 1978)
hainan_ref <- china_edit_data %>%
  mutate(region = "Hainan") 
china_fixed <- region_data %>%
  filter(GID_0 == "CHN") %>%
  bind_rows(hainan_ref) %>%
  arrange(region, year)

# Russia
# Ingushetia has existed as a separate administrative unit only since 1994,
# and Chechnya since 2005. Before these dates, the GDP of both regions was
# reported together under a single administrative unit.
russia_edit_data <- region_data %>% 
  filter(region == "Ingushetia" & year < 2005)
chechen_ref <- russia_edit_data %>%
  mutate(region = "Chechen") 
russia_fixed <- region_data %>%
  filter(GID_0 == "RUS") %>%
  bind_rows(chechen_ref) %>%
  arrange(region, year)

# Australia has GDP data available from 1970 onward.
australia_fixed <- region_data %>%
  filter(GID_0 == "AUS" & year > 1953)

# Combining data for administrative units.
# Australia and the USA have data starting from the same year, so no adjustments
# are needed for these countries.
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
  mutate(log_gdp = log(gdp_ppp_2015)) %>%
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

# populating the macro-region
ref_meta <- df_economy_bigsix_all %>%
  select(GID_0, macroregion) %>% 
  distinct()

# Creating a time series for the final year of each administrative unit up to 2025
last_rows_extended <- region_data_all %>%
  group_by(GID_0) %>%
  filter(year == max(year)) %>% # get the last row for each country
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

data_socieconomic_temporal_raw <- bind_rows(
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

data_socieconomic_temporal_raw %>%
  group_by(name_en) %>%
  summarise(
    countrycode       = first(countrycode),
    first_year = min(year),
    last_year  = max(year),
    n_points   = n(),
  ) %>% View()

# We need to harmonize the names of administrative units across the
# socioeconomic databases and our own occurrence database.
setdiff(occ_year$name_en,
        data_socieconomic_temporal_raw$name_en)

setdiff(data_socieconomic_temporal_raw$name_en,
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

data_socieconomic_temporal <- data_socieconomic_temporal_raw %>%
  mutate(
    name_en = recode_factor(name_en, !!!harmonize_names)
  ) 

setdiff(occ_year$name_en,
        data_socieconomic_temporal$name_en)

data <- left_join(
  occ_year,
  data_socieconomic_temporal,
  by = c("name_en", "year")) 

countries_without_gdp <- data %>%
  filter(is.na(gdp_pc_final)) %>% 
  distinct(name_en) 

# Create national GDP dataset
# Identify the six countries with administrative-level data.
bigsix_countries <- data_socieconomic_temporal %>%
  filter(countrycode %in% c("BRA", "RUS", "CAN", "USA", "CHN", "AUS")) %>%
  select(countrycode, name_en) %>%
  distinct()

# Retrieve country-level economic data from 1950 onward.
bigsix_country_data <- df_economy_all %>%
  filter(
    countrycode %in% c("BRA", "RUS", "CAN", "USA", "CHN", "AUS"),
    year >= 1950
  ) %>%
  select(-name_en)

# Associate country names with country-level data.
bigsix_country_post1950 <- left_join(
  bigsix_countries,
  bigsix_country_data,
  by = "countrycode"
)

# Retain administrative-level data before 1950 and all other countries
# throughout the entire time series.
bigsix_regional_pre1950 <- data_socieconomic_temporal %>%
  filter(
    !(countrycode %in% c("BRA", "RUS", "CAN", "USA", "CHN", "AUS") &
        year >= 1950)
  )

# Combine regional data before 1950 and country-level data from 1950 onward.
df_economy_country_comparable <- bind_rows(
  bigsix_country_post1950,
  bigsix_regional_pre1950
) %>%
  arrange(countrycode, year)

# Calculate the original annual GDP rankings.
df_economy_original_ranks <- data_socieconomic_temporal %>%
  group_by(year) %>%
  mutate(
    gdp_rank = rank(log_gdp_final) / n(),
    gdp_percentile = percent_rank(log_gdp_final)
  ) %>%
  ungroup()

# Calculate annual GDP rankings using the harmonized comparison units.
df_economy_comparable_ranks <- df_economy_country_comparable %>%
  group_by(year) %>%
  mutate(
    gdp_rank_country = rank(log_gdp_final) / n(),
    gdp_percentile_country = percent_rank(log_gdp_final)
  ) %>%
  ungroup() %>%
  select(
    name_en,
    year,
    gdp_percentile_country,
    gdp_rank_country
  )

# Compare original and harmonized GDP rankings.
data_socieconomic_temporal <- left_join(
  df_economy_original_ranks,
  df_economy_comparable_ranks,
  by = c("name_en", "year")
)

data_socieconomic_temporal %>%
  filter(year < 1950) %>%
  summarise(
    rank_identical = all(gdp_rank == gdp_rank_country, na.rm = TRUE),
    percentile_identical = all(gdp_percentile == gdp_percentile_country, na.rm = TRUE)
  )

# Save temporal series
save(
  data_socieconomic_temporal,
  countries_without_gdp,
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_temporal.RData")
)

## Human density over time ----
load(
  file.path(
    local_directory,
    "00_raw_data",
    "df_humandensity.RData"
  )
)

# Aggregate human density data by species and year (median of the polygons)
df_humandensity_agg <- df_humandensity %>%
  mutate(median = ifelse(is.na(median), 0, median),
         mean = ifelse(is.na(mean), 0, mean)) %>%
  group_by(speciesKey, year) %>%
  summarise(
    median = median(median, na.rm = TRUE),
    mean = mean(mean, na.rm = TRUE),
    n_polygons = n(),  
    .groups = "drop"
  )

data_humandensity_temporal <- df_humandensity_agg %>%
  group_by(speciesKey) %>%
  group_modify(~{
    years_full <- tibble(
      year = seq(min(.x$year), 2025, by = 1)  # vai até 2025
    )
    df_full <- years_full %>%
      left_join(.x, by = "year") %>%
      arrange(year)
    valor_2017 <- df_full %>%
      filter(year <= 2017, !is.na(median)) %>%
      arrange(desc(year)) %>%
      slice(1) %>%
      pull(median)
    if(length(valor_2017) > 0) {
      df_full <- df_full %>%
        mutate(
          median = case_when(
            year >= 2017 & is.na(median) ~ valor_2017,
            TRUE ~ median
          )
        )
    }
    interp_vals <- approx(
      x = df_full$year[!is.na(df_full$median)],
      y = df_full$median[!is.na(df_full$median)],
      xout = df_full$year
    )$y
    df_full %>%
      mutate(
        median_interp = interp_vals,
        is_interpolated = if_else(is.na(median), 1, 0)
      )
  }) %>%
  ungroup() %>%
  rename(HumanDensityTime = median_interp)
# 5217144 check

save(
  data_humandensity_temporal,
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_humandensity_temporal.RData")
)