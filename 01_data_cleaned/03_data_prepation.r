# Functions and packages ----
build_wallacean_time <- function(
    data_wallacean_all,
    completeness,                 # speciesKey + ua_total (n_distinct(name_en) por espécie)
    data_socieconomic_temporal,
    data_humandensity_temporal,
    countries_without_gdp = NULL, # data.frame com coluna name_en (opcional)
    min_year   = 1500,
    end_date   = as.POSIXct("2026-01-01", tz = "UTC"),
    origin_cal = as.POSIXct("1756-01-01", tz = "UTC"),
    seed       = NULL             # para reprodutibilidade do desempate aleatório
) {
  if (!is.null(seed)) set.seed(seed)
  yr_secs <- 86400 * 365.25
  
  # completeness scenarios
  comp <- completeness %>%
    distinct(speciesKey, ua_total) %>%
    mutate(
      ua_total_75 = as.integer(ceiling(ua_total * 0.75)),
      ua_total_50 = as.integer(ceiling(ua_total * 0.50))
    )
  
  # occurences + variables
  dat <- data_wallacean_all %>%
    filter(year > min_year)
  
  if (!is.null(countries_without_gdp)) {
    dat <- dat %>% filter(!name_en %in% countries_without_gdp$name_en)
  }
  
  dat <- dat %>%
    # espécie só pode ser encontrada depois de descrita
    mutate(year_modified = if_else(year < YearOfDescription, YearOfDescription, year)) %>%
    distinct(speciesKey, year_modified, name_en, .keep_all = TRUE) %>%
    mutate(
      year_modified = as.integer(year_modified),
      month = if_else(is.na(month), 1L, as.integer(month)),
      day   = if_else(is.na(day),   1L, as.integer(day)),
      date  = as.POSIXct(make_datetime(year_modified, month, day), tz = "UTC") + seconds(1)
    ) %>%
    left_join(comp, by = "speciesKey") %>%
    left_join(data_socieconomic_temporal, by = c("name_en", "year_modified" = "year")) %>%
    left_join(
      data_humandensity_temporal %>% select(speciesKey, year, HumanDensityTime),
      by = c("speciesKey", "year_modified" = "year")
    ) %>%
    select(speciesKey, scientificName, Class, Order, Family,
           name_en, date, year_modified, log_gdp_final, gdp_rank, gdp_percentile,
           gdp_percentile_country,
           BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
           HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus,
           ua_total, ua_total_75, ua_total_50,
           YearOfDescription, RangeSize, origin_of_data)
  
  # Chave aleatória só ordena os empates; o "empurrão" de +1s é o mesmo do
  # loop original, só que vetorizado: c_i = max(d_i, c_{i-1} + 1)
  dat <- dat %>%
    mutate(.rand = runif(n())) %>%
    arrange(speciesKey, date, .rand) %>%
    group_by(speciesKey) %>%
    distinct(name_en, date, .keep_all = TRUE) %>%
    mutate(
      .idx = row_number() - 1L,
      .sec = cummax(as.numeric(date) - .idx) + .idx,
      date = as.POSIXct(.sec, origin = "1970-01-01", tz = "UTC")
    ) %>%
    ungroup() %>%
    select(-.rand, -.idx, -.sec)
  
  # tabela de tempo até o evento
  tt <- dat %>%
    group_by(speciesKey) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      t.stop.date  = date,
      t.start.date = lag(t.stop.date),
      t.start.date = if_else(is.na(t.start.date),
                             as.POSIXct(make_date(YearOfDescription, 1, 1), tz = "UTC"),
                             t.start.date),
      t.stop.date  = if_else(is.na(t.stop.date), end_date, t.stop.date),
      event        = if_else(!duplicated(name_en), 1L, 0L),
      ua_acumulada = cumsum(!duplicated(name_en)),
      duracao_evento = as.numeric(difftime(t.stop.date, t.start.date, units = "secs")),
      duracao_evento = if_else(duracao_evento <= 0, 1, duracao_evento),
      t.start = cumsum(lag(duracao_evento, default = 0)) / yr_secs,
      t.stop  = t.start + duracao_evento / yr_secs,
      t.start.year = as.numeric(difftime(t.start.date, origin_cal, units = "secs")) / yr_secs + year(origin_cal),
      t.stop.year  = as.numeric(difftime(t.stop.date,  origin_cal, units = "secs")) / yr_secs + year(origin_cal)
    ) %>%
    ungroup()
  
  # corte após o evento terminal, por cenário
  cut_scenario <- function(df, ua_col) {
    df %>%
      filter(!is.na(.data[[ua_col]])) %>%
      group_by(speciesKey) %>%
      mutate(
        WallaceCompletude = event == 1L & ua_acumulada == .data[[ua_col]],
        # mantém linhas até (e incluindo) o primeiro evento terminal
        .antes = cumsum(lag(WallaceCompletude, default = FALSE))
      ) %>%
      filter(.antes == 0) %>%
      ungroup() %>%
      mutate(WallaceCompletude = as.integer(WallaceCompletude)) %>%
      select(speciesKey, Class, Order, Family, scientificName, name_en, date,
             t.start.date, t.stop.date, t.start, t.stop, t.start.year, t.stop.year,
             event, WallaceCompletude, log_gdp_final, gdp_rank, gdp_percentile,
             gdp_percentile_country,
             BodyLength_mm, BodyMass_g, Nocturnality, Verticality,
             HumanDensityTime, Latitude, Elevation, ETA50K, AssessedStatus, RangeSize,
             ua_total, ua_total_75, ua_total_50, origin_of_data, YearOfDescription)
  }
  
  list(
    c100 = cut_scenario(tt, "ua_total"),
    c75  = cut_scenario(tt, "ua_total_75"),
    c50  = cut_scenario(tt, "ua_total_50")
  )
}

install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    #renv::restore()
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "sf",
  "tidyverse",
  "arrow"
)

# executing install & load for each package
sapply(package_vec, install.load.package)

# Set directory 
rm(list = setdiff(ls(), c("build_wallacean_time",
                          "cut_scenario"))); gc()
local_directory <- file.path("C:",
                             "Users", "Mathias-190925-2", "OneDrive",
                             "Pesquisa", "wallacean_time") 

# Load data ----
# occurences data
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "dataset_occurences.RData")
)

rm(list = setdiff(ls(), c("local_directory",
                          "build_wallacean_time",
                          "list_per_admunit",
                          "data_wallacean_unnested"))); gc()
# trait data
load(file.path(
  local_directory,
  "00_raw_data",
  "trait_data.RData")
)

# socieconomic data
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_temporal.RData")
)
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_humandensity_temporal.RData")
)

# species list
load(
  file.path(
    local_directory,
    "00_raw_data",
    "tetrapods_list_sensitivity.RData"
  )
)

# Adm. unit completeness list ----
# Amphibia
rich_admunit_expected_amphibia <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Amphibia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_amphibia <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Amphibia") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_amphibia <- left_join(
  rich_admunit_expected_amphibia,
  rich_admunit_observed_amphibia,
  by = "name_en"
)
# Squamata
rich_admunit_expected_reptilia <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Reptilia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_reptilia <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Reptilia") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_reptilia <- left_join(
  rich_admunit_expected_reptilia,
  rich_admunit_observed_reptilia,
  by = "name_en"
)
# Aves
rich_admunit_expected_aves <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Aves") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_aves <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Aves") %>%
  filter(observed == TRUE) %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(observed_richness = n_distinct(speciesKey))

richness_completude_aves <- left_join(
  rich_admunit_expected_aves,
  rich_admunit_observed_aves,
  by = "name_en"
)
# Mammals
rich_admunit_expected_mammalia <- species_admin_all %>%
  select(name_en, all_species) %>%
  unnest(all_species) %>%
  filter(Class == "Mammalia") %>%
  select(name_en, speciesKey) %>%
  group_by(name_en) %>%
  summarise(polygons_richness = n_distinct(speciesKey))

rich_admunit_observed_mammalia <- species_admin_all %>%
  select(name_en, joined_df) %>%
  unnest(joined_df) %>%
  filter(Class == "Mammalia") %>%
  filter(observed == TRUE) %>%
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
  file = file.path(
    local_directory,
    "02_data_analysis",
    "richness_completude.RData")
)

## Species-level completeness ----
# Species-level completeness is based on the number of administrative units
# expected to contain each species according to its distribution polygons.
# For each species, we count the distinct administrative units where its
# polygon intersects the unit (all polygons) or reaches the specified minimum
# coverage threshold (0.1%, 0.5%, or 1%). These counts represent the expected
# geographic distribution of each species across administrative units and
# serve as the denominator for calculating species-level completeness.
completeness_all <- species_admin_occ_polygons %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  select(name_en, occ_polygons) %>%
  unnest(occ_polygons) %>%
  select(speciesKey, name_en) %>%
  distinct() %>%
  group_by(speciesKey) %>%
  summarise(
    ua_total = n_distinct(name_en),
    .groups = "drop"
  )

completeness01 <- species_admin_01perc %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  select(name_en, perc_01) %>%
  unnest(perc_01) %>%
  select(speciesKey, name_en) %>%
  distinct() %>%
  group_by(speciesKey)  %>%
  summarise(
    ua_total = n_distinct(name_en),
    .groups = "drop"
  )

completeness05 <- species_admin_05perc %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  select(name_en, perc_05) %>%
  unnest(perc_05) %>%
  select(speciesKey, name_en) %>%
  distinct() %>%
  group_by(speciesKey)  %>%
  summarise(
    ua_total = n_distinct(name_en),
    .groups = "drop"
  )

completeness1 <- species_admin_1perc %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  select(name_en, perc_1) %>%
  unnest(perc_1) %>%
  select(speciesKey, name_en) %>%
  distinct() %>%
  group_by(speciesKey)  %>%
  summarise(
    ua_total = n_distinct(name_en),
    .groups = "drop"
  ) 

# Species list filtering ----
# Without filtering, the occurrences could be counted in an administrative
# unit that is not part of the species' expected distribution according to
# its distribution polygons. Therefore, temporal completeness at the species
# level should be calculated only across administrative units considered
# expected. This prevents occurrences in unexpected units from artificially
# contributing to the completeness of the expected distribution.
#
# For example: 
# SpeciesKey 5231640 (Carduelis tristis, American Goldfinch): occurrence records
# were originally assigned to New Mexico, but the current administrative
# boundary layer assigns them to Oklahoma. The records are approximately
# 68–72 m from the New Mexico boundary.
expected_admin_all <- species_admin_occ_polygons %>%
  select(name_en, occ_polygons) %>%
  unnest(occ_polygons) %>%
  select(speciesKey, name_en) %>%
  distinct()

expected_admin_01 <- species_admin_01perc %>%
  select(name_en, perc_01) %>%
  unnest(perc_01) %>%
  select(speciesKey, name_en) %>%
  distinct()

expected_admin_05 <- species_admin_05perc %>%
  select(name_en, perc_05) %>%
  unnest(perc_05) %>%
  select(speciesKey, name_en) %>%
  distinct()

expected_admin_1 <- species_admin_1perc %>%
  select(name_en, perc_1) %>%
  unnest(perc_1) %>%
  select(speciesKey, name_en) %>%
  distinct()

# Create time-to-event table ----
## Uncut dataset ----
data_wallacean_all <- data_wallacean_unnested %>%
  # Remove countries without economic data
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  # Retain occurrences within expected administrative units
  semi_join(expected_admin_all, by = c("speciesKey", "name_en"))

df_time_wallacean <- build_wallacean_time(
  data_wallacean_all, 
  completeness_all,
  data_socieconomic_temporal,
  data_humandensity_temporal,
  countries_without_gdp = countries_without_gdp,
  seed = 42
)

save(df_time_wallacean,
     file = file.path(
       local_directory,
       "02_data_analysis",
       "data_wallacean_time.RData")
)

## Cut 0.1% ----
rm(list = setdiff(ls(), c("local_directory",
                          "countries_without_gdp",
                          "completeness01",
                          "expected_admin_01",
                          "completeness05",
                          "expected_admin_05",
                          "completeness1",
                          "expected_admin_1",
                          "build_wallacean_time",
                          "data_socieconomic_temporal",
                          "data_humandensity_temporal",
                          "data_wallacean_unnested"))); gc()

data_wallacean_01 <- data_wallacean_unnested %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  semi_join(expected_admin_01, by = c("speciesKey", "name_en"))

df_time_wallacean01 <- build_wallacean_time(
  data_wallacean_01, 
  completeness01,
  data_socieconomic_temporal,
  data_humandensity_temporal,
  countries_without_gdp = countries_without_gdp,
  seed = 42
)

save(df_time_wallacean01,
     file = file.path(
       local_directory,
       "02_data_analysis",
       "data_wallacean_time01.RData")
)

## Cut 0.5% ----
rm(list = setdiff(ls(), c("local_directory",
                          "countries_without_gdp",
                          "completeness05",
                          "expected_admin_05",
                          "completeness1",
                          "expected_admin_1",
                          "data_socieconomic_temporal",
                          "data_humandensity_temporal",
                          "data_wallacean_unnested",
                          "build_wallacean_time"))); gc()

data_wallacean_05 <- data_wallacean_unnested %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  semi_join(expected_admin_05, by = c("speciesKey", "name_en"))

df_time_wallacean05 <- build_wallacean_time(
  data_wallacean_05, 
  completeness05,
  data_socieconomic_temporal,
  data_humandensity_temporal,
  countries_without_gdp = countries_without_gdp,
  seed = 42
)

save(df_time_wallacean05,
     file = file.path(
       local_directory,
       "02_data_analysis",
       "data_wallacean_time05.RData")
)

## Cut 1% ----
rm(list = setdiff(ls(), c("local_directory",
                          "countries_without_gdp",
                          "completeness1",
                          "expected_admin_1",
                          "data_socieconomic_temporal",
                          "data_humandensity_temporal",
                          "data_wallacean_unnested",
                          "build_wallacean_time"))); gc()

data_wallacean_1 <- data_wallacean_unnested %>%
  filter(!name_en %in% countries_without_gdp$name_en) %>%
  semi_join(expected_admin_1, by = c("speciesKey", "name_en"))

df_time_wallacean1 <- build_wallacean_time(
  data_wallacean_1, 
  completeness1,
  data_socieconomic_temporal,
  data_humandensity_temporal,
  countries_without_gdp = countries_without_gdp,
  seed = 42
)

save(df_time_wallacean1,
     file = file.path(
       local_directory,
       "02_data_analysis",
       "data_wallacean_time1.RData")
)