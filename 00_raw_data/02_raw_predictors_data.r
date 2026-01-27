library(tidyverse)         # Data handling
library(CoordinateCleaner) # Access to institution
options(scipen = 999, digits = 15) 

# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time")

# Data loading
# V-Dem provides a multidimensional and disaggregated dataset that reflects the 
# complexity of the concept of democracy as a system of rule that goes beyond 
# the simple presence of elections. We distinguish between five high-level 
# principles of democracy: electoral, liberal, participatory, deliberative, 
# and egalitarian, and collect data to measure these principles
politic_data <- readRDS(file.path(
  local_directory, "00_raw_data",
  "V-Dem-CY-FullOthers-v15_rds", "V-Dem-CY-Full+Others-v15.rds")) %>%
  select("country_name", "country_text_id", "year", 
         "e_peedgini",    # Educational inequality, Gini (E)
         "e_pop",         # Populacao
         "e_gdppc",       # GPD per capta
         "e_miinterc",    # Armed conflict, internal
         "v2x_polyarchy", # Electoral democracy index
         "v2x_api",       # Additive polyarchy index 
         "v2x_mpi",       # Multiplicative polyarchy index
         "v2clacfree")    # Freedom of academic and cultural expression

# DOSE é um dataset global separado por unidades administrativas, no caso, estamos
# usando as unidades administrativas dos 6 países do mundo que tem cobertura de terra
# maior que 5%, o que aqui chamaremos de 'big six'. Com isso, esperamos identificar
# desigualdades regionais dentro desses países muito grandes. Do DOSE, iremos
# usar apenas os dados dos Big Six.
dose_data <- read.csv2(file.path(
  local_directory, "00_raw_data", "DOSE", "DOSE_V2.11.csv"), sep = ',') %>%
  mutate(pop_abs = as.numeric(pop),
         gdp_2015 = as.numeric(grp_pc_usd_2015)) %>%
  select(country, region, GID_0, GID_1, year, pop_abs, gdp_2015) 

# Maddison project
# The Maddison Project Database provides information on comparative economic 
# growth and income levels over the very long run.
maddison_data <- readxl::read_xlsx(file.path(
  local_directory, "00_raw_data", "Maddison", "mpd2023_web.xlsx"),
  sheet = "Full data") %>%
  filter(year > 1755) %>%
  rename(
    gdp_2011 = gdppc,
    macroregion = region
  ) 

# The Quality of Government Institute
qog_data <- read.csv2(file.path(
  local_directory, "00_raw_data", "The_QoG_Institute",
  "qog_std_ts_jan25.csv"), sep = ",") %>%
  select("cname", "ccodealp", "year", "ht_colonial", "gpi_gpi")

# A global gazetteer for biodiversity institutions from various sources, 
# including zoos, museums, botanical gardens, GBIF contributors, herbaria,
# university collections. 'institutions' data from CoordinateCleaner
data(institutions)

# Save socieconomic dataset
save(
  politic_data,
  dose_data,
  maddinson_data,
  qog_data,
  institutions,
  file = file.path(
    local_directory,
    "00_raw_data",
    "data_socieconomic.RData"
  )
)