## Description
This repository contains the data and scripts used to investigate the time to understand the geographical distribution of tetrapods in global scale. The data were collected from three main databases: the Global Biodiversity Information Facility (GBIF), speciesLink, and BioTime, supplemented by functional trait data from TetrapodTraits.

## Contents
- **`\00_raw_data`**: Contains R scripts to obtain the raw data used in the analyses.
- **`\01_data_cleaned`**: Contains R scripts used for data processing, cleaning, and join datasets.
- **`\02_data_analysis`**: In development.

## Data Collection
Occurrence data for terrestrial tetrapods were obtained from the following sources:
1. **Global Biodiversity Information Facility (GBIF)** [www.gbif.org](https://www.gbif.org/)
2. **speciesLink** [https://specieslink.net/](https://specieslink.net/)
3. **BioTime** (Dornelas et al., 2018)

To ensure data quality and taxonomic consistency across different sources, we used the Darwin Core taxonomy framework available via the `rgbif` package API (Chamberlain et al., 2024). Taxonomic matching was validated with a confidence interval greater than 95%.

## Methodology
The data were selected and standardized according to the following criteria:
- **Species selection**: We used TetrapodTraits (Moura et al., 2024) as the reference to select species, ensuring that they met the inclusion criteria: (1) included in global tetrapod phylogenies; (2) described at least 10 years ago; and (3) had functional traits available (e.g., body size, habitat, circadian activity).
- **Taxonomic standardization**: Species were homogenized using the `rgbif` package, specifically Darwin Core taxonomic backbone, to ensure name matching across different taxonomic sources, guaranteeing that occurrences were associated with the same taxonomic entity.
  
### Data Cleaning
The data underwent a cleaning process to remove records with imprecise or missing coordinates, following these steps:
1. Removal of records without coordinates or year.
2. Removal of records dated prior to 1900 due to low coordinate accuracy.
3. Application of the `CoordinateCleaner` pipeline (Zizka et al., 2019) to eliminate records with coordinates at country and/or capital centroids, located in research institutions, at sea, with equal absolute longitude and latitude, at GBIF headquarters
4. Exclusion of records with coordinate precision lower than 100 km or without precision data.
5. Duplicated coordinates

## How to Use the Data
In development

## License
This repository is not released for use. 
