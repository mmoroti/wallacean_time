library(tidyverse)
library(sf)
library(stringi)
library(betareg)
library(cowplot)
library(ggstats)
library(ggrepel)
library(ggpubr)
library(ggdist) # for shadeable density slabs
library(gghalves) # for half-half geoms
library(car)
library(broom)
library(marginaleffects)
library(broom.helpers)
library(tidyr)
library(GGally)

rm(list=ls()); gc() # clean local enviroment

# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 
drive <- file.path("G:", "Meu Drive", "Artigos",
                   "wallace_time", "Figures")

## Function ----
plot_bivariate_map <- function(geographic_shape_data,
                               richness_completude_amphibia,
                               title = NULL) {
  
  # 1) porcentagem de completude
  richness_completude <- richness_completude_amphibia %>%
    mutate(percent_wallace = observed_richness * 100 / polygons_richness)
  
  # 2) join com shapefile
  dados_join <- geographic_shape_data %>%
    left_join(richness_completude, by = c("name_en" = "name_en")) %>%
    st_transform(6933)
  
  # 3) quantis (usar só >0 e não NA para richness)
  rich_valid <- dados_join$polygons_richness[
    !is.na(dados_join$polygons_richness) & dados_join$polygons_richness > 0
  ]
  
  rich_q_raw <- quantile(rich_valid, probs = seq(0, 1, length.out = 5))
  # converter para inteiros e garantir monotonicidade
  rich_breaks <- as.integer(floor(rich_q_raw))
  # garantir monotonicidade estrita
  for(i in 2:length(rich_breaks)){
    if(rich_breaks[i] <= rich_breaks[i-1]){
      rich_breaks[i] <- rich_breaks[i-1] + 1
    }
  }
  rich_breaks[1] <- min(rich_valid)
  rich_breaks[length(rich_breaks)] <- max(rich_valid)
  
  # percentuais (continua contínuo)
  perc_breaks <- quantile(dados_join$percent_wallace,
                          probs = seq(0, 1, length.out = 5),
                          na.rm = TRUE)
  
  # funções de label
  make_labels_int <- function(breaks, suffix = "") {
    paste0(
      breaks[-length(breaks)],
      "–",
      breaks[-1],
      suffix
    )
  }
  make_labels_pct <- function(breaks, suffix = "") {
    paste0(
      round(breaks[-length(breaks)], 1),
      "–",
      round(breaks[-1], 1),
      suffix
    )
  }
  
  rich_labs <- make_labels_int(rich_breaks, " spp")
  perc_labs <- make_labels_pct(perc_breaks, "%")
  
  # 4) criar categorias (NA para sem dado / zero)
  dados_bi <- dados_join %>%
    mutate(
      richness_cat = case_when(
        is.na(polygons_richness) | polygons_richness == 0 ~ NA_character_,
        TRUE ~ as.character(cut(polygons_richness,
                                breaks = rich_breaks,
                                include.lowest = TRUE,
                                labels = rich_labs))
      ),
      percent_cat = case_when(
        is.na(percent_wallace) ~ NA_character_,
        TRUE ~ as.character(cut(percent_wallace,
                                breaks = perc_breaks,
                                include.lowest = TRUE,
                                labels = perc_labs))
      ),
      bi_class = ifelse(is.na(richness_cat) | is.na(percent_cat),
                        NA_character_,
                        paste0(richness_cat, " | ", percent_cat))
    )
  
  # 5) paleta 4x4 — evitar coluna com mesmo nome 'color' para não sobrescrever
  rich_labs_factor <- factor(rich_labs, levels = rich_labs, ordered = TRUE)
  perc_labs_factor <- factor(perc_labs, levels = perc_labs, ordered = TRUE)
  
  bi_pal <- expand.grid(
    richness_cat = rich_labs_factor,
    percent_cat  = perc_labs_factor,
    stringsAsFactors = FALSE
  ) %>%
    arrange(richness_cat, percent_cat) %>%
    mutate(
      bi_class = paste0(richness_cat, " | ", percent_cat),
      fill = c(
        "#4e3d43", "#535e68", "#577e8b", "#5b9cad",
        "#724448", "#786970", "#7e8d96", "#84afba",
        "#944b4d", "#9c7478", "#a59ba0", "#acc1c6",
        "#b65252", "#c07f7f", "#caaaaa", "#d3d3d3"
      )
    )
  
  # join (bi_class é chave)
  dados_bi <- dados_bi %>%
    left_join(bi_pal %>% select(bi_class, fill), by = "bi_class")
  
  # 6) definir cor padrão (cinza claro) para quem não entrou na classificação
  dados_bi <- dados_bi %>%
    mutate(fill = ifelse(is.na(fill), "#ffffff", fill))
  
  # 7) MAPA
  bbox <- st_bbox(dados_bi)
  bbox_poly <- st_as_sfc(bbox)
  
  map_plot <- ggplot(dados_bi) +
    geom_sf(aes(fill = fill), color = "black", size = 0.15) +
    geom_sf(data = bbox_poly, fill = NA, color = "black", linewidth = 0.3) +
    scale_fill_identity() +
    labs(title = title) +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"]),
      expand = FALSE
    ) +
    theme_minimal() +
    theme(
      text = element_text(color = "black"),
      legend.position = "none",
      panel.border = element_blank(),
      axis.text = element_text(size = 8, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      axis.ticks.length = unit(2, "mm"),
      panel.grid = element_line(color = "grey80", linewidth = 0.2)
    )
  # 8) LEGENDA
  legend_data <- bi_pal
  
  legend_plot <- ggplot(legend_data,
                        aes(x = percent_cat,
                            y = richness_cat,
                            fill = fill)) +
    geom_tile() +
    scale_fill_identity() +
    labs(
      x = "Completeness of \ngeographic distribution (%)",
      y = "Species richness"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 5, color = "black"),
      axis.text.x = element_text(size = 5, angle = 45, hjust = 0.8, color = "black"),
      axis.text.y = element_text(size = 5, color = "black"),
      panel.grid = element_blank()
    )
  # 9) combinar mapa e legenda
  final_plot <- ggdraw() +
    draw_plot(map_plot) +
    draw_plot(legend_plot,
              x = 0.085, y = 0.23,
              width = 0.18, height = 0.25)
  
  return(final_plot)
}

get_r2 <- function(model) {
  summary(model)$pseudo.r.squared
}

delta_r2 <- function(var_name, data, full_model) {
  f_full <- formula(full_model)
  f_reduced <- update(f_full, paste(". ~ . -", var_name))
  m_reduced <- betareg(f_reduced, data = data)
  
  r2_full <- get_r2(full_model)
  r2_reduced <- get_r2(m_reduced)
  
  data.frame(
    variable = var_name,
    r2_full = r2_full,
    r2_reduced = r2_reduced,
    unique_r2 = r2_full - r2_reduced
  )
}

marginal_r2 <- function(var_name, data, response) {
  f <- as.formula(paste(response, "~", var_name))
  m <- betareg(f, data = data)
  get_r2(m)
}

## Load data ----
load(file.path(
  local_directory,
  "02_data_analysis",
  "richness_completude.RData"))

load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData"))

load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_agregatte.RData")
)

load(file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_list_sensitivity.RData"
))

# Preparing data ----
# All terrestrial tetrapods
richness_all <- richness_completude_amphibia %>%
  left_join(richness_completude_reptilia, by = "name_en") %>%
  left_join(richness_completude_aves, by = "name_en") %>%
  left_join(richness_completude_mammalia, by = "name_en") %>%
  mutate(
    polygons_richness = rowSums(select(., starts_with("polygons_richness")), na.rm = TRUE),
    observed_richness = rowSums(select(., starts_with("observed_richness")), na.rm = TRUE)
  ) %>%
  select("name_en", "polygons_richness", "observed_richness") %>%
  mutate(
    percent_wallace = observed_richness * 100 / polygons_richness,
    percent_wallace_prop = percent_wallace / 100
  ) 

df_all <- left_join(
  richness_all,
  data_sociopolitic,
  by = "name_en") %>%
  left_join(geographic_shape_data %>%
              select(name_en, total_road_km) %>%
              st_drop_geometry(), by = "name_en") %>% 
  select(percent_wallace_prop, polygons_richness, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions, total_road_km, ISO3) %>%
  remove_missing() %>%
  mutate(
    n_institutions_scale = scale(log(n_institutions+1/area_km2)),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness/area_km2))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    #academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    road_dens_scaled = scale(log(as.numeric(total_road_km/area_km2))),
    #area_km2_scale = scale(log(as.numeric(area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
    #mobilization_effort_scaled = scale(log(as.numeric(mobilization_effort)))
  ) %>%
  # sensitivity analysis
  left_join(species_admin_01perc %>% 
              select(name_en, completeness_all), by = "name_en") %>% 
  rename(perc_01 = completeness_all) %>%
  left_join(species_admin_05perc %>% 
              select(name_en, completeness_all), by = "name_en") %>% 
  rename(perc_05 = completeness_all) %>%
  left_join(species_admin_1perc %>% 
              select(name_en, completeness_all), by = "name_en") %>% 
  rename(perc_1 = completeness_all) %>%
  left_join(species_admin_occ_polygons %>% 
              select(name_en, completeness_all), by = "name_en") %>% 
  rename(perc_poly = completeness_all)

vars_to_plot <- df_all %>% 
  select(
    polygons_richness_scaled,
    democracy_index_scaled,
    dens_pop_scaled,
    road_dens_scaled,
    mean_grp_scaled,
    n_institutions_scale
  ) %>%
  mutate(across(everything(), ~ as.numeric(.)))

ggpairs(data.frame(vars_to_plot))

# Beta regression ----
table(df_all$GlobalNorth)
names(df_all)

modelo_beta <- betareg(percent_wallace_prop ~ 
                         polygons_richness_scaled + 
                         dens_pop_scaled + 
                         mean_grp_scaled +
                         democracy_index_scaled +
                         colonial_origin +
                         n_institutions_scale +
                         road_dens_scaled,
                       data = df_all)
summary(modelo_beta)
car::vif(modelo_beta)

modelo_beta_occ <- betareg(perc_poly ~ 
                         polygons_richness_scaled + 
                         dens_pop_scaled + 
                         mean_grp_scaled +
                         democracy_index_scaled +
                         colonial_origin +
                         n_institutions_scale +
                         road_dens_scaled,
                       data = df_all)
summary(modelo_beta_occ)
car::vif(modelo_beta_occ)

modelo_beta_01 <- betareg(perc_01 ~ 
                             polygons_richness_scaled + 
                             dens_pop_scaled + 
                             mean_grp_scaled +
                             democracy_index_scaled +
                             colonial_origin +
                             n_institutions_scale +
                             road_dens_scaled,
                           data = df_all)
summary(modelo_beta_01)
car::vif(modelo_beta_01)

modelo_beta_05 <- betareg(perc_05 ~ 
                            polygons_richness_scaled + 
                            dens_pop_scaled + 
                            mean_grp_scaled +
                            democracy_index_scaled +
                            colonial_origin +
                            n_institutions_scale +
                            road_dens_scaled,
                          data = df_all)
summary(modelo_beta_05)
car::vif(modelo_beta_05)

modelo_beta_1 <- betareg(perc_1 ~ 
                            polygons_richness_scaled + 
                            dens_pop_scaled + 
                            mean_grp_scaled +
                            democracy_index_scaled +
                            colonial_origin +
                            n_institutions_scale +
                            road_dens_scaled,
                          data = df_all)
summary(modelo_beta_1)
car::vif(modelo_beta_1)

# TODO Pensar sobre
modelo_hier <- glmmTMB::glmmTMB(
  percent_wallace_prop ~
    polygons_richness_scaled + 
    dens_pop_scaled + 
    mean_grp_scaled +
    democracy_index_scaled +
    colonial_origin +
    n_institutions_scale +
    road_dens_scaled +
    (1 | ISO3),
  family = glmmTMB::beta_family(),
  data = df_all
)
summary(modelo_hier)

# Variation partitioning ----
vars <- c(
  "polygons_richness_scaled",
  "dens_pop_scaled",
  "mean_grp_scaled",
  "democracy_index_scaled",
  "colonial_origin",
  "n_institutions_scale",
  "road_dens_scaled"
)

results_all <- do.call(rbind, lapply(vars, function(v) {
  
  unique_part <- delta_r2(v, df_all, modelo_beta)
  total_part <- marginal_r2(v, df_all, "percent_wallace_prop")
  
  data.frame(
    variable = v,
    unique_r2 = unique_part$unique_r2,
    total_r2 = total_part,
    shared_r2 = total_part - unique_part$unique_r2
  )
})) %>% 
  mutate(model = "Species polygon cover >0% of the admin. unit area")

results_occ <- do.call(rbind, lapply(vars, function(v) {
  
  unique_part <- delta_r2(v, df_all, modelo_beta_occ)
  total_part <- marginal_r2(v, df_all, "perc_poly")
  
  data.frame(
    variable = v,
    unique_r2 = unique_part$unique_r2,
    total_r2 = total_part,
    shared_r2 = total_part - unique_part$unique_r2
  )
})) %>% 
  mutate(model = "Only species with polygons & occurences")

results_01 <- do.call(rbind, lapply(vars, function(v) {
  
  unique_part <- delta_r2(v, df_all, modelo_beta_01)
  total_part <- marginal_r2(v, df_all, "perc_01")
  
  data.frame(
    variable = v,
    unique_r2 = unique_part$unique_r2,
    total_r2 = total_part,
    shared_r2 = total_part - unique_part$unique_r2
  )
})) %>% 
  mutate(model = "Species polygon cover ≥0.1% of the admin. unit area")

results_05 <- do.call(rbind, lapply(vars, function(v) {
  
  unique_part <- delta_r2(v, df_all, modelo_beta_05)
  total_part <- marginal_r2(v, df_all, "perc_05")
  
  data.frame(
    variable = v,
    unique_r2 = unique_part$unique_r2,
    total_r2 = total_part,
    shared_r2 = total_part - unique_part$unique_r2
  )
})) %>%  
  mutate(model = "Species polygon cover ≥0.5% of the admin. unit area")

results_1 <- do.call(rbind, lapply(vars, function(v) {
  
  unique_part <- delta_r2(v, df_all, modelo_beta_1)
  total_part <- marginal_r2(v, df_all, "perc_1")
  
  data.frame(
    variable = v,
    unique_r2 = unique_part$unique_r2,
    total_r2 = total_part,
    shared_r2 = total_part - unique_part$unique_r2
  )
})) %>%  
  mutate(model = "Species polygon cover ≥1% of the admin. unit area")

results_r2 <- rbind(
  results_all,
  results_occ,
  results_01,
  results_05,
  results_1
)  %>%
  mutate(
    unique_r2 = pmax(unique_r2, 0),
    shared_r2 = pmax(shared_r2, 0))

# Forest plot ----
new_names <- tibble::tribble(
  ~term, ~term_novo,
  "polygons_richness_scaled", "Species\nrichness\nper km²",
  "dens_pop_scaled", "Human\npopulation\nper km²",
  "mean_grp_scaled", "GDP-PPP\nper capita",
  "democracy_index_scaled", "Liberal\ndemocracy\nindex",
  "colonial_origin", "Colonial\norigin",
  "n_institutions_scale", "Research\ninstitutions\nper km²",
  "road_dens_scaled", "Road\nlength\nper km²",
) #  "mobilization_effort_scaled", "Mobilization effort"
map <- setNames(new_names$term_novo, new_names$term)

# Extract avg_slopes of beta regression
coef_list <- tidy_avg_slopes(modelo_beta) %>%
  #filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(variable, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover >0% of the admin. unit area")

coef_list_occ <- tidy_avg_slopes(modelo_beta_occ) %>%
  #filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(variable, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Only species with polygons & occurences")

coef_list_01 <- tidy_avg_slopes(modelo_beta_01) %>%
  #filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(variable, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover ≥0.1% of the admin. unit area")

coef_list_05 <- tidy_avg_slopes(modelo_beta_05) %>%
  #filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(variable, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover ≥0.5% of the admin. unit area")

coef_list_1 <- tidy_avg_slopes(modelo_beta_1) %>%
  #filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(variable, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover ≥1% of the admin. unit area")

coef_all <- bind_rows(
  coef_list,
  coef_list_occ,
  coef_list_1,
  coef_list_05,
  coef_list_01
) %>%
  mutate(
    term = factor(term, levels = c(
      "Species\nrichness\nper km²",
      "Colonial\norigin",
      "Human\npopulation\nper km²",
      "Research\ninstitutions\nper km²",
      "GDP-PPP\nper capita",
      "Road\nlength\nper km²",
      "Liberal\ndemocracy\nindex"
    )),
    model = factor(model, levels = c(
      "Only species with polygons & occurences",
      "Species polygon cover >0% of the admin. unit area",
      "Species polygon cover ≥0.1% of the admin. unit area",
      "Species polygon cover ≥0.5% of the admin. unit area",
      "Species polygon cover ≥1% of the admin. unit area")
    ))

# Definir parâmetros comuns
base_size <- 8
border_size <- 0.3
tick_size <- 0.2
tick_length <- unit(2, "mm")
legend_key_size <- unit(0.25, "cm")

# PLOT 1: Coef plot
pd <- position_dodge(width = 0.8) # evitar sobreposição de barras e pontos

plot <- ggplot(coef_all, aes(x = estimate, y = term)) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high, group = model),
    height = 0.4,
    color = "black",
    position = pd,
    linewidth = 0.3  # Padronizar
  ) +
  geom_point(
    aes(shape = model, fill = significant, group = model),
    size = 1.5,
    color = "black",
    stroke = 0.5,
    position = pd
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = "gray90", "TRUE" = "black")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  labs(x = "Average marginal effects (CI 95%)", y = NULL, shape = "", fill = "Significant") +
  coord_cartesian(xlim = c(-0.15, 0.15)) +  # Ajuste conforme seus dados
  theme_minimal(base_size = base_size) +
  theme(
    legend.position = c(0.02, 0.8),
    legend.justification = c(0, 0),
    legend.text = element_text(size = base_size),
    legend.title = element_text(size = base_size),
    legend.key.size = legend_key_size,
    legend.background = element_rect(fill = "transparent", color = NA),
    
    axis.title = element_text(size = base_size),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(size = base_size, color = "black"),
    axis.text.x = element_text(size = base_size, color = "black"),
    
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = border_size),
    
    plot.margin = margin(5, 2, 5, 5),
    axis.ticks = element_line(color = "black", linewidth = tick_size),
    axis.ticks.length = tick_length,
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    axis.line = element_blank()  # Remove linha do eixo (usamos panel.border)
  ) + guides(
    shape = guide_legend(
      override.aes = list(fill = "black"),
      reverse = TRUE  # <--- INVERTE A ORDEM DA LEGENDA
    ),
    fill = "none"
  ); plot

# PLOT 2: R²
df_long <- results_r2 %>%
  pivot_longer(cols = c(unique_r2, shared_r2),
               names_to = "type",
               values_to = "value") %>%
  mutate(type = dplyr::recode(type,
                              unique_r2 = "Unique Pseudo-R²",
                              shared_r2 = "Shared Pseudo-R²"),
         variable = case_when(
           variable == "democracy_index_scaled" ~ "Liberal democracy index",
           variable == "colonial_origin"        ~ "Colonial origin",
           variable == "mean_grp_scaled"        ~ "GDP-PPP per capita",
           variable == "dens_pop_scaled"        ~ "Human population per km²",
           variable == "road_dens_scaled"       ~ "Road length per km²",
           variable == "n_institutions_scale"   ~ "Research institutions per km²",
           variable == "polygons_richness_scaled" ~ "Species richness per km²",
           TRUE ~ NA_character_
         ),
         variable = factor(variable, levels = c(
           "Species richness per km²",
           "Colonial origin",
           "Human population per km²",
           "Research institutions per km²",
           "GDP-PPP per capita",
           "Road length per km²",
           "Liberal democracy index"
         )),
         model = factor(model, levels = c(
           "Only species with polygons & occurences",
           "Species polygon cover >0% of the admin. unit area",
           "Species polygon cover ≥0.1% of the admin. unit area",
           "Species polygon cover ≥0.5% of the admin. unit area",
           "Species polygon cover ≥1% of the admin. unit area")
         ))

df_long <- df_long %>%
  arrange(variable, model, desc(type))
df_long$type <- factor(df_long$type, levels = c("Unique Pseudo-R²",
                                                "Shared Pseudo-R²"))

r2_plot <- ggplot(df_long, aes(x = variable, y = value, fill = type)) +
  geom_col(
    aes(group = model),
    position = position_dodge(preserve = "single", width = 0.8),
    width = 0.8,
    color = "white",
    linewidth = 0.1
  ) +
  #geom_text(aes(label = model, group = model),
  #          position = position_dodge(width = 0.7),
  #          hjust = -0.1,  # posição à direita da barra
  #          size = 2,
  #          angle = 0) +
  scale_fill_manual(values = c("Unique Pseudo-R²" = "#858F94", "Shared Pseudo-R²" = "#bdbdbd")) +
  labs(x = NULL, y = "Pseudo R²", fill = NULL) +
  scale_y_continuous(
    expand = c(0, 0), 
    limits = c(0, 0.4)
  ) +
  coord_flip() +
  theme_classic(base_size = base_size) +
  theme(
    legend.position = c(0.98, 0),
    legend.justification = c(1, 0),
    legend.key.size = legend_key_size,
    legend.key.width = legend_key_size,
    legend.key.height = legend_key_size,
    legend.text = element_text(size = base_size),
    legend.background = element_rect(fill = "transparent", color = NA),
    
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = base_size, color = "black"),
    axis.title.x = element_text(size = base_size, margin = margin(t = 10)),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),  # <-- AQUI!
    #axis.line = element_line(color = "black", linewidth = tick_size),
    axis.line = element_blank(),
    #panel.border = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = border_size),
    axis.ticks = element_line(color = "black", linewidth = tick_size),
    axis.ticks.length = tick_length,
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(),
    
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(5, 5, 5, 2)
  ); r2_plot

# Unindo os plots
coefplot_r2 <- plot + r2_plot + 
  plot_layout(widths = c(1.2, 0.9)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  ); coefplot_r2

# https://www.nature.com/documents/natrev-figure-guidelines-v1.pdf
# largura maxima é 180mm
# altura maxima é 215mm
# fonte size, de 8pt
ggsave(filename="Figure2_EffectSize.pdf", 
       plot= coefplot_r2,
       width=180, 
       height=160, 
       units="mm",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf) 

# Plot compare North vs South----
# >0%
wtest_0 <- wilcox.test(percent_wallace_prop ~ GlobalNorth, data = df_all)
p_label_0 <- paste0(
  "W = ", round(wtest_0$statistic, 1)," ","p = ", format.pval(wtest_0$p.value, digits = 2, eps = .001)
)
# Polygons&Occ
wtest_occ <- wilcox.test(perc_poly ~ GlobalNorth, data = df_all)
p_label_occ <- paste0(
  "W = ", round(wtest_occ$statistic, 1)," ","p = ", format.pval(wtest_occ$p.value, digits = 2, eps = .001)
)
# 0.1%
wtest_01 <- wilcox.test(perc_01 ~ GlobalNorth, data = df_all)
p_label_01 <- paste0(
  "W = ", round(wtest_01$statistic, 1), " ","p = ", format.pval(wtest_01$p.value, digits = 2, eps = .001)
)
# 0.5%
wtest_05 <- wilcox.test(perc_05 ~ GlobalNorth, data = df_all)
p_label_05 <- paste0(
  "W = ", round(wtest_05$statistic, 1), " ","p = ", format.pval(wtest_05$p.value, digits = 2, eps = .001)
)
# 1%
wtest_1 <- wilcox.test(perc_1 ~ GlobalNorth, data = df_all)
p_label_1 <- paste0(
  "W = ", round(wtest_1$statistic, 1), " ","p = ", format.pval(wtest_1$p.value, digits = 2, eps = .001)
)

# ordem desejada
ordem_threshold <- c(
  "Only species with \npolygons & occurences",
  "Species polygon\ncover >0% of the\nadmin. unit area",
  "Species polygon \ncover ≥0.1% of the \nadmin. unit area",
  "Species polygon \ncover ≥0.5% of the \nadmin. unit area",
  "Species polygon \ncover ≥1% of the \nadmin. unit area"
)

df_all_pivot <- df_all %>%
  mutate(regiao = if_else(GlobalNorth == 1, "North", "South"),
         regiao = as.factor(regiao)) %>%
  select("GlobalNorth", "regiao", "percent_wallace_prop","perc_poly",
         "perc_01","perc_05","perc_1", "name_en") %>%
  pivot_longer(
    cols = c(percent_wallace_prop, perc_poly, perc_01, perc_05, perc_1),
    names_to = "threshold",
    values_to = "valor"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "percent_wallace_prop" ~ "Species polygon\ncover >0% of the\nadmin. unit area",
      threshold == "perc_poly" ~ "Only species with \npolygons & occurences",
      threshold == "perc_01" ~ "Species polygon \ncover ≥0.1% of the \nadmin. unit area",
      threshold == "perc_05" ~ "Species polygon \ncover ≥0.5% of the \nadmin. unit area",
      threshold == "perc_1"  ~ "Species polygon \ncover ≥1% of the \nadmin. unit area"
    ),
    threshold = factor(threshold, levels = ordem_threshold)
  ) %>%
  select(GlobalNorth, regiao, name_en, threshold, valor)

boxplot <- ggplot(data = df_all_pivot, 
                  aes(y = valor, x = threshold, fill = regiao)) +
  
  # density slab (mais fino)
  stat_slab(
    side = "right", 
    scale = 0.55,  # tamanho dos histogramas
    position = position_dodge(width = .65), # distancia dos histogramas
    aes(fill_ramp = stat(level)),
    .width = c(.50, .95, 1)
  ) +
  
  # dots (menores)
  stat_dots(
    side = "left", 
    scale = 0.55,  # chuva de pontos
    #binwidth = 0.05,   # menor = mais pontos / mais compacto
    #stackratio = 0.05,    # controla compactação vertical
    show.legend = TRUE,
    position = position_dodge(width = .6),
    aes(color = regiao)
  ) +
  
  # whisker (mais fino)
  stat_summary(
    fun = median,
    fun.min = function(z) median(z) - sd(z),
    fun.max = function(z) median(z) + sd(z),
    geom = "errorbar",
    width = 0.12,     # tamanho das barras do error bar
    color = "black",
    size = 0.32,       #  expessura da linha do error bar
    position = position_dodge(width = .6)
  ) +
  
  # ponto mediana (menor)
  stat_summary(
    aes(shape = threshold, group = interaction(regiao, threshold)),
    fun = median,
    geom = "point",
    color = "black",
    fill = "black",
    size = 0.8,       # tamanho do simbolo
    stroke = 0.8,     # expessura do simbolo
    position = position_dodge(width = .6)
  ) +
  
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  
  scale_fill_manual(values = c("North" = "#2DA5E8", "South"="#c40e3e")) +
  scale_color_manual(values = c("North" = "#2DA5E8", "South"="#c40e3e")) +
  
  guides(fill_ramp = "none") +
  
  labs(
    x = "",
    y = "Geographic distribution completeness (%)",
    fill = ""
  ) +
  
  coord_flip() +
  
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.3),
    
    axis.text = element_text(color = "black", size = 7),  # 🔧 menor
    axis.title.x = element_text(size = 7),
    
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    
    legend.position = "none",
    
    plot.margin = margin(10, 10, 10, 10) # 🔧 menos margem
  ); boxplot

global_plot <- ggplot(geographic_shape_data) +
  geom_sf(aes(fill = as.factor(GlobalNorth)), 
          color = "black", size = 0.05) +  # 🔧 linha um pouco mais visível
  
  scale_fill_manual(
    values = c("1" = "#2DA5E8",
               "0" = "#c40e3e")
  ) +
  
  coord_sf(datum = NA) +  # 🔧 importante para proporção correta
  
  theme_void() +
  theme(
    legend.position = "none"
  )

# plot 
cairo_pdf("figure.pdf", 
          width = 7.09, 
          height = 4.73, 
          pointsize = 8,
          family = "Arial")
par(
  lwd = 0.5,
  xaxs = "r",
  yaxs = "r"
)
plot(tetrapods_100,
     se = TRUE,
     legend = FALSE,
     ylim = c(0, 30),
     xlim = c(1740, 2030),
     xlab = "Year",
     ylab = "Marginal mean",
     col = c("#c40e3e", "#2DA5E8"),
     xaxt = "n",
     yaxt = "n"   # 🔧 controla eixo Y também
)

# eixo X fino
axis(1,
     at = anos_25,
     labels = ifelse(anos_25 %in% anos_50, anos_25, ""),
     las = 1,
     lwd = 0.3,
     lwd.ticks = 0.3
)

# eixo Y fino
axis(2,
     lwd = 0.3,
     lwd.ticks = 0.3)

# linhas verticais a cada 25 anos
abline(v = anos_25,
       col = "gray70",
       lty = 3,
       lwd = 0.4)

# mapa menor (menos achatado)
vp <- grid::viewport(
  x = 0.70, y = 0.74,
  width = 0.40,   # 🔧 antes 0.45
  height = 0.55   # 🔧 antes 0.60
)
print(global_plot, vp = vp)

# boxplot menor e mais compacto
vp2 <- grid::viewport(
  x = 0.28, y = 0.58,
  width = 0.45,   # 🔧 antes 0.45
  height = 0.60  # 🔧 antes 0.65
)
print(boxplot, vp = vp2)
dev.off()

#Largura 7.09  - Altura final	4.73

# export 12x8 pdf
# largura maxima 180mm - Width 7.09 in 
# altura maxima 215 - Heigh 8.47 in
# font size = 8 mm

# Plot completeness map ----
plot <- plot_bivariate_map(
  geographic_shape_data,
  richness_all,
  title = ""); plot

ggsave(filename= "Figure2_TerrestrialVertebrates.pdf",
       plot= plot,
       width=7.09,
       height=4.73,
       units="in",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf)

# Countries and democracy index ----
# Define colors for each biogeographic realm (Pastel1 from RColorBrewer):
summary(df_all$percent_wallace_prop)
summary(df_all$democracy_index_mean)

# abreviacao das und. adm
states_ref <- rnaturalearth::ne_states(returnclass = "sf") %>%
  st_drop_geometry() %>%
  filter(adm0_a3 %in% c("BRA", "RUS", "CHN", "USA", "CAN", "AUS")) %>%
  select(name_en, iso_3166_2) 

df_all_abs <- df_all %>%
  left_join(states_ref, by = "name_en") %>%
  mutate(
    ID = coalesce(iso_3166_2, ISO3)
  ) %>%
  mutate(
    ID = ifelse(
      grepl("-", ID),
      paste0(ISO3, "-", sub("^[A-Z]{2}-", "", ID)),
      ISO3
    )
  )

plot <- ggplot(df_all_abs, aes(x = democracy_index_mean,
                           y = percent_wallace_prop)) +
  
  geom_point(aes(color = as.factor(GlobalNorth)),
             size = 2,
             alpha = 0.7) +
  
  geom_smooth(method = "lm",
              se = TRUE,
              level = 0.95,
              color = "black",
              linetype = "solid",
              alpha = 0.2) +
  
  geom_label_repel(
    aes(label = ID,
        color = as.factor(GlobalNorth)),
    size = 2,
    max.overlaps = 45,
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.2,
    min.segment.length = 0.1,
    force = 1,
    fill = "white",      # fundo branco
    label.size = 0.2     # borda fina
  ) +
  
  scale_color_manual(
    name = "",
    values = c("0" = "#c40e3e", "1" = "#2DA5E8"),
    labels = c("0" = "Global South", "1" = "Global North")
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # começa em 0, máximo automático
    expand = expansion(mult = c(0, 0.05))  # remove espaço extra na base
  ) +
  
  labs(
    x = "Liberal Democracy index",
    y = "Geographic knowledge completeness (%)"
  ) +
  
  theme_classic() +
  theme(
    legend.direction = "horizontal", 
    legend.position = c(1, 0.01), 
    legend.justification = c(1, 0),    
    legend.background = element_blank(),
    axis.title.x = element_text(size = base_size),
    axis.title.y = element_text(size = base_size),
    axis.text.x  = element_text(size = base_size),
    axis.text.y  = element_text(size = base_size),
    axis.ticks = element_line(color = "black", linewidth = tick_size),
    axis.ticks.length = tick_length,
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    axis.line = element_blank(),  # Remove linha do eixo (usamos panel.border)
    legend.text = element_text(size = 8),
    panel.border = element_rect(color = "black",
                                fill = NA,
                                linewidth = border_size)
  ); plot

ggsave(filename= "FigS3CompletenessLDI.pdf", 
       plot= plot,
       width=180, 
       height=130, 
       units="mm",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf) 

# Sensitivity analysis ----
load(file.path(
  local_directory,
  "00_raw_data",
  "tetrapods_list_sensitivity.RData"
))

# Boxplot comparing units between north and south
sp_north <- df_wallacean_100 %>%
  filter(GlobalNorth == 1) %>%
  distinct(speciesKey)
nrow(sp_north)

sp_south <- df_wallacean_100 %>%
  filter(GlobalNorth == 0) %>%
  distinct(speciesKey)
nrow(sp_south)

only_south <- sp_south %>%
  anti_join(sp_north, by = "speciesKey") %>%
  left_join(sp_count, by = "speciesKey") %>%
  mutate(only_north = 0)
nrow(only_south) # only south 10912

only_north <- sp_north %>%
  anti_join(sp_south, by = "speciesKey") %>%
  left_join(sp_count, by = "speciesKey") %>%
  mutate(only_north = 1)
nrow(only_north) # only norte 1414

both <- sp_north %>%
  inner_join(sp_south, by = "speciesKey")
nrow(both) # 2944 com occ no norte e no sul

df <- bind_rows(
  only_south,
  only_north
)

hist(df$ua_total)

plot <- ggplot(df, aes(
  x = factor(only_north, labels = c("Global South", "Global North")),
  y = ua_total,
  fill = factor(only_north)
)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "N. of occupied \nadministrative units"
  ) +
  scale_fill_manual(
    values = c("0" = "#c40e3e", "1" = "#2DA5E8"),
    labels = c("0" = "Global South", "1" = "Global North")
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    legend.position = "none",              # remove legenda
    axis.text.x = element_text(size = 14), # aumenta fonte eixo X
    axis.text.y = element_text(size = 14), # aumenta fonte eixo Y
    axis.title.y = element_text(size = 16) # aumenta título do eixo Y
  ); plot

ggsave(file.path(drive, "FigS4.tiff"), 
       plot = plot, 
       width = 8, 
       height = 7, 
       dpi = 300, 
       bg = "white")

