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
  map_plot <- ggplot(dados_bi) +
    geom_sf(aes(fill = fill), color = "black", size = 0.15) +
    scale_fill_identity() +
    theme_minimal() +
    labs(title = title) +
    theme(legend.position = "none")
  
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
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )
  
  # 9) combinar mapa e legenda
  final_plot <- ggdraw() +
    draw_plot(map_plot) +
    draw_plot(legend_plot,
              x = 0.05, y = 0.22,
              width = 0.18, height = 0.25)
  
  return(final_plot)
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
         colonial_origin, n_institutions, total_road_km) %>%
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

# Amphibia
#richness_amphibia <- richness_completude_amphibia %>%
#  select("name_en", "polygons_richness", "observed_richness") %>%
#  mutate(
#    percent_wallace = observed_richness * 100 / polygons_richness,
#    percent_wallace_prop = percent_wallace / 100
#  ) 
#
#df_amphibia <- left_join(
#  richness_amphibia,
#  data_sociopolitic,
#  by = "name_en") %>%
#  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
#         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
#         colonial_origin, n_institutions) %>%
#  mutate(
#    n_institutions_scale = scale(n_institutions),
#    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
#    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
#    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
#    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
#    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
#  ) %>% 
#  remove_missing()
#
## Reptilia
#richness_reptilia <- richness_completude_reptilia %>%
#  select("name_en", "polygons_richness", "observed_richness") %>%
#  mutate(
#    percent_wallace = observed_richness * 100 / polygons_richness,
#    percent_wallace_prop = percent_wallace / 100
#  ) 
#
#df_reptilia <- left_join(
#  richness_reptilia,
#  data_sociopolitic,
#  by = "name_en") %>%
#  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
#         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
#         colonial_origin, n_institutions) %>%
#  mutate(
#    n_institutions_scale = scale(n_institutions),
#    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
#    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
#    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
#    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
#    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
#  ) %>% 
#  remove_missing()
#
## Aves
#richness_aves <- richness_completude_aves %>%
#  select("name_en", "polygons_richness", "observed_richness") %>%
#  mutate(
#    percent_wallace = observed_richness * 100 / polygons_richness,
#    percent_wallace_prop = percent_wallace / 100
#  ) 
#
#df_aves <- left_join(
#  richness_aves,
#  data_sociopolitic,
#  by = "name_en") %>%
#  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
#         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
#         colonial_origin, n_institutions) %>%
#  mutate(
#    n_institutions_scale = scale(n_institutions),
#    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
#    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
#    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
#    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
#    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
#  ) %>% 
#  remove_missing()
#
## Mammalia
#richness_mammalia <- richness_completude_mammalia %>%
#  select("name_en", "polygons_richness", "observed_richness") %>%
#  mutate(
#    percent_wallace = observed_richness * 100 / polygons_richness,
#    percent_wallace_prop = percent_wallace / 100
#  ) 
#
#df_mammalia <- left_join(
#  richness_mammalia,
#  data_sociopolitic,
#  by = "name_en") %>%
#  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
#         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
#         colonial_origin, n_institutions) %>%
#  mutate(
#    n_institutions_scale = scale(n_institutions),
#    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
#    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
#    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
#    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
#    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
#  ) %>% 
#  remove_missing()

# Plot maps ----
drive <- file.path("G:", "Meu Drive", "Artigos", "wallace_time")

plot <- plot_bivariate_map(geographic_shape_data,
                   richness_all,
                   title = "Terrestrial vertebrates"); plot

ggsave(filename= 
         file.path(drive,
                   "Figures",
                   "Figure2_TerrestrialVertebrates.pdf"),
       plot= plot,
       width=12,
       height=8,
       units="in",
       bg="white",
       limitsize=F)

# Class
plot_amphibia <- plot_bivariate_map(geographic_shape_data,
                                    richness_completude_amphibia,
                                    title = "Amphibia")
plot_reptilia <- plot_bivariate_map(geographic_shape_data,
                                    richness_completude_reptilia,
                                    title = "Reptilia")
plot_aves <- plot_bivariate_map(geographic_shape_data,
                                richness_completude_aves,
                                title = "Reptilia")
plot_mammalia <- plot_bivariate_map(geographic_shape_data,
                                    richness_completude_mammalia,
                                    title = "Mammalia")

ggsave(filename="Figures/Figure1_Amphibia.pdf",
       plot= plot_amphibia,
       width=12,
       height=8,
       units="in",
       bg="white",
       limitsize=F)

ggsave(filename="Figures/Figure1_Reptilia.pdf",
       plot= plot_reptilia,
       width=12,
       height=8,
       units="in",
       bg="white",
       limitsize=F)

ggsave(filename="Figures/Figure1_Aves.pdf",
       plot= plot_aves,
       width=12,
       height=8,
       units="in",
       bg="white",
       limitsize=F)

ggsave(filename="Figures/Figure1_Mammalia.pdf",
       plot= plot_mammalia,
       width=12,
       height=8,
       units="in",
       bg="white",
       limitsize=F)

# Beta regression ----
table(df_all$GlobalNorth)
table(df_amphibia$GlobalNorth)
table(df_reptilia$GlobalNorth)
table(df_aves$GlobalNorth)
table(df_mammalia$GlobalNorth)

names(df_all)
names(species_admin_05perc)

# Ajustar beta regression (percent_wallace já em proporção 0-1)
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

# Forest plot ----
new_names <- tibble::tribble(
  ~term, ~term_novo,
  "polygons_richness_scaled", "Species\nrichness density",
  "dens_pop_scaled", "Population density",
  "mean_grp_scaled", "GDP per capta",
  "democracy_index_scaled", "Liberal \ndemocracy index",
  "colonial_origin", "Colonial origin",
  "n_institutions_scale", "Research\ninstitutions density",
  "road_dens_scaled", "Accessibility\ndensity"
) #  "mobilization_effort_scaled", "Mobilization effort"
map <- setNames(new_names$term_novo, new_names$term)

coef_list <- tidy(modelo_beta) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover >0% of the admin. unit area")

coef_list_occ <- tidy(modelo_beta_occ) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Only species with polygons & occurences")

coef_list_01 <- tidy(modelo_beta_01) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover ≥0.1% of the admin. unit area")

coef_list_05 <- tidy(modelo_beta_05) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)),
         model = "Species polygon cover ≥0.5% of the admin. unit area")

coef_list_1 <- tidy(modelo_beta_1) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map)) %>%
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
    term = factor(term, levels = rev(c(
      "Species\nrichness density",
      "Colonial origin",
      "Population density",
      "Research\ninstitutions density",
      "GDP per capta",
      "Accessibility\ndensity",
      "Liberal \ndemocracy index"
    ))),
    model = factor(model, levels = rev(c(
      "Only species with polygons & occurences",
      "Species polygon cover >0% of the admin. unit area",
      "Species polygon cover ≥0.1% of the admin. unit area",
      "Species polygon cover ≥0.5% of the admin. unit area",
      "Species polygon cover ≥1% of the admin. unit area")
    )))

# Para evitar sobreposição de barras e pontos
pd <- position_dodge(width = 0.8)

plot <- ggplot(coef_all, aes(x = estimate, y = term)) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error,
        group = model),
    height = 0.4,
    color = "black",
    position = pd,
    size = 0.7
  ) +
  geom_point(
    aes(shape = model, fill = significant, group = model),
    size = 2,
    color = "black",
    stroke = 1,
    position = pd
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  scale_fill_manual(values = c("FALSE" = "gray90", "TRUE" = "black")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  labs(
    x = "Model coefficients (CI 95%)",
    y = "",
    shape = "Model",
    fill = "Significant"
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.05, 0.08),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 12),     # 👈 aumenta aqui
    legend.title = element_text(size = 13),
    #legend.background = element_rect(fill = "white"),
    #legend.key = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(margin = margin(t = 14)),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) + guides(
    shape = guide_legend(
      override.aes = list(fill = "black")
    ),
    fill = "none"
  ); plot

# Compare North vs South Completeness (Boxplot) ----
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

# pivot data.frame
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
  
  # density slab
  stat_slab(
    side = "right", scale = 0.4,
    position = position_dodge(width = .8),
    aes(fill_ramp = stat(level)),
    .width = c(.50, .95, 1)
  ) +
  # dots
  stat_dots(
    side = "left", scale = 0.4, show.legend = TRUE,
    position = position_dodge(width = .8),
    aes(color = regiao)
  ) +
  # whisker preto (mediana ± SD)
  stat_summary(
    fun = median,
    fun.min = function(z) median(z) - sd(z),
    fun.max = function(z) median(z) + sd(z),
    geom = "errorbar",
    width = 0.2,
    color = "black",
    size = 0.7,
    position = position_dodge(width = .8)
  ) +
  # ponto preto da mediana
  stat_summary(
    aes(shape = threshold, group = interaction(regiao, threshold)),
    fun = median,
    geom = "point",
    color = "black",
    fill = "black",
    size = 2,
    stroke = 1,
    position = position_dodge(width = .8)
  ) +
  scale_shape_manual(values = rev(c(21, 22, 23, 24, 25))) +
  # styling
  scale_fill_manual(values = c("North" = "#2DA5E8", "South"="#c40e3e")) +
  scale_color_manual(values = c("North" = "#2DA5E8", "South"="#c40e3e")) +
  guides(fill_ramp = "none") +
  # Rótulos dos eixos e título
  labs(
    x = "",
    y = "Geographic distribution completeness (%)",
    title = "",
    fill = ""  # Título vazio para legenda de fill
  ) +
  coord_flip() +
  # Tema personalizado
  theme(
    panel.background = element_rect(
      fill = "white",
    ),
    plot.background = element_rect(
      fill = "white",
    ),
    #panel.background = element_blank(),
    #plot.background = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(size = 12, margin = margin(t = 0)),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, margin = margin(r = 10)),
    axis.ticks.y = element_line(color = "black", size = 0.5),   # aqui
    axis.ticks.x = element_line(color = "black", size = 0.5),   # opcional
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  ) + # export 6.5 x 6
  annotate(
    "text",
    x = 2.5,                 
    y = 0.9,
    label = "***",
    hjust = 0.5,             
    vjust = -20.5,           
    size = 3.5,
    fontface = "italic"
  ) +
  annotate(
    "text",
    x = 2.5,                 
    y = 0.9,
    label = "***",
    hjust = 0.5,             
    vjust = -13.5,           
    size = 3.5,
    fontface = "italic"
  ) +
  annotate(
    "text",
    x = 2.5,                 
    y = 0.9,
    label = "***",
    hjust = 0.5,             
    vjust = -6.0,            
    size = 3.5,
    fontface = "italic"
  ) +
  annotate(
    "text",
    x = 2.5,                 
    y = 0.9,
    label = "***",
    hjust = 0.5,             
    vjust = 1.5,             
    size = 3.5,
    fontface = "italic"
  ) +
  annotate(
    "text",
    x = 2.5,                 
    y = 0.9,
    label = "***",
    hjust = 0.5,             
    vjust = 8.7,             
    size = 3.5,
    fontface = "italic"); boxplot

# Plotando as curvas de acumulação + NorthSouth Global map
load(
  file.path(
    local_directory, "02_data_analysis", "tetrapods_model.RData")
)

global_plot <- ggplot(geographic_shape_data) +
  geom_sf(aes(fill = as.factor(GlobalNorth)), color = "black", size = 0.01) +
  scale_fill_manual(
    values = c("1" = "#2DA5E8",   # North
               "0" = "#c40e3e")   # South
  ) +
  theme_void() +  # remove tudo: eixos, ticks, grid, fundo
  theme(
    legend.position = "none",
    panel.background = element_blank(),   # fundo transparente
    plot.background  = element_blank()    # fundo transparente
  ); global_plot

dev.off()
# sequência de 25 em 25
anos_25 <- seq(1750, 2025, by = 25)
# sequência de 50 em 50 (labels)
anos_50 <- seq(1750, 2025, by = 50)
par(cex.main = 2,
    cex.lab  = 1.5,
    cex.axis = 1.3)

plot(tetrapods_100,
     se = TRUE,
     legend = FALSE,
     ylim = c(0, 30),
     xlim = c(1750, 2025),
     xlab = "Year",
     ylab = "Marginal mean",
     col = c("#c40e3e", "#2DA5E8"),
     xaxt = "n"
)
# ticks a cada 25 anos, labels só nos anos de 50
axis(1,
     at = anos_25,
     labels = ifelse(anos_25 %in% anos_50, anos_25, ""),  # só mostra labels nos anos de 50
     las = 1)
# linhas verticais a cada 25 anos
abline(v = anos_25,
       col = "gray70",
       lty = 3,
       lwd = 0.5)

# sobrepor mapa
vp <- grid::viewport(x = 0.70, y = 0.75, width = 0.45, height = 0.60)
print(global_plot, vp = vp)
# sobrepor boxplot
vp2 <- grid::viewport(x = 0.27, y = 0.62, width = 0.45, height = 0.65)
print(boxplot, vp = vp2)
# export 12x8 pdf

# Countries and democracy index ----
# Define colors for each biogeographic realm (Pastel1 from RColorBrewer):
myColors<-c("#2171b5",  "#08306b", "#4292c6", "#6baed6","#fc9272", "#fb6a4a", "#ef3b2c", "#a50f15", "#67000d")
names(myColors)<-c("Australia and New Zealand",
                   "Europe",
                   "Northern America",
                   "Russian Federation",
                   "Central, East and South Asia",
                   "Latin America and Caribbean", 
                   "Near East and Northern Africa",
                   "Sub-Saharan Africa",
                   "Southeast Asia and Pacific Islands")

table(df_all$GlobalNorth)

summary(df_all$percent_wallace_prop)
summary(df_all$democracy_index_mean)

# Criar o gráfico base
plot <- ggplot(df_all, aes(x = democracy_index_mean, y = percent_wallace_prop)) +
  
  geom_point(aes(color = as.factor(GlobalNorth)),
             size = 2,
             alpha = 0.7) +
  
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black",
              linetype = "solid",
              alpha = 0.2) +
  
  geom_text_repel(aes(label = name_en,
                      color = as.factor(GlobalNorth)),
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.2,
                  min.segment.length = 0.1,
                  force = 1) +
  
  scale_color_manual(
    name = "",
    values = c("0" = "#c40e3e", "1" = "#2DA5E8"),
    labels = c("0" = "Global south", "1" = "Global north")
  ) +
  
  labs(
    x = "Liberal Democracy index",
    y = "Completeness (%)"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom",
    
    # Aumentar títulos dos eixos
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    # Aumentar labels dos eixos
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    
    # Aumentar legenda
    legend.text = element_text(size = 14)
  ); plot

ggsave(file.path(drive, "FigS3.tiff"), 
       plot = plot, 
       width = 12, 
       height = 8, 
       dpi = 300, 
       bg = "white")

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

