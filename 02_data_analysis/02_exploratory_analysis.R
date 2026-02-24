library(tidyverse)
library(sf)
library(stringi)
library(betareg)
library(cowplot)
library(ggstats)
library(ggrepel)
library(ggpubr)
library(car)
library(broom)
library(tidyr)

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
  select(percent_wallace_prop, polygons_richness, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>%
  remove_missing() %>%
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
    #mobilization_effort_scaled = scale(log(as.numeric(mobilization_effort)))
  ) 

# Amphibia
richness_amphibia <- richness_completude_amphibia %>%
  select("name_en", "polygons_richness", "observed_richness") %>%
  mutate(
    percent_wallace = observed_richness * 100 / polygons_richness,
    percent_wallace_prop = percent_wallace / 100
  ) 

df_amphibia <- left_join(
  richness_amphibia,
  data_sociopolitic,
  by = "name_en") %>%
  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>%
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  ) %>% 
  remove_missing()

# Reptilia
richness_reptilia <- richness_completude_reptilia %>%
  select("name_en", "polygons_richness", "observed_richness") %>%
  mutate(
    percent_wallace = observed_richness * 100 / polygons_richness,
    percent_wallace_prop = percent_wallace / 100
  ) 

df_reptilia <- left_join(
  richness_reptilia,
  data_sociopolitic,
  by = "name_en") %>%
  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>%
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  ) %>% 
  remove_missing()

# Aves
richness_aves <- richness_completude_aves %>%
  select("name_en", "polygons_richness", "observed_richness") %>%
  mutate(
    percent_wallace = observed_richness * 100 / polygons_richness,
    percent_wallace_prop = percent_wallace / 100
  ) 

df_aves <- left_join(
  richness_aves,
  data_sociopolitic,
  by = "name_en") %>%
  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>%
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  ) %>% 
  remove_missing()

# Mammalia
richness_mammalia <- richness_completude_mammalia %>%
  select("name_en", "polygons_richness", "observed_richness") %>%
  mutate(
    percent_wallace = observed_richness * 100 / polygons_richness,
    percent_wallace_prop = percent_wallace / 100
  ) 

df_mammalia <- left_join(
  richness_mammalia,
  data_sociopolitic,
  by = "name_en") %>%
  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>%
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  ) %>% 
  remove_missing()

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
                         n_institutions_scale,
                       data = df_all)
summary(modelo_beta)
car::vif(modelo_beta)

modelo_beta_amphibia <- betareg(percent_wallace_prop ~ 
                         polygons_richness_scaled + 
                         dens_pop_scaled + 
                         mean_grp_scaled +
                         academic_freedom_scaled +
                         colonial_origin +
                         n_institutions_scale,
                       data = df_amphibia)
summary(modelo_beta_amphibia)
car::vif(modelo_beta_amphibia)

modelo_beta_reptilia <- betareg(percent_wallace_prop ~ 
                                  polygons_richness_scaled + 
                                  dens_pop_scaled + 
                                  mean_grp_scaled +
                                  academic_freedom_scaled +
                                  colonial_origin +
                                  n_institutions_scale,
                                data = df_reptilia)
summary(modelo_beta_reptilia)
car::vif(modelo_beta_reptilia)

modelo_beta_aves <- betareg(percent_wallace_prop ~ 
                                  polygons_richness_scaled + 
                                  dens_pop_scaled + 
                                  mean_grp_scaled +
                                  academic_freedom_scaled+
                                  colonial_origin +
                                  n_institutions_scale,
                                data = df_aves)
summary(modelo_beta_aves)
car::vif(modelo_beta_aves)

modelo_beta_mammalia <- betareg(percent_wallace_prop ~ 
                              polygons_richness_scaled + 
                              dens_pop_scaled + 
                              mean_grp_scaled +
                              academic_freedom_scaled +
                              colonial_origin +
                              n_institutions_scale,
                            data = df_mammalia)
summary(modelo_beta_mammalia)
car::vif(modelo_beta_mammalia)

# Forest plot ----
new_names <- tibble::tribble(
  ~term, ~term_novo,
  "polygons_richness_scaled", "Species richness",
  "dens_pop_scaled", "Population density",
  "mean_grp_scaled", "GDP per capta",
  "democracy_index_scaled", "Liberal \ndemocracy index",
  "colonial_origin", "Colonial origin",
  "n_institutions_scale", "N. of research \ninstitutions"
) #  "mobilization_effort_scaled", "Mobilization effort"
map <- setNames(new_names$term_novo, new_names$term)

coef_list <- tidy(modelo_beta) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map)) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)))

ggplot(coef_list, aes(x = estimate, y = reorder(term, estimate))) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error, 
        xmax = estimate + 1.96 * std.error),
    height = 0.1, color = "black", linewidth = 1
  ) +
  geom_point(aes(fill = significant),
             size = 3,
             shape = 21,
             color = "black",
             stroke = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50",
             linewidth = 0.5) +
  scale_fill_manual(
    values = c("FALSE" = "gray90", "TRUE" = "black")
  ) +
  labs(
    x = "Model coefficients (CI 95%)",
    y = ""
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

model_list <- list(
  "Terrestrial tetrapods" = modelo_beta,
  "Amphibia" = modelo_beta_amphibia,
  "Reptilia" = modelo_beta_reptilia,
  "Aves" = modelo_beta_aves,
  "Mammalia" = modelo_beta_mammalia
)

coef_all <- bind_rows(
  lapply(names(model_list), function(model_name) {
    tidy(model_list[[model_name]], conf.int = TRUE) %>%
      mutate(modelo = model_name)
  })
) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))

# Criar plot combinado
ggplot(coef_all, aes(x = estimate, y = term, color = modelo)) +
  # Linha de referência
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  
  # Intervalos de confiança
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.1,
    position = position_dodge(width = 0.7),
    size = 0.8
  ) +
  
  # Pontos
  geom_point(
    size = 3,
    position = position_dodge(width = 0.7)
  ) +

  # Cores por grupo
  scale_color_brewer(palette = "Set1", name = "Grupo") +
  
  # Tema
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "bottom",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  
  labs(
    title = "",
    subtitle = "",
    x = "Model coefficients (CI 95%)",
    y = ""
  )

# Compare North vs South Completeness (Boxplot) ----
wtest <- wilcox.test(percent_wallace_prop ~ GlobalNorth, data = df_all)
p_label <- paste0(
  "Wilcoxon test\n",
  "W = ", round(wtest$statistic, 1),
  "\np = ", format.pval(wtest$p.value, digits = 2, eps = .001)
)

boxplot <- ggplot(df_all, aes(x = factor(GlobalNorth),
                   y = percent_wallace_prop)) +
  # Fundo branco primeiro
  theme_minimal() +
  # Boxplot com fill
  geom_boxplot(
    aes(fill = as.factor(GlobalNorth)),
    outlier.shape = NA  # Remover outliers padrão
  ) +
  # Adicionar pontos individuais (jitter) bem sutis
  geom_jitter(
    aes(color = as.factor(GlobalNorth)),
    width = 0.2,        # Pouco espalhamento horizontal
    height = 0,         # Sem espalhamento vertical
    alpha = 0.15,       # Muito transparente
    size = 1.5          # Pequeno
  ) +
  # Cores dos boxplots
  scale_fill_manual(
    values = c("0" = "#c40e3e", "1" = "#2DA5E8"),
    labels = c("0" = "Global South", "1" = "Global North")
  ) +
  # Cores dos pontos (mesmas cores, mas muito transparentes)
  scale_color_manual(
    values = c("0" = "#c40e3e", "1" = "#2DA5E8"),
    guide = "none"  # Sem legenda para cores
  ) +
  # Rótulos do eixo X
  scale_x_discrete(
    labels = c("0" = "Global South", "1" = "Global North")
  ) +
  # Rótulos dos eixos e título
  labs(
    x = "",
    y = "Geographic distribution \ncompleteness (%)",
    title = "",
    fill = ""  # Título vazio para legenda de fill
  ) +
  # Tema personalizado
  theme(
    panel.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    ),
    plot.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    ),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(size = 15, margin = margin(t = 5)),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    axis.ticks.y = element_line(color = "black", size = 0.5),   # aqui
    axis.ticks.x = element_line(color = "black", size = 0.5),   # opcional
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  ) +
  annotate(
    "text",
    x = 2.5,                     # lado direito (Global North)
    y = 0.42,
    label = p_label,
    hjust = 1.05,              # encosta no canto
    vjust = 1.2,               # encosta no topo
    size = 4.5,
    fontface = "italic"
  )

# Plotando as curvas de acumulação + boxplot
load(
  file.path(
    local_directory, "02_data_analysis", "tetrapods_model.RData")
)

par(cex.main = 2,   # título principal
    cex.lab  = 1.5, # rótulos dos eixos
    cex.axis = 1.3) # números dos eixos

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

anos <- c(1756, seq(1775, 2025, by = 25))
axis(1, at = anos, labels = anos, las = 2)  # labels verticais
abline(v = anos, 
       col = "gray70",  
       lty = 3,         
       lwd = 0.5)
# sobrepor boxplot
vp <- grid::viewport(x = 0.28, y = 0.59, width = 0.43, height = 0.60)
print(boxplot, vp = vp)
# export 12x8 pdf

# Countries and democracy index ----
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

df_01perc <- left_join(
  species_admin_01perc,
  data_sociopolitic,
  by = "name_en") %>%
  select(completeness_all, expected_rich, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>% 
  remove_missing() %>%
  ungroup() %>% 
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(expected_rich))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  )

df_05perc <- left_join(
  species_admin_05perc,
  data_sociopolitic,
  by = "name_en") %>%
  select(completeness_all, expected_rich, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>% 
  remove_missing() %>%
  ungroup() %>% 
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(expected_rich))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  )

df_1perc <- left_join(
  species_admin_1perc,
  data_sociopolitic,
  by = "name_en") %>%
  select(completeness_all, expected_rich, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>% 
  remove_missing() %>%
  ungroup() %>% 
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(expected_rich))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  ) 

df_occ_polygons <- left_join(
  species_admin_occ_polygons,
  data_sociopolitic,
  by = "name_en") %>%
  select(completeness_all, expected_rich, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions) %>% 
  remove_missing() %>%
  ungroup() %>% 
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(expected_rich))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)))
  ) 

modelo_beta01 <- betareg(completeness_all ~ 
                           polygons_richness_scaled + 
                           dens_pop_scaled + 
                           mean_grp_scaled +
                           democracy_index_scaled +
                           colonial_origin +
                           n_institutions_scale,
                         data = df_01perc)
summary(modelo_beta01)

modelo_beta05 <- betareg(completeness_all ~ 
                           polygons_richness_scaled + 
                           dens_pop_scaled + 
                           mean_grp_scaled +
                           democracy_index_scaled +
                           colonial_origin +
                           n_institutions_scale,
                         data = df_05perc)
summary(modelo_beta05)

modelo_beta1 <- betareg(completeness_all ~ 
                          polygons_richness_scaled + 
                          dens_pop_scaled + 
                          mean_grp_scaled +
                          democracy_index_scaled +
                          colonial_origin +
                          n_institutions_scale,
                        data = df_1perc)
summary(modelo_beta1)

modelo_betaocc <- betareg(completeness_all ~ 
                            polygons_richness_scaled + 
                            dens_pop_scaled + 
                            mean_grp_scaled +
                            democracy_index_scaled +
                            colonial_origin +
                            n_institutions_scale,
                          data = df_occ_polygons)
summary(modelo_betaocc)

coef_beta01 <- tidy(modelo_beta01) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map))

coef_beta05 <- tidy(modelo_beta05) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map))

coef_beta1 <- tidy(modelo_beta1) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map))

coef_occ <- tidy(modelo_betaocc) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map))

coef_beta01$model  <- "Threshold = 0.01%"
coef_beta05$model  <- "Threshold = 0.05%"
coef_beta1$model   <- "Threshold = 1%"
coef_occ$model     <- "Occurrence&Polygon"

coef_list <- bind_rows(
  coef_beta01,
  coef_beta05,
  coef_beta1,
  coef_occ
) %>%
  mutate(significant = factor(p.value < 0.05,
                              levels = c(TRUE, FALSE)))

coef_list$model <- factor(
  coef_list$model,
  levels = c("Occurrence&Polygon",
             "Threshold = 0.01%",
             "Threshold = 0.05%",
             "Threshold = 1%")
)

plot <- ggplot(coef_list, aes(x = estimate, 
                      y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "gray50", 
             linewidth = 0.5) +
  geom_point(aes(fill = significant),
             size = 3,
             shape = 21,
             color = "black",
             stroke = 1) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error),
    height = 0.1,
    color = "black",
    linewidth = 1
  ) +
  facet_wrap(~ model, ncol = 2) +
  labs(
    x = "Model coefficients (CI 95%)",
    y = ""
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  scale_fill_manual(
    values = c("TRUE" = "black",
               "FALSE" = "white"),
    guide = "none"
  ); plot

ggsave(file.path(drive, "FigS2.tiff"), 
       plot = plot, 
       width = 8, 
       height = 12, 
       dpi = 300, 
       bg = "white")

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

