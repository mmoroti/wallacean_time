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
  left_join(mobilization_effort, by = "name_en") %>%
  select(percent_wallace_prop, polygons_richness, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth, name_en, area_km2,
         colonial_origin, n_institutions, mobilization_effort) %>%
  remove_missing() %>%
  mutate(
    n_institutions_scale = scale(n_institutions),
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    dens_pop_scaled = scale(log(as.numeric(pop/area_km2))),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean))),
    mobilization_effort_scaled = scale(log(as.numeric(mobilization_effort)))
  ) 

glimpse(df_all)
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

# Plot dos coeficientes 
# Criar uma tabela de mapeamento para nomes mais legíveis
new_names <- tibble::tribble(
  ~term, ~term_novo,
  "polygons_richness_scaled", "Species richness",
  "dens_pop_scaled", "Population density",
  "mean_grp_scaled", "GDP per capta",
  "democracy_index_scaled", "Democracy index",
  "colonial_origin", "Colonial origin",
  "n_institutions_scale", "N. of research \ninstitutions"
) #  "mobilization_effort_scaled", "Mobilization effort"
map <- setNames(new_names$term_novo, new_names$term)

coef_list <- tidy(modelo_beta) %>%
  filter(!term %in% c("(Intercept)", "(phi)", "Log(nu)"))  %>%
  mutate(term = dplyr::recode(term, !!!map))

ggplot(coef_list, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error, 
        xmax = estimate + 1.96 * std.error),
    height = 0.1, color = "black", linewidth = 1
  ) +
  labs(
    x = "Model coefficients (CI 95%)",
    y = "",
    title = "",
    subtitle = ""
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  # Alterações para aumentar e colocar em negrito as variáveis
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(face = "bold", size = 12),  # Aumentar título do eixo X
    axis.title.x = element_text(margin = margin(t = 10)),  # Espaçamento para o título
    axis.text.y = element_text(
      size = 11,          # Tamanho aumentado das variáveis
      face = "bold",      # Texto em negrito
      color = "black"     # Cor preta
    ),
    axis.text.x = element_text(
      size = 10,          # Tamanho dos números no eixo X
      color = "black"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # Adicionar borda preta fina ao redor do painel
    panel.border = element_rect(
      color = "black",     # Cor da borda
      fill = NA,           # Não preencher
      linewidth = 0.5      # Espessura da borda
    ),
    # Ajustar margens se necessário
    plot.margin = margin(10, 10, 10, 10)  # Topo, direita, baixo, esquerda
  )

# Grupos combinados ----
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

ggplot(df_all, aes(x = factor(GlobalNorth),
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
    # Fundo branco
    panel.background = element_rect(
      fill = "white",
      color = "black",        # Borda preta
      linewidth = 0.5         # Linha fina
    ),
    plot.background = element_rect(
      fill = "white",
      color = "black",        # Borda preta ao redor de todo o gráfico
      linewidth = 0.5
    ),
    # Remover grid de fundo (opcional - para deixar mais limpo)
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    # Eixos
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(
      #face = "bold",          # Negrito nos rótulos
      size = 15,
      margin = margin(t = 5)  # Espaço acima
    ),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(
      size = 15,
      margin = margin(r = 10)  # Espaço à direita
    ),
    # Legenda (se quiser mantê-la, mas discreta)
    legend.position = "none",  # Remove completamente a legenda
    # Margens
    plot.margin = margin(15, 15, 15, 15)  # Espaço interno da borda
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

