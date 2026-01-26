library(tidyverse)
library(sf)
library(stringi)
library(betareg)
library(cowplot)
library(ggstats)
library(ggrepel)
library(ggpubr)
library(car)

## Function ----
plot_bivariate_map <- function(geographic_shape_data,
                               richness_completude_amphibia,
                               title = NULL) {
  
  # Calcular a porcentagem de completude Wallace
  # Aqui criamos uma nova coluna "percent_wallace" = (riqueza observada / riqueza total) * 100
  richness_completude <- richness_completude_amphibia %>%
    mutate(percent_wallace = observed_richness * 100 / polygons_richness)
  
  # Juntar com o shapefile (geometrias)
  # Fazemos um left_join pelo nome da unidade administrativa
  # Em seguida transformamos para o sistema de coordenadas projetadas EPSG:6933
  dados_join <- geographic_shape_data %>%
    left_join(richness_completude, by = c("name_en" = "name_en")) %>%
    st_transform(6933)
  
  # Criar classes de riqueza e completude em 4 categorias usando quantis
  # Aqui usamos 'cut' para transformar valores contínuos em categorias ordenadas
  dados_bi <- dados_join %>%
    mutate(
      richness_cat = cut(observed_richness,
                         breaks = quantile(observed_richness, probs = seq(0, 1, length.out = 5), na.rm = TRUE),
                         include.lowest = TRUE,
                         labels = c("very low", "low", "medium", "high")),
      percent_cat = cut(percent_wallace,
                        breaks = quantile(percent_wallace, probs = seq(0, 1, length.out = 5), na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = c("very low", "low", "medium", "high"))
    ) %>%
    mutate(
      bi_class = paste0(richness_cat, "-", percent_cat)
    )
  
  # Criar paleta de cores para cada combinação 4x4
  bi_pal <- tibble(
    bi_class = c(
      "very low-very low", "low-very low", "medium-very low", "high-very low",
      "very low-low",      "low-low",      "medium-low",      "high-low",
      "very low-medium",   "low-medium",   "medium-medium",   "high-medium",
      "very low-high",     "low-high",     "medium-high",     "high-high"
    ),
    color = c(
      "#f0f0f0", "#cfe7e7", "#9fd0d0", "#5ababa",
      "#e6cfe6", "#c5cfe0", "#94b8c9", "#4fa3b3",
      "#d7a6d7", "#b69fcf", "#8688b8", "#4373a1",
      "#c179c1", "#9f73b9", "#6f5ca2", "#2d478b"
    )
  )
  
  # Adicionar a cor correspondente em cada polígono
  dados_bi <- dados_bi %>%
    left_join(bi_pal, by = "bi_class")
  
  # Criar dados para a legenda 4x4
  # Aqui também transformamos as categorias em números para posicionar no grid
  legend_data <- expand.grid(
    richness_cat = c("very low", "low", "medium", "high"),
    percent_cat  = c("very low", "low", "medium", "high")
  ) %>%
    mutate(
      bi_class = paste0(richness_cat, "-", percent_cat)
    ) %>%
    left_join(bi_pal, by = "bi_class") %>%
    mutate(
      richness_num = case_when(
        richness_cat == "very low" ~ 1,
        richness_cat == "low"      ~ 2,
        richness_cat == "medium"   ~ 3,
        richness_cat == "high"     ~ 4
      ),
      percent_num = case_when(
        percent_cat == "very low" ~ 1,
        percent_cat == "low"      ~ 2,
        percent_cat == "medium"   ~ 3,
        percent_cat == "high"     ~ 4
      )
    )
  
  # Criar o mapa principal
  map_plot <- ggplot(dados_bi) +
    geom_sf(aes(fill = color), color = "black", size = 0.2) +
    scale_fill_identity() +
    theme_minimal() +
    labs(title = title) +
    theme(legend.position = "none")
  
  # Criar a legenda bivariada
  legend_plot <- ggplot(legend_data, aes(x = percent_num, y = richness_num, fill = color)) +
    geom_tile() +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = 1:4,
      labels = c("Very Low", "Low", "Medium", "High"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = 1:4,
      labels = c("Very Low", "Low", "Medium", "High"),
      expand = c(0, 0)
    ) +
    labs(x = "Knowledge of \ngeographic distribution", y = "Richness") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1), # rotação do eixo X
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
    )
  
  # ️Combinar mapa e legenda no mesmo layout
  final_plot <- ggdraw() +
    draw_plot(map_plot) +
    draw_plot(legend_plot,
              x = 0.07, y = 0.1, 
              width = 0.20, height = 0.35)
  
  # Retornar o gráfico final
  return(final_plot)
}

# with numbers
## Load data ----
# Lista de países/unidades do Sul Global (exemplo baseado em América Latina, África, Sul/Sudeste da Ásia, Oceania não-ocidental)
load(file.path(
  local_directory,
  "02_data_analysis",
  "richness_completude.RData"))

load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData"))

# Plots
plot_bivariate_map(geographic_shape_data,
                   richness_completude_amphibia,
                   title = "Amphibia")

plot_bivariate_map(geographic_shape_data,
                   richness_completude_reptilia,
                   title = "Reptilia")

plot_bivariate_map(geographic_shape_data, 
                   richness_completude_aves,
                   title = "Aves")

plot_bivariate_map(geographic_shape_data,
                   richness_completude_mammalia,
                   title = "Mammalia")

# correlacoes entre os preditores e range size
# Carregar variaveis socieconomicas
load(
  file = file.path(
    local_directory,
    "01_data_cleaned",
    "data_socieconomic_agregatte.RData")
)

richness_all <- richness_completude_amphibia %>%
  left_join(richness_completude_amphibia, by = "name_en") %>%
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

View(richness_all)
#df_completeness_adj <- richness_completude_reptilia %>%
#  mutate(
#    name_en = stri_trans_general(tolower(name_en), "Latin-ASCII"))
#
df_all <- left_join(
  richness_all,
  data_sociopolitic,
  by = "name_en") 

df <- df_all %>%
  select(polygons_richness, percent_wallace_prop, democracy_index_mean,
         academic_freedom_mean, gdp_mean, pop, GlobalNorth) %>%
  mutate(
    polygons_richness_scaled = scale(log(as.numeric(polygons_richness))),
    democracy_index_scaled = scale(log(as.numeric(democracy_index_mean))),
    academic_freedom_scaled = scale(as.numeric(academic_freedom_mean)),
    mean_pop_scaled = scale(log(as.numeric(pop)+1)),
    mean_grp_scaled = scale(log(as.numeric(gdp_mean)+1))
  ) %>% 
  remove_missing()
table(df$GlobalNorth)

# ajustar glm binomial
model_glm <- glm(percent_wallace_prop ~ 
                   polygons_richness_scaled + 
                   mean_pop_scaled + 
                   mean_grp_scaled +
                   democracy_index_scaled, data = df)
summary(model_glm)
performance::r2(model_glm) # R2: 0.419

ggstats::ggcoef_model(model_glm)

table(df$GlobalNorth)
vif(model_glm)

# Ajustar beta regression (percent_wallace já em proporção 0-1)
modelo_beta <- betareg(percent_wallace_prop ~ 
                       polygons_richness_scaled + 
                       mean_pop_scaled + 
                       mean_grp_scaled +
                       democracy_index_scaled,
                       data = df)
summary(modelo_beta)
car::vif(modelo_beta)

par(mfrow = c(1, 2))
plot(modelo_beta)
dev.off()

## sem perder muitos paises
#df <- df_all %>%
#  mutate(
#    latitude_scaled = scale(abs(latitude)),
#    polygons_richness_scaled = scale(log(polygons_richness/area_km2)),
#    percent_wallace = observed_richness * 100 / polygons_richness,
#    percent_wallace_prop = percent_wallace / 100
#  ) %>%
#  select(name_en, latitude_scaled, polygons_richness_scaled,
#         percent_wallace, percent_wallace_prop) %>%
#  remove_missing()

df %>%
  ggplot(aes(x = latitude_scaled, y = percent_wallace_prop, color = global_region)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              alpha = 0.5, color = "black") +  # linha de tendência preta
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = name_en), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Norte Global" = "#2DA5E8", "Sul Global" = "#c40e3e")) +
  theme_bw() +
  labs(
    x = "Latitude absoluta",
    y = "Proporção do Conhecimento Wallaceano (%)",
    color = ""
    )

richness_all %>%
  ggplot(aes(x = GlobalNorth, y = percent_wallace_prop, fill = GlobalNorth)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, aes(group = GlobalNorth)) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  scale_fill_manual(values = c("1" = "#2DA5E8", "0" = "#c40e3e")) +
  theme_bw() +
  labs(
    x = "Região geopolítica",
    y = "Proporção conhecimento Wallaceano (%)"
    )  

ggplot(df, aes(x = factor(GlobalNorth),
               y = democracy_index_scaled)) +
  geom_boxplot(aes(fill = as.factor(GlobalNorth))) +
  scale_fill_manual(values = c("1" = "#2DA5E8", "0" = "#c40e3e")) +
  labs(x = "Global North",
       y = "Democracy",
       title = "") +
  theme_minimal()

# Teste de Wilcoxon
wilcox.test(democracy_index_scaled ~ GlobalNorth, data = df)

hist(df$latitude_scaled)
hist(df$polygons_richness_scaled)
hist(df$mean_pop_scaled)
hist(df$mean_grp_scaled)

# Ajustar beta regression (percent_wallace já em proporção 0-1)
modelo_beta <- betareg(percent_wallace_prop ~ latitude_scaled + polygons_richness_scaled + 
                         mean_pop_scaled + mean_grp_scaled, data = df)
summary(modelo_beta)
car::vif(modelo_beta)

par(mfrow = c(1, 2))
plot(modelo_beta)
dev.off()

# ajustar glm binomial

model_glm <- glm(percent_wallace_prop ~ latitude_scaled + polygons_richness_scaled + 
                   mean_pop_scaled + mean_grp_scaled, data = df)
summary(model_glm)
performance::r2(model_glm) # R2: 0.419

ggstats::ggcoef_model(model_glm)
par(mfrow = c(2, 2))
plot(model_glm)



# simples
model_glm <- glm(percent_wallace_prop ~ 
                   latitude_scaled + polygons_richness_scaled, data = df) 
summary(model_glm)
