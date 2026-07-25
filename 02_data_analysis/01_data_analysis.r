# Packages and functions ----
library(tidyverse)
library(mets)
library(timereg)
library(sf)
library(cowplot)

rm(list=ls()); gc() # clean local enviroment

# Set directory 
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# Functions
extrair_exp_coefficients_simples <- function(output_text) {
  # Encontrar a linha com "exp(coeffients):"
  start_idx <- which(output_text == "exp(coeffients):")[1]
  
  if (is.na(start_idx)) {
    stop("Cabeçalho 'exp(coeffients):' não encontrado.")
  }
  
  # As linhas de dados começam 2 linhas após o cabeçalho
  data_lines <- output_text[(start_idx + 2):length(output_text)]
  
  # Encontrar onde termina a seção (linha vazia)
  empty_lines <- which(data_lines == "" | grepl("^\\s*$", data_lines))
  
  if (length(empty_lines) > 0) {
    data_lines <- data_lines[1:(min(empty_lines[1], length(data_lines)) - 1)]
  }
  
  # Processar cada linha
  results <- list()
  
  for (line in data_lines) {
    # Dividir por múltiplos espaços
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    
    if (length(parts) >= 4) {
      # Nome da variável (pode ser concatenado se houver múltiplas palavras)
      var_name <- paste(parts[1:(length(parts)-3)], collapse = " ")
      
      # Valores numéricos (últimos 3 elementos)
      estimate <- as.numeric(parts[length(parts)-2])
      lower_ci <- as.numeric(parts[length(parts)-1])
      upper_ci <- as.numeric(parts[length(parts)])
      
      results[[length(results) + 1]] <- data.frame(
        variable = var_name,
        hr = estimate,
        lower = lower_ci,
        upper = upper_ci,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(results) == 0) {
    return(NULL)
  }
  
  return(do.call(rbind, results))
}

plot_gof_custom <- function(x,
                            col = 3,
                            type = NULL,
                            titles = TRUE,
                            bg.col = NULL,
                            xlab = "",
                            ylab = "",
                            cex.lab = 1,
                            cex.axis = 1,
                            show.xlab = TRUE,
                            show.ylab = TRUE,
                            show.xaxis = TRUE,
                            show.yaxis = TRUE,
                            show.xticks = TRUE,
                            show.yticks = TRUE,
                            tick_seq = seq(1750, 2025, by = 25),
                            label_seq = seq(1750, 2025, by = 50),
                            variable = NULL,
                            show.varname = FALSE,
                            var.labels = NULL,
                            grid.lines = TRUE,
                            grid.col = "gray",
                            grid.lty = 2,
                            grid.lwd = 0.5,
                            tick.length = -0.02,  # novo parâmetro: comprimento dos ticks (negativo = para fora)
                            ...) {
  
  if (is.null(type)) {
    if (x$type == "prop") type <- "time"
    if (x$type == "modelmatrix") type <- "modelmatrix"
    if (x$type == "Zmodelmatrix") type <- "z"
  }
  
  xlab_final <- if (show.xlab) xlab else ""
  ylab_final <- if (show.ylab) ylab else ""
  
  get_varname <- function(i){
    
    nm <- rownames(x$res)[i]
    
    if(!is.null(var.labels)){
      if(nm %in% names(var.labels)){
        nm <- var.labels[[nm]]
      }
    }
    
    nm
  }
  
  if (type == "time" || type == "modelmatrix") {
    
    p <- ncol(x$score)
    
    vars <- if(is.null(variable)) seq_len(p) else variable
    
    for(i in vars){
      
      simU <- x$simUt[, (0:49) * p + i]
      
      rsU <- max(abs(simU))
      rsU <- max(rsU, abs(x$score[, i]))
      
      plot(x$jumptimes,
           x$score[, i],
           type = "n",
           ylim = c(-rsU, rsU),
           xlim = range(tick_seq),
           xlab = xlab_final,
           ylab = ylab_final,
           cex.lab = cex.lab,
           cex.axis = cex.axis,
           xaxt = "n",
           yaxt = "n",
           ...)
      
      if(grid.lines) {
        abline(v = seq(1750, 2025, by = 50), 
               col = grid.col, 
               lty = grid.lty, 
               lwd = grid.lwd)
      }
      
      if(show.xaxis || show.xticks){
        
        axis(
          1,
          at = tick_seq,
          labels = if(show.xaxis){
            ifelse(tick_seq %in% label_seq, tick_seq, "")
          } else {
            rep("", length(tick_seq))
          },
          tick = show.xticks,
          cex.axis = cex.axis,
          tcl = tick.length  # adicionado: controla o comprimento dos ticks
        )
        
      }
      
      if(show.yaxis || show.yticks){
        
        axis(
          2,
          labels = show.yaxis,
          tick = show.yticks,
          cex.axis = cex.axis,
          tcl = tick.length  # adicionado: controla o comprimento dos ticks
        )
        
      }
      
      box()
      
      if(titles)
        title(main = get_varname(i))
      
      if(show.varname)
        mtext(
          get_varname(i),
          side = 2,
          line = 2.8,
          las = 0,
          font = 2,
          cex = 0.9
        )
      
      matlines(
        x$jumptimes,
        simU,
        type = "s",
        lwd = 0.3,
        col = if(!is.null(bg.col)) bg.col else col
      )
      
      lines(
        x$jumptimes,
        x$score[, i],
        type = "s",
        lwd = 1.5
      )
      
    }
    
  } else {
    
    vars <- if(is.null(variable)) seq_along(x$Zres) else variable
    
    for(i in vars){
      
      xr <- x$Zres[[i]]
      
      obsz <- tail(xr$score, 1)
      
      times <- xr$xaxs
      
      rsU <- max(
        max(abs(obsz)),
        max(abs(xr$simUtlast[1:50, ]))
      )
      
      plot(times,
           obsz,
           type = "n",
           ylim = c(-rsU, rsU),
           xlim = range(tick_seq),
           xlab = xlab_final,
           ylab = ylab_final,
           cex.lab = cex.lab,
           xaxt = "n",
           yaxt = "n",
           ...)
      
      if(grid.lines) {
        abline(v = seq(1750, 2025, by = 50), 
               col = grid.col, 
               lty = grid.lty, 
               lwd = grid.lwd)
      }
      
      if(show.xaxis || show.xticks){
        
        axis(
          1,
          at = tick_seq,
          labels = if(show.xaxis){
            ifelse(tick_seq %in% label_seq, tick_seq, "")
          } else {
            rep("", length(tick_seq))
          },
          tick = show.xticks,
          cex.axis = cex.axis,
          tcl = tick.length  # adicionado: controla o comprimento dos ticks
        )
        
      }
      
      if(show.yaxis || show.yticks){
        
        axis(
          2,
          labels = show.yaxis,
          tick = show.yticks,
          cex.axis = cex.axis,
          tcl = tick.length  # adicionado: controla o comprimento dos ticks
        )
        
      }
      
      box()
      
      if(titles)
        title(main = get_varname(i))
      
      if(show.varname)
        mtext(
          get_varname(i),
          side = 2,
          line = 2.8,
          las = 0,
          font = 2,
          cex = 0.9
        )
      
      matlines(
        times,
        t(xr$simUtlast[1:50, ]),
        type = "l",
        lwd = 0.3,
        col = if(!is.null(bg.col)) bg.col else col
      )
      lines(
        times,
        obsz,
        lwd = 2,
        col = 1
      )
    }
  }
}

plot_surv_panel <- function(fit_object,
                            comps = 3:10,
                            show_main = FALSE,
                            show_xlab = FALSE,
                            show_ylab = FALSE,
                            show_years = TRUE,
                            year_seq = seq(1750, 2025, by = 50),
                            tick_seq = seq(1750, 2025, by = 25),
                            label_seq = seq(1750, 2025, by = 50),
                            coef_col = "blue",
                            coef_lwd = 2,
                            shade_alpha = 0.15,
                            var_index = NULL,
                            tick.length = -0.02,
                            xlim_range = c(1750, 2025),
                            cex.axis = 1) {
  
  if(!is.null(var_index)) {
    comps <- comps[var_index]
  }
  
  for(i in seq_along(comps)) {
    comp <- comps[i]
    
    time <- fit_object$cum[, 1]
    est <- fit_object$cum[, comp]
    sd <- sqrt(fit_object$robvar.cum[, comp])
    upper <- est + 2*sd
    lower <- est - 2*sd
    
    xlab <- ""
    ylab <- ""
    main <- ""
    
    ylim <- range(c(lower, upper, 0))
    
    plot(time, est, type = "n", xlab = xlab, ylab = ylab,
         main = main, ylim = ylim, xlim = xlim_range,
         xaxt = "n", yaxt = "n")
    
    # 1. PRIMEIRO: Grid de fundo (a cada 50 anos)
    abline(v = label_seq, col = "gray90", lty = 1, lwd = 0.5)
    
    # 2. SEGUNDO: Linha horizontal em zero
    abline(h = 0, col = "gray50", lty = 2)
    
    # 3. TERCEIRO: Polígono (shade) - por cima do grid
    polygon(c(time, rev(time)),
            c(upper, rev(lower)),
            col = adjustcolor(coef_col, alpha.f = shade_alpha),
            border = NA)
    
    # 4. QUARTO: Linha central - por cima do polígono
    lines(time, est, col = coef_col, lwd = coef_lwd)
    
    # 5. QUINTO: Eixos
    axis(1, 
         at = tick_seq,
         labels = ifelse(tick_seq %in% label_seq, tick_seq, ""),
         las = 2, 
         cex.axis = cex.axis,
         tcl = tick.length)
    
    axis(2, 
         cex.axis = cex.axis,
         tcl = tick.length)
    
    box()
  }
}

plot_bivariate_map <- function(geo_plot,
                               var_x,
                               var_y,
                               title = NULL,
                               x_label = var_x,
                               y_label = var_y) {
  
  # 1) transforma no CRS desejado
  dados_join <- geo_plot %>%
    st_transform(6933)
  
  # garantir que as colunas existem
  if (!var_x %in% names(dados_join)) {
    stop(paste("A coluna", var_x, "não existe no dataframe."))
  }
  if (!var_y %in% names(dados_join)) {
    stop(paste("A coluna", var_y, "não existe no dataframe."))
  }
  
  x_vals <- dados_join[[var_x]]
  y_vals <- dados_join[[var_y]]
  
  # 2) quantis da variável X
  x_valid <- x_vals[!is.na(x_vals) & x_vals > 0]
  
  if (length(x_valid) < 2) {
    stop(paste("A coluna", var_x, "não tem valores suficientes maiores que zero."))
  }
  
  x_q_raw <- quantile(x_valid, probs = seq(0, 1, length.out = 5), na.rm = TRUE)
  x_breaks <- as.integer(floor(x_q_raw))
  
  for(i in 2:length(x_breaks)){
    if(x_breaks[i] <= x_breaks[i-1]){
      x_breaks[i] <- x_breaks[i-1] + 1
    }
  }
  
  x_breaks[1] <- min(x_valid)
  x_breaks[length(x_breaks)] <- max(x_valid)
  
  # 3) quantis da variável Y
  y_valid <- y_vals[!is.na(y_vals)]
  
  if (length(y_valid) < 2) {
    stop(paste("A coluna", var_y, "não tem valores suficientes não-NA."))
  }
  
  y_breaks <- quantile(y_valid, probs = seq(0, 1, length.out = 5), na.rm = TRUE)
  
  # 4) funções de label - SEM CASAS DECIMAIS
  make_labels_int <- function(breaks, suffix = "") {
    paste0(
      round(breaks[-length(breaks)], 0),  # arredonda para inteiro
      "–",
      round(breaks[-1], 0),               # arredonda para inteiro
      suffix
    )
  }
  
  make_labels_pct <- function(breaks, suffix = "") {
    paste0(
      round(breaks[-length(breaks)], 0),  # arredonda para inteiro
      "–",
      round(breaks[-1], 0),               # arredonda para inteiro
      suffix
    )
  }
  
  x_labs <- make_labels_int(x_breaks, "")
  y_labs <- make_labels_pct(y_breaks, "")
  
  # 5) categorias
  dados_bi <- dados_join %>%
    mutate(
      x_cat = case_when(
        is.na(.data[[var_x]]) | .data[[var_x]] == 0 ~ NA_character_,
        TRUE ~ as.character(cut(.data[[var_x]],
                                breaks = x_breaks,
                                include.lowest = TRUE,
                                labels = x_labs))
      ),
      y_cat = case_when(
        is.na(.data[[var_y]]) ~ NA_character_,
        TRUE ~ as.character(cut(.data[[var_y]],
                                breaks = y_breaks,
                                include.lowest = TRUE,
                                labels = y_labs))
      ),
      bi_class = ifelse(is.na(x_cat) | is.na(y_cat),
                        NA_character_,
                        paste0(x_cat, " | ", y_cat))
    )
  
  # 6) paleta 4x4
  x_labs_factor <- factor(x_labs, levels = x_labs, ordered = TRUE)
  y_labs_factor <- factor(y_labs, levels = y_labs, ordered = TRUE)
  
  bi_pal <- expand.grid(
    x_cat = x_labs_factor,
    y_cat = y_labs_factor,
    stringsAsFactors = FALSE
  ) %>%
    arrange(x_cat, y_cat) %>%
    mutate(
      bi_class = paste0(x_cat, " | ", y_cat),
      fill = c(
        "#4e3d43", "#535e68", "#577e8b", "#5b9cad",
        "#724448", "#786970", "#7e8d96", "#84afba",
        "#944b4d", "#9c7478", "#a59ba0", "#acc1c6",
        "#b65252", "#c07f7f", "#caaaaa", "#d3d3d3"
      )
    )
  
  # 7) juntar cor
  dados_bi <- dados_bi %>%
    left_join(bi_pal %>% select(bi_class, fill), by = "bi_class") %>%
    mutate(fill = ifelse(is.na(fill), "#ffffff", fill))
  
  # 8) mapa
  bbox <- st_bbox(dados_bi)
  bbox_poly <- st_as_sfc(bbox)
  
  map_plot <- ggplot(dados_bi) +
    geom_sf(aes(fill = fill), color = "black", size = 0.15) +
    geom_sf(data = bbox_poly, fill = NA, color = "black", linewidth = 0.3) +
    scale_fill_identity() +
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
      axis.text = element_text(size = 6, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      axis.ticks.length = unit(2, "mm"),
      panel.grid = element_line(color = "grey80", linewidth = 0.2),
      plot.margin = margin(5, 15, 5, 5)
    )
  
  # 9) legenda
  legend_data <- bi_pal
  
  legend_plot <- ggplot(legend_data,
                        aes(x = x_cat,
                            y = y_cat,
                            fill = fill)) +
    geom_tile() +
    scale_fill_identity() +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 7, color = "black"),
      axis.text.x = element_text(size = 5, angle = 45, hjust = 0.8, color = "black"),
      axis.text.y = element_text(size = 5, color = "black"),
      panel.grid = element_blank()
    )
  
  # 10) combinar
  final_plot <- ggdraw() +
    draw_plot(map_plot) +
    draw_plot(legend_plot,
              x = 0.065, y = 0.28,
              width = 0.20, height = 0.30)
  
  return(final_plot)
}

normalize <- function(x){
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Load data ----
load(file = file.path(
  local_directory,
  "02_data_analysis",
  "data_wallacean_time.RData")
)

# shapefile to crop adm unit
load(file.path(
  local_directory,
  "00_raw_data",
  "geographic_shape_data.RData"))

# Richness completeness
load(file.path(
  local_directory,
  "02_data_analysis",
  "richness_completude.RData"))

# Preparing data ----
# Thresholds of completeness 
df_wallacean_01 <- df_wallacean_01 %>%
  filter(ua_total01 > 2) %>%
  mutate(status = event + WallaceCompletude01) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  )

df_wallacean_05 <- df_wallacean_05 %>%
  filter(ua_total05 > 2) %>%
  mutate(status = event + WallaceCompletude05) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  )

df_wallacean_1 <- df_wallacean_1 %>%
  filter(ua_total1 > 2) %>%
  mutate(status = event + WallaceCompletude1) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  )

# Each group
df_tetrapods_100 <- df_wallacean_100 %>%
  filter(ua_total > 2) %>%
  mutate(status = event + WallaceCompletude) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  ) %>%
  mutate(status = event + WallaceCompletude) %>%
  mutate(BodySize = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensityTime+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation)),
         gdp_rank_z = scale(gdp_rank),
         gdp_percentile_z = scale(gdp_percentile),
         gdp = scale(log_gdp_final))

df_amphibia_100 <- df_wallacean_100 %>%
  filter(ua_total > 2) %>%
  mutate(status = event + WallaceCompletude) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  ) %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Amphibia") %>%
  mutate(BodySize = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensityTime+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation)),
         gdp_rank_z = scale(gdp_rank),
         gdp_percentile_z = scale(gdp_percentile),
         gdp = scale(log_gdp_final))

df_reptilia_100 <- df_wallacean_100 %>%
  filter(ua_total > 2) %>%
  mutate(status = event + WallaceCompletude) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  ) %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Reptilia") %>%
  mutate(BodySize = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensityTime+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation)),
         gdp_rank_z = scale(gdp_rank),
         gdp_percentile_z = scale(gdp_percentile),
         gdp = scale(log_gdp_final))

df_aves_100 <- df_wallacean_100 %>%
  filter(ua_total > 2) %>%
  mutate(status = event + WallaceCompletude) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  ) %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Aves") %>%
  mutate(BodySize = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensityTime+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation)),
         gdp_rank_z = scale(gdp_rank),
         gdp_percentile_z = scale(gdp_percentile),
         gdp = scale(log_gdp_final))

df_mammalia_100 <- df_wallacean_100 %>%
  filter(ua_total > 2) %>%
  mutate(status = event + WallaceCompletude) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  ) %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Mammalia") %>%
  mutate(BodySize = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensityTime+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation)),
         gdp_rank_z = scale(gdp_rank),
         gdp_percentile_z = scale(gdp_percentile),
         gdp = scale(log_gdp_final))

# Recurrent event model ----
# All species
tetrapods_100 <- recurrentMarginal(
  Event(t.start.year, t.stop.year, status==1)~
    strata(GlobalNorth)+
    cluster(speciesKey),
    cause=1,
    death.code=2,
  data=df_wallacean_100 %>%
    mutate(status = event + WallaceCompletude) %>%
    left_join(geographic_shape_data %>%
        select(GlobalNorth, Georegion, name_en) %>%
        sf::st_drop_geometry(), by = "name_en"
    ))

save(
  tetrapods_100,
  file = file.path(
    local_directory, "02_data_analysis", "tetrapods_model.RData")
) # plot in script "02_exploratory_analysis.R"

# For widely distributed species (at least 3 administrative units)
# timereg aalen estimator
#fit.tetrapods.surv <- aalen(
#  formula = Surv(t.start.year, t.stop.year, status == 1) ~
#    BodySize +
#    Verticality +
#    Nocturnality +
#    RangeSize +
#    Elevation +
#    HumanDensity +
#    Latitude +
#    gdp_percentile_z +
#    cluster(speciesKey)+
#    strata(Class),
#  data= df_tetrapods_100,
#  start.time = 1760,
#  max.time     = 2025,
#  #residuals = 1,
#  robust       = 1,       # variância robusta
#  n.sim        = 100,     # número de simulações para CI
#  resample.iid = 1       # para funções de residuais/CI i.i.d.
#)
#
#fit.tetrapods.completeness <- aalen(
#  formula = Surv(t.start.year, t.stop.year, status == 2) ~
#    BodySize +
#    Verticality +
#    Nocturnality +
#    RangeSize +
#    Elevation +
#    HumanDensity +
#    Latitude +
#    gdp_percentile_z +
#    cluster(speciesKey)+
#    strata(Class),
#  data= df_tetrapods_100,
#  start.time = 1760,
#  max.time     = 2025,
#  #residuals = 1,
#  robust       = 1,       # variância robusta
#  n.sim        = 100,     # número de simulações para CI
#  resample.iid = 1       # para funções de residuais/CI i.i.d.
#)
#
#fit.tetrapods <- recurrent.marginal.mean(fit.tetrapods.surv,
#                                         fit.tetrapods.completeness)
#

# mets package - cox estimator
tetrapods_event_100 <- phreg(Surv(t.start.year, t.stop.year, status==1)~
                               BodySize +
                               Verticality +
                               Nocturnality +
                               RangeSize +
                               Elevation +
                               HumanDensity +
                               Latitude +
                               gdp_percentile_z +
                               cluster(speciesKey)+
                               strata(Class),
                             data=df_tetrapods_100)

tetrapods_terminal_100 <- phreg(Surv(t.start.year,
                                     t.stop.year,
                                     status==2)~
                                  BodySize +
                                  Verticality +
                                  Nocturnality +
                                  RangeSize +
                                  Elevation +
                                  HumanDensity +
                                  Latitude +
                                  gdp_percentile_z +
                                  cluster(speciesKey)+
                                  strata(Class),
                                data=df_tetrapods_100)

tetrapods_100 <- recurrentMarginal(Event(t.start.year,
                                         t.stop.year,
                                         status==1)~
                                     BodySize +
                                     Verticality +
                                     Nocturnality +
                                     RangeSize +
                                     Elevation +
                                     HumanDensity +
                                     Latitude +
                                     gdp_percentile_z +
                                     strata(Class)+
                                     cluster(speciesKey),
                                   data= df_tetrapods_100,
                                   cause=1,
                                   death.code=2)

## Amphibia ----
# Cox model
amphibia_event_100 <- phreg(Surv(t.start.year, t.stop.year, status==1)~
                              BodySize +
                              Verticality +
                              Nocturnality +
                              RangeSize +
                              Elevation +
                              HumanDensity +
                              Latitude +
                              gdp_percentile_z +
                              cluster(speciesKey),
                            data=df_amphibia_100)

summary(amphibia_event_100)

gof.amphibia <- gof(amphibia_event_100)

amphibia_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~
                                 BodySize +
                                 Verticality +
                                 Nocturnality +
                                 RangeSize +
                                 Elevation +
                                 HumanDensity +
                                 Latitude +
                                 gdp_percentile_z +
                                 cluster(speciesKey),
                               data=df_amphibia_100)
summary(amphibia_terminal_100)

gof.terminal.amphibia <- gof(amphibia_terminal_100)

# Aalen model
fit.amphibia.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_amphibia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)
plot(fit.amphibia.surv)

# Extrair a lista de resíduos
#residuos_do_modelo <- fit.amphibia.surv$residuals
#str(residuos_do_modelo)
# Passar o objeto diretamente para a função de validação
#resultado_validacao <- cum.residuals(fit.amphibia.surv,
#                                     data = df_amphibia_100,
#                                     cum.resid = 1,
#                                     n.sim = 100)
# 3. Visualizar os resultados
#plot(resultado_validacao, score = 2)
# A lista geralmente contém:
# - dM: os incrementos do martingale estimados (resíduos)
# - time: o vetor de tempo correspondente a cada incremento

fit.amphibia.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_amphibia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)
plot(fit.amphibia.completeness)

fit.amphibia <- recurrent.marginal.mean(
  fit.amphibia.surv, 
  fit.amphibia.completeness
)

## Reptilia ----
# Cox model
reptilia_event_100 <- phreg(Surv(t.start.year,
                                 t.stop.year,
                                 status==1)~
                              BodySize +
                              Verticality +
                              Nocturnality +
                              RangeSize +
                              Elevation +
                              HumanDensity +
                              Latitude +
                              gdp_percentile_z +
                              cluster(speciesKey),
                            data=df_reptilia_01)
summary(reptilia_event_100)
gof.reptilia <- gof(reptilia_event_100)

reptilia_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~
                                 BodySize +
                                 Verticality +
                                 Nocturnality +
                                 RangeSize +
                                 Elevation +
                                 HumanDensity +
                                 Latitude +
                                 gdp_percentile_z +
                                 cluster(speciesKey),
                               data=df_reptilia_01)
summary(reptilia_terminal_100)
gof.terminal.reptilia <- gof(reptilia_terminal_100)

# Aalen model
fit.reptilia.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_reptilia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.reptilia.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_reptilia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.reptilia <- recurrent.marginal.mean(
  fit.reptilia.surv,
  fit.reptilia.completeness
)

## Aves ----
aves_event_100 <- phreg(Surv(t.start.year,
                             t.stop.year,
                             status==1)~
                          BodySize +
                          Verticality +
                          Nocturnality +
                          RangeSize +
                          Elevation +
                          HumanDensity +
                          Latitude +
                          gdp_percentile_z +
                          cluster(speciesKey),
                        data=df_aves_100)
summary(aves_event_100)
gof.aves <- gof(aves_event_100)

aves_terminal_100 <- phreg(Surv(t.start.year,
                                t.stop.year,
                                status==2)~
                             BodySize +
                             Verticality +
                             Nocturnality +
                             RangeSize +
                             Elevation +
                             HumanDensity +
                             Latitude +
                             gdp_percentile_z +
                             cluster(speciesKey),
                           data=df_aves_100)
summary(aves_terminal_100)
gof.terminal.aves <- gof(aves_terminal_100)

# Aalen model
fit.aves.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_aves_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.aves.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_aves_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.aves <- recurrent.marginal.mean(
  fit.aves.surv,
  fit.aves.completeness
)

## Mammals ----
mammalia_event_100 <- phreg(Surv(t.start.year,
                                 t.stop.year,
                                 status==1)~
                              BodySize +
                              Verticality +
                              Nocturnality +
                              RangeSize +
                              Elevation +
                              HumanDensity +
                              Latitude +
                              gdp_percentile_z +
                              cluster(speciesKey),
                            data=df_mammalia_100)
summary(mammalia_event_100)
gof.mammalia <- gof(mammalia_event_100)

mammalia_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~                              
                                 BodySize +
                                 Verticality +
                                 Nocturnality +
                                 RangeSize +
                                 Elevation +
                                 HumanDensity +
                                 Latitude +
                                 gdp_percentile_z +
                                 cluster(speciesKey),
                               data=df_mammalia_100)
summary(mammalia_terminal_100)
gof.terminal.mammalia <- gof(mammalia_terminal_100)

# Aalen model
fit.mammalia.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_mammalia_100 %>% filter(!Order %in% c("Chiroptera")),
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)
plot(fit.mammalia.surv)

fit.mammalia.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    BodySize +
    Verticality +
    Nocturnality +
    RangeSize +
    Elevation +
    HumanDensity +
    Latitude +
    gdp_percentile_z +
    cluster(speciesKey),
  data= df_mammalia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.mammalia <- recurrent.marginal.mean(
  fit.mammalia.surv,
  fit.mammalia.completeness
)

# Save models ----
# Cox models
save(
  amphibia_event_100,
  gof.amphibia,
  reptilia_event_100,
  gof.reptilia,
  aves_event_100,
  gof.aves,
  mammalia_event_100,
  gof.mammalia,
  tetrapods_event_100,
  tetrapods_100,
  file = file.path(
    local_directory,
    "02_data_analysis",
    "results",
    "cox_event_models.RData")
)
save(
  amphibia_terminal_100,
  gof.terminal.amphibia,
  reptilia_terminal_100,
  gof.terminal.reptilia,
  aves_terminal_100,
  gof.terminal.aves,
  mammalia_terminal_100,
  gof.terminal.mammalia,
  tetrapods_terminal_100,
  file = file.path(
    local_directory,
    "02_data_analysis",
    "results",
    "cox_terminal_models.RData")
)

# Aalen models
save(
  fit.amphibia.surv,
  fit.reptilia.surv,
  fit.aves.surv,
  fit.mammalia.surv,
  file = file.path(
    local_directory,
    "02_data_analysis",
    "results",
    "aalen_event_models.RData"))
save(
  fit.amphibia.completeness,
  fit.reptilia.completeness,
  fit.aves.completeness,
  fit.mammalia.completeness, 
  file = file.path(
    local_directory,
    "02_data_analysis",
    "results",
    "aalen_terminal_models.RData"))

#save(
#  fit.tetrapods,
#  fit.tetrapods.surv,
#  fit.tetrapods.completeness,
#  file = file.path(
#    local_directory,
#    "02_data_analysis",
#    "results",
#    "aalen_hazard_models.RData")
#)

# Results plots ----
# Exploratory data
length(unique(df_amphibia_100$speciesKey)) # 1955 spp
length(unique(df_reptilia_100$speciesKey)) # 3401 spp
length(unique(df_aves_100$speciesKey))     # 6073 spp
length(unique(df_mammalia_100$speciesKey)) # 2954 spp

table(df_amphibia_100$event)
table(df_reptilia_100$event)
table(df_aves_100$event)
table(df_mammalia_100$event)

load(file.path(
  local_directory,
  "02_data_analysis",
  "results",
  "cox_event_models.RData")
)

load(file.path(
  local_directory,
  "02_data_analysis",
  "results",
  "cox_terminal_models.RData")
)

#load(
#  file = file.path(local_directory,
#                   "02_data_analysis",
#                   "results",
#                   "aalen_hazard_models.RData")  
#)

cores <- c("Amphibia" = "#ed3a8b",
           "Reptilia" = "#00cc00",
           "Aves" = "#7a67ed",
           "Mammalia" = "#ffa400")

## Cumulative hazard plots ----
dev.off()
par(mfrow=c(1,3),
    cex.lab = 1.5, 
    cex.axis = 1.5)  # Aumentar tamanho dos labels dos eixos

ordem <- tetrapods_event_100$strata.level
cores <- cores[ordem]

plot(tetrapods_event_100, se = TRUE, col = cores,
     legend = FALSE,
     ylim = c(0, 35),
     xlim = c(1750, 2025),
     xlab = "",  # Omitir título do eixo x
     cex.main = 1.5)  # Aumentar tamanho do título
title(main = "Geodiscovery", line = 0.5, cex.main = 1.5)  # Ajustar posição do título
legend("topleft", legend = names(cores), col = unname(cores), 
       lwd = 2, bty = "n", title = "", cex = 1.5)  # Aumentar tamanho da legenda

anos <- seq(1750, 2025, by = 50)
abline(v = anos, col = "gray70", lty = 2, lwd = 1)

anos_25 <- seq(1750, 2025, by = 25)
anos_50 <- seq(1750, 2025, by = 50)
axis(1, at = anos_25, labels = FALSE, tck = -0.015)  # ticks menores
axis(1, at = anos_50, labels = anos_50, las = 1)     # labels só nos anos de 50

ordem <- tetrapods_terminal_100$strata.level
cores <- cores[ordem]
plot(tetrapods_terminal_100, se = TRUE, col = cores, 
     legend = FALSE,
     ylim = c(0, 0.25),
     xlim = c(1750, 2025),
     xlab = "",  # Omitir título do eixo x
     cex.main = 1.5)
title(main = "Completeness", line = 0.5, cex.main = 1.5)
anos <- seq(1750, 2025, by = 50)
abline(v = anos, col = "gray70", lty = 2, lwd = 1)
anos_25 <- seq(1750, 2025, by = 25)
anos_50 <- seq(1750, 2025, by = 50)
axis(1, at = anos_25, labels = FALSE, tck = -0.015)  # ticks menores
axis(1, at = anos_50, labels = anos_50, las = 1)     # labels só nos anos de 50

# Mean between event + terminal
plot(tetrapods_100, se = TRUE, col = cores, 
     legend = FALSE,
     ylim = c(0, 20),
     xlim = c(1750, 2025),
     xlab = "",
     cex.main = 1.5)
title(main = "Geographic distribution knowledge", line = 0.5, cex.main = 1.5)
anos <- seq(1750, 2025, by = 50)
abline(v = anos, col = "gray70", lty = 2, lwd = 1)
anos_25 <- seq(1750, 2025, by = 25)
anos_50 <- seq(1750, 2025, by = 50)
axis(1, at = anos_25, labels = FALSE, tck = -0.015)  # ticks menores
axis(1, at = anos_50, labels = anos_50, las = 1)     # labels só nos anos de 50

## Coefficient plot ----
# Preparar dados - exemplo para répteis eventos
amphibia_event <- extrair_exp_coefficients_simples(capture.output(summary(amphibia_event_100)))
amphibia_terminal <- extrair_exp_coefficients_simples(capture.output(summary(amphibia_terminal_100)))
amphibia_event$modelo <- "Geodiscovery"
amphibia_terminal$modelo <- "Completeness"
amphibia_combinado <- rbind(amphibia_event, amphibia_terminal)
amphibia_combinado$Class <- "Amphibia"

# Reptiles
reptilia_event <- extrair_exp_coefficients_simples(capture.output(summary(reptilia_event_100)))
reptilia_terminal <- extrair_exp_coefficients_simples(capture.output(summary(reptilia_terminal_100)))
reptilia_event$modelo <- "Geodiscovery"
reptilia_terminal$modelo <- "Completeness"
reptilia_combinado <- rbind(reptilia_event, reptilia_terminal)
reptilia_combinado$Class <- "Reptilia"

# Aves
aves_event <- extrair_exp_coefficients_simples(capture.output(summary(aves_event_100)))
aves_terminal <- extrair_exp_coefficients_simples(capture.output(summary(aves_terminal_100)))
aves_event$modelo <- "Geodiscovery"
aves_terminal$modelo <- "Completeness"
aves_combinado <- rbind(aves_event, aves_terminal)
aves_combinado$Class <- "Aves"

# Mammalia
mammalia_event <- extrair_exp_coefficients_simples(capture.output(summary(mammalia_event_100)))
mammalia_terminal <- extrair_exp_coefficients_simples(capture.output(summary(mammalia_terminal_100)))
mammalia_event$modelo <- "Geodiscovery"
mammalia_terminal$modelo <- "Completeness"
mammalia_combinado <- rbind(mammalia_event, mammalia_terminal)
mammalia_combinado$Class <- "Mammalia"

coef_plot <- bind_rows(
  amphibia_combinado,
  reptilia_combinado,
  aves_combinado,
  mammalia_combinado
)

# Ordem das variáveis
ordem_variaveis <- c(
  "BodySize",
  "Verticality",
  "Nocturnality",
  "RangeSize",
  "Elevation",
  "HumanDensity",
  "Latitude",
  "gdp_percentile_z"
)

# Ordem das classes
ordem_grupos <- c(
  "Amphibia",
  "Reptilia",
  "Aves",
  "Mammalia"
)

# Background
background_data <- expand.grid(
  Class = factor(ordem_grupos, levels = ordem_grupos),
  variable = factor(ordem_variaveis, levels = ordem_variaveis)
) %>%
  mutate(fill = cores[as.character(Class)])

# Organizar dados
coef_plot_ordenado <- coef_plot %>%
  mutate(
    variable = factor(variable, levels = ordem_variaveis),
    Class = factor(Class, levels = ordem_grupos),
    modelo = factor(
      modelo,
      levels = c("Geodiscovery", "Completeness")
    ),
    significant = lower > 1 | upper < 1
  )

background_data_ordenado <- background_data %>%
  mutate(Class = factor(Class, levels = ordem_grupos))

# Mesmo dodge para pontos e barras
pd <- position_dodge(width = 0.6, reverse = TRUE)

coefplot <- ggplot(
  coef_plot_ordenado,
  aes(
    x = hr,
    y = variable,
    colour = modelo,
    group = modelo
  )
) +
  # Background
  geom_rect(
    data = background_data_ordenado,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = fill
    ),
    inherit.aes = FALSE,
    alpha = 0.05,
    colour = NA
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "gray60",
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper,
      linetype = significant
    ),
    position = pd,
    height = 0.25,
    linewidth = 0.8
  ) +
  geom_point(
    aes(shape = significant),
    position = pd,
    size = 3.5,
    fill = "black"
  ) +
  scale_shape_manual(
    values = c(
      `TRUE` = 16,
      `FALSE` = 1
    )
  ) +
  scale_linetype_manual(
    values = c(
      `TRUE` = "solid",
      `FALSE` = "dashed"
    )
  ) +
  
  scale_color_manual(
    values = c(
      Completeness = "#696969",
      Geodiscovery = "#000000"
    )
  ) +
  
  scale_fill_identity() +
  facet_wrap(
    ~Class,
    ncol = 4
  ) +
  scale_y_discrete(
    limits = rev,
    labels = c(
      BodySize = "Body size",
      Verticality = "Verticality",
      Nocturnality = "Nocturnality",
      RangeSize = "Range size",
      Elevation = "Elevation",
      HumanDensity = "Human density",
      Latitude = "Latitude",
      gdp_percentile_z = "GDP-PPP"
    )
  ) +
  labs(
    x = "Hazard Ratio (HR)",
    y = "",
    colour = "Model"
  ) +
  
  geom_hline(
    data = data.frame(
      yintercept = seq(
        1.5,
        length(ordem_variaveis) - 0.5,
        by = 1
      )
    ),
    aes(yintercept = yintercept),
    inherit.aes = FALSE,
    colour = "grey70",
    linewidth = 0.2
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    axis.title.x = element_text(
      size = 15,
      #face = "bold",
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(
      size = 15,
      colour = "black"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    strip.text = element_text(
      #face = "bold",
      size = 15
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "grey80"
    ),
    panel.spacing.x = unit(0.3, "lines"),
    panel.spacing.y = unit(0.1, "lines"),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "grey80",
      linewidth = 0.4
    ),
    plot.margin = margin(15, 15, 15, 15)
  ); coefplot

ggsave(file.path("G:", "Meu Drive", "Artigos", "wallace_time", "Figures",
                 "Fig2.CoefPlot.pdf"),
       plot=coefplot, width=15, height=8, units="in",
       dpi = 'print', cairo_pdf)

## GOF plot ----
# Event
var.labels <- c(
  BodySize = "Body size",
  Verticality = "Verticality",
  Nocturnality = "Nocturnality",
  RangeSize = "Range size",
  Elevation = "Elevation",
  HumanDensity = "Human density",
  Latitude = "Latitude",
  gdp_percentile_z = "GDP-PPP"
)

par(mfrow=c(8,4),
    mar=c(1.2,1.7,1.2,1.1),   # margens internas (entre painéis)
    oma=c(1,1,0.5,1),         # margens externas (distância da página)
    mgp=c(1.7,0.3,0),
    xaxs="i",
    yaxs="i",
    cex.axis=1,
    cex.lab=5)
par(las=2)

for(i in 1:8){
  
  plot_gof_custom(
    gof.amphibia,
    variable = i,
    show.varname = TRUE,
    var.labels = var.labels,
    bg.col = cores[["Amphibia"]],
    titles = FALSE,
    tick.length = -0.3
  )
  
  plot_gof_custom(
    gof.reptilia,
    variable = i,
    show.varname = FALSE,
    var.labels = var.labels,
    bg.col = cores[["Reptilia"]],
    titles = FALSE,
    tick.length = -0.3
  )
  
  plot_gof_custom(
    gof.aves,
    variable = i,
    show.varname = FALSE,
    var.labels = var.labels,
    bg.col = cores[["Aves"]],
    titles = FALSE,
    tick.length = -0.3
  )
  
  plot_gof_custom(
    gof.mammalia,
    variable = i,
    show.varname = FALSE,
    var.labels = var.labels,
    bg.col = cores[["Mammalia"]],
    titles = FALSE,
    tick.length = -0.3
  )
}

# Terminal 
dev.off()
par(mfrow=c(8,4),
    mar=c(1.2,1.7,1.2,1.1),   # margens internas (entre painéis)
    oma=c(1,1,0.5,1),         # margens externas (distância da página)
    mgp=c(1.7,0.3,0),
    xaxs="i",
    yaxs="i",
    cex.axis=1,
    cex.lab=5)
par(las=2)

for(i in 1:8){
  
  plot_gof_custom(
    gof.terminal.amphibia,
    variable = i,
    show.varname = TRUE,
    var.labels = var.labels,
    bg.col = cores[["Amphibia"]],
    titles = FALSE,
    tick.length = -0.3
  )
  
  plot_gof_custom(
    gof.terminal.reptilia,
    variable = i,
    show.varname = FALSE,
    var.labels = var.labels,
    bg.col = cores[["Reptilia"]],
    titles = FALSE,
    tick.length = -0.3
  )
  
  plot_gof_custom(
    gof.terminal.aves,
    variable = i,
    show.varname = FALSE,
    var.labels = var.labels,
    bg.col = cores[["Aves"]],
    titles = FALSE,
    tick.length = -0.3
  )
  
  plot_gof_custom(
    gof.terminal.mammalia ,
    variable = i,
    show.varname = FALSE,
    var.labels = var.labels,
    bg.col = cores[["Mammalia"]],
    titles = FALSE,
    tick.length = -0.3
  )
}

## Time varying effects ----
load(
  file.path(
    local_directory,
    "02_data_analysis",
    "results",
    "aalen_event_models.RData")
  )

load(
  file.path(
    local_directory,
    "02_data_analysis",
    "results",
    "aalen_terminal_models.RData")
  )

dev.off()
par(mfrow = c(8, 4),
    mar = c(2.0, 2.0, 1.2, 1.1),  # AUMENTEI as margens inferior e esquerda
    oma = c(1.5, 1.5, 0.5, 0.5),  # Ajuste nas margens externas
    mgp = c(1.7, 0.5, 0),         # Ajuste para posição dos labels
    xaxs = "i",
    yaxs = "i",
    cex.axis = 0.8)               
par(las = 2)

# Definir sequências de tempo
tick_seq <- seq(1750, 2025, by = 25)
label_seq <- seq(1750, 2025, by = 50)

# Definir os grupos
grupos <- list(
  Amphibia = fit.amphibia.surv,
  Reptilia = fit.reptilia.surv,
  Aves = fit.aves.surv,
  Mammalia = fit.mammalia.surv
)

grupos <- list(
  Amphibia = fit.amphibia.completeness,
  Reptilia = fit.reptilia.completeness,
  Aves = fit.aves.completeness,
  Mammalia = fit.mammalia.completeness
)

# Plotar linha por linha
for(i in 1:8) {
  for(j in 1:4) {
    fit_atual <- grupos[[j]]
    nome_grupo <- names(grupos)[j]
    
    plot_surv_panel(fit_object = fit_atual,
                    comps = 3:10,
                    show_main = FALSE,
                    show_xlab = FALSE,
                    show_ylab = FALSE,
                    show_years = TRUE,
                    year_seq = label_seq,
                    tick_seq = tick_seq,
                    label_seq = label_seq,
                    coef_col = cores[[nome_grupo]],
                    coef_lwd = 2,
                    shade_alpha = 0.15,
                    var_index = i,
                    tick.length = -0.3,
                    xlim_range = c(1750, 2025))
  }
}

## Spatiotemporal map ----
richness_all <- richness_completude_amphibia %>%
  left_join(richness_completude_reptilia, by = "name_en") %>%
  left_join(richness_completude_aves, by = "name_en") %>%
  left_join(richness_completude_mammalia, by = "name_en") %>%
  mutate(
    polygons_richness = rowSums(select(., starts_with("polygons_richness")),
                                na.rm = TRUE),
    observed_richness = rowSums(select(., starts_with("observed_richness")),
                                na.rm = TRUE)
  ) %>%
  select("name_en", "polygons_richness", "observed_richness") %>%
  mutate(
    percent_wallace = observed_richness * 100 / polygons_richness,
    percent_wallace_prop = percent_wallace / 100
  ) 

event_total <- df_wallacean_100 %>%
  filter(event == 1) %>%
  count(name_en, name = "n_events")

event_recent <- df_wallacean_100 %>%
  mutate(year = year(date)) %>%
  filter(year > 1974) %>%
  filter(event == 1) %>%
  count(name_en, name = "n_events") %>%
  mutate(geo_taxa = n_events/50) %>%
  select(name_en, geo_taxa)

df_geo <- df_wallacean_100 %>%
  select(date, event, name_en, GlobalNorth) %>%
  filter(event == 1) %>%
  mutate(year = year(date)) %>%
  arrange(name_en, year) %>%
  left_join(event_total, by = "name_en") %>%
  left_join(richness_all, by = "name_en") %>%
  left_join(event_recent, by = "name_en") %>%
  group_by(name_en) %>%
  mutate(
    GlobalNorth = first(GlobalNorth),
    event_order = cumsum(event),
    prop = event_order / n_events,
    geo_taxa = first(geo_taxa)
  )

milestones <- df_geo %>%
  group_by(name_en) %>%
  summarise(
    GlobalNorth = first(GlobalNorth),
    PercentWallace = first(percent_wallace),
    PropWallace = first(percent_wallace_prop),
    geo_taxa = first(geo_taxa),
    year25 = min(year[prop >= 0.25]),
    year50 = min(year[prop >= 0.50]),
    year75 = min(year[prop >= 0.75]),
    year90 = min(year[prop >= 0.90])
  ) 

geo_plot <- geographic_shape_data %>%
  left_join(milestones %>%
              mutate(
                IQR = year75 - year25
              ), by = "name_en")

# a acumulação foi rápida?
# Quando 75% das geodescobertas ja haviam ocorrido?
plot_time <- plot_bivariate_map(
  geo_plot,
  "IQR",
  "year75",
  title = NULL,
  x_label = "Timing of geodiscovery",
  y_label = "Duration of geodiscovery"); plot_time

plot_spatialtime <- plot_bivariate_map(
  geo_plot,
  "PercentWallace",
  "year75",
  title = NULL,
  x_label = "Inventory completeness (%)",
  y_label = "Year when 75% of\ngeodiscoveries were accumulated"); plot_spatialtime

plot_spatialiqr <- plot_bivariate_map(
  geo_plot,
  "PercentWallace",
  "IQR",
  title = NULL,
  x_label = "Duration of geodiscovery",
  y_label = "Inventory completeness (%)"); plot_spatialiqr

plot_rate <- plot_bivariate_map(
  geo_plot,
  "PercentWallace",
  "geo_taxa",
  title = NULL,
  x_label = "Inventory completeness (%)",
  y_label = "Geographic expansion rate"); plot_rate

ggsave(filename= "FigureS1.pdf",
       plot= plot_time,
       width=7.09,
       height=4.73,
       units="in",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf)

ggsave(filename= "FigureS2.pdf",
       plot= plot_spatialtime,
       width=7.09,
       height=4.73,
       units="in",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf)

ggsave(filename= "FigureS3.pdf",
       plot= plot_spatialiqr,
       width=7.09,
       height=4.73,
       units="in",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf)

ggsave(filename= "FigureS4.pdf",
       plot= plot_rate,
       width=7.09,
       height=4.73,
       units="in",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf)

# Explorando mapa tricolor e PCA
library(tricolore)
library(ggplot2)
library(sf)
library(ggtern)
library(ggtern)
library(grid)

# 1. Calcular a cor para cada região com base nas 3 variáveis
data <- geo_plot %>%
  mutate(
    KnowledgeAge = normalize(2025 - year75), 
    year75 = as.numeric(normalize(year75)),
    PercentWallace = as.numeric(normalize(PercentWallace)),
    IQR = as.numeric(normalize(IQR))) %>% 
  select(name_en, ISO3, 
         GlobalNorth.y, year75, IQR, PercentWallace, KnowledgeAge) %>%
  remove_missing() %>%
  st_transform(6933)  # ou o CRS que você preferir

# analises
vars <- data %>%
  st_drop_geometry() %>%
  select(ISO3, GlobalNorth.y, PercentWallace, IQR, year75)

cor(vars[,-1])

# PCA
pca <- prcomp(vars[, c("KnowledgeAge",
                       "PercentWallace",
                       "IQR")],
              scale. = TRUE)

scores <- as.data.frame(pca$x) %>%
  bind_cols(vars[, c("ISO3", "GlobalNorth.y")])

ggplot(scores,
       aes(PC1, PC2,
           label = ISO3,
           color = factor(GlobalNorth.y))) +
  geom_point(size = 2.5) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_color_manual(values = c("red", "blue")) +
  theme_classic()


# Mapa tricolor
tric <- Tricolore(data, 
                  p1 = 'year75', 
                  p2 = 'PercentWallace', 
                  p3 = 'IQR',
                  breaks = Inf)

# Adicionar a cor calculada ao seu dataframe
data$rgb <- tric$rgb

# Bounding box para borda
bbox <- st_bbox(data)
bbox_poly <- st_as_sfc(bbox)

plot <- ggplot(data) +
  geom_sf(aes(fill = rgb, geometry = geometry), 
          color = "black", 
          size = 0.15) +
  geom_sf(data = bbox_poly, 
          fill = NA, 
          color = "black", 
          linewidth = 0.3) +
  scale_fill_identity() +
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
    axis.text = element_text(size = 6, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    axis.ticks.length = unit(2, "mm"),
    panel.grid = element_line(color = "grey80", linewidth = 0.2),
    plot.margin = margin(5, 15, 5, 5)
  )

# 1. Aumentar o tamanho da legenda
key_grob <- ggplotGrob(tric$key)
key_grob$widths <- key_grob$widths * 2.0   
key_grob$heights <- key_grob$heights * 2.0 

# 2. Descobrir os limites do seu mapa para posicionar a legenda
# Você precisa saber os valores reais de x e y do seu mapa
bbox <- st_bbox(data)  # Pega os limites geográficos

# 3. Posicionar a legenda usando coordenadas REAIS do mapa
# Use valores baseados no bbox do seu mapa
plot_final <- plot + annotation_custom(
  ggplotGrob(
    tric$key),
  xmin = bbox$xmin + (bbox$xmax - bbox$xmin) * 0.02,  # 2% da esquerda
  xmax = bbox$xmin + (bbox$xmax - bbox$xmin) * 0.30,  # 30% da esquerda
  ymin = bbox$ymin + (bbox$ymax - bbox$ymin) * 0.12,  # 2% da base
  ymax = bbox$ymin + (bbox$ymax - bbox$ymin) * 0.45   # 35% da base
); plot_final

ggsave(filename= "FigureTest.pdf",
       plot= plot_final,
       width=16.09,
       height=9.73,
       units="in",
       dpi = 300,
       bg="white", 
       limitsize=F,
       device = cairo_pdf)
