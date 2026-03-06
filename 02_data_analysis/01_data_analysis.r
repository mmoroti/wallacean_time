# Packages and functions ----
library(tidyverse)
library(mets)
library(timereg)

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
                            year_seq = seq(1750, 2025, by = 25),  # <<< NOVO
                            ...) {
  
  if (is.null(type)) {
    if (x$type == "prop") type <- "time"
    if (x$type == "modelmatrix") type <- "modelmatrix"
    if (x$type == "Zmodelmatrix") type <- "z"
  }
  
  xlab_final <- if (show.xlab) xlab else ""
  ylab_final <- if (show.ylab) ylab else ""
  
  if (type == "time" || type == "modelmatrix") {
    
    p <- ncol(x$score)
    
    for (i in 1:p) {
      
      simU <- x$simUt[, (0:49) * p + i]
      rsU <- max(abs(simU))
      rsU <- max(rsU, abs(x$score[, i]))
      
      # FIXAR LIMITES DO EIXO X
      plot(x$jumptimes, x$score[, i],
           type = "n",
           ylim = c(-rsU, rsU),
           xlim = range(year_seq),   # <<< garante mesmo eixo
           xlab = xlab_final,
           ylab = ylab_final,
           cex.lab = cex.lab,
           cex.axis = cex.axis,
           xaxt = "n",
           yaxt = "n",
           ...)
    
      if (show.xaxis || show.xticks) {
        axis(1,
             at = year_seq,
             labels = if (show.xaxis) year_seq else rep("", length(year_seq)),
             tick = show.xticks,
             cex.axis = cex.axis)
      }

      if (show.yaxis || show.yticks) {
        axis(2,
             labels = show.yaxis,
             tick = show.yticks,
             cex.axis = cex.axis)
      }
      
      box()
      
      if (titles) title(main = rownames(x$res)[i])
      
      matlines(x$jumptimes, simU,
               type = "s",
               lwd = 0.3,
               col = if (!is.null(bg.col)) bg.col else col)
      
      lines(x$jumptimes, x$score[, i],
            type = "s",
            lwd = 1.5)
    }
    
  } else {
    
    for (i in 1:length(x$Zres)) {
      
      xr <- x$Zres[[i]]
      obsz <- c(tail(xr$score, 1))
      times <- xr$xaxs
      
      rsU <- max(max(abs(obsz)),
                 max(abs(xr$simUtlast[1:50, ])))
      
      plot(times, obsz,
           type = "n",
           ylim = c(-rsU, rsU),
           xlim = range(year_seq),   # <<< garante mesmo eixo
           xlab = xlab_final,
           ylab = ylab_final,
           cex.lab = cex.lab,
           xaxt = "n",
           yaxt = "n",
           ...)
      
      if (show.xaxis || show.xticks) {
        axis(1,
             at = year_seq,
             labels = if (show.xaxis) year_seq else rep("", length(year_seq)),
             tick = show.xticks,
             cex.axis = cex.axis)
      }
      
      if (show.yaxis || show.yticks) {
        axis(2,
             labels = show.yaxis,
             tick = show.yticks,
             cex.axis = cex.axis)
      }
      
      box()
      
      if (titles) title(main = rownames(x$res)[i])
      
      matlines(times,
               t(xr$simUtlast[1:50, ]),
               type = "l",
               lwd = 0.3,
               col = if (!is.null(bg.col)) bg.col else col)
      
      lines(times,
            obsz,
            lwd = 2,
            col = 1)
    }
  }
}

plot_surv_panel <- function(fit_object,
                            comps = 3:10,
                            show_main = TRUE,
                            show_xlab = TRUE,
                            show_ylab = TRUE,
                            show_years = TRUE,       # controla labels, não ticks
                            year_seq = seq(1750, 2025, by = 50),
                            coef_col = "blue",
                            coef_lwd = 2,
                            shade_alpha = 0.15) {
  
  for(i in seq_along(comps)) {
    comp <- comps[i]
    
    # extrai dados do componente
    time <- fit_object$cum[, 1]
    est <- fit_object$cum[, comp]
    sd <- sqrt(fit_object$robvar.cum[, comp])
    upper <- est + 2*sd
    lower <- est - 2*sd
    
    # labels
    xlab <- ifelse(show_xlab, "Year", "")
    ylab <- ifelse(show_ylab, paste0("Coef ", colnames(fit_object$cum)[comp]), "")
    main <- ifelse(show_main, colnames(fit_object$cum)[comp], "")
    
    # plot vazio com limites
    ylim <- range(c(lower, upper, 0))
    plot(time, est, type = "n", xlab = xlab, ylab = ylab,
         main = main, ylim = ylim, xaxt = "n")
    
    # linha no zero
    abline(h = 0, col = "gray50", lty = 2)
    
    # shade
    polygon(c(time, rev(time)),
            c(upper, rev(lower)),
            col = adjustcolor(coef_col, alpha.f = shade_alpha),
            border = NA)
    
    # linha central
    lines(time, est, col = coef_col, lwd = coef_lwd)
    
    # grid
    abline(v = year_seq, col = "gray80", lty = 3, lwd = 0.5)
    
    # eixo X: controla só os labels
    axis(1, at = year_seq,
         labels = if(show_years) year_seq else rep("", length(year_seq)),
         las = 2, cex.axis = 1.0)
  }
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

# Preparing data ----
df_wallacean_100 <- df_wallacean_100 %>%
  mutate(status = event + WallaceCompletude) %>%
  left_join(
    geographic_shape_data %>%
      select(GlobalNorth, Georegion, name_en) %>%
      sf::st_drop_geometry(), by = "name_en"
  )

df_amphibia_100 <- df_wallacean_100 %>%
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
tetrapods_100 <- recurrentMarginal(Event(t.start.year, t.stop.year, status==1)~
                                     strata(GlobalNorth)+
                                     cluster(speciesKey),
                                   data=df_wallacean_100,
                                   cause=1,
                                   death.code=2)
save(
  tetrapods_100,
  file = file.path(
    local_directory, "02_data_analysis", "tetrapods_model.RData")
) # plot in script "02_exploratory_analysis.R"

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
                            data=df_reptilia_100)
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
                               data=df_reptilia_100)
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

## Cumulative hazard plots ----
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

cores <- c("Amphibia" = "#75C6FF",
           "Reptilia" = "#86A94D",
           "Aves" = "#7F5C3B",
           "Mammalia" = "#D66FFF")

dev.off()
par(mfrow=c(1,2),
    cex.lab = 1.5, 
    cex.axis = 1.3)  # Aumentar tamanho dos labels dos eixos

# Event
plot(amphibia_event_100, se = TRUE, col = cores["Amphibia"],
     legend = FALSE,
     ylim = c(0, 35),
     xlim = c(1750, 2025),
     xlab = "",  # Omitir título do eixo x
     cex.main = 1.5)  # Aumentar tamanho do título
plot(reptilia_event_100, se = TRUE, add = TRUE, col = cores["Reptilia"])
plot(aves_event_100, se = TRUE, add = TRUE, col = cores["Aves"])
plot(mammalia_event_100, se = TRUE, add = TRUE, col = cores["Mammalia"])
title(main = "Geodiscovery", line = 0.5, cex.main = 1.5)  # Ajustar posição do título
legend("topleft", legend = names(cores), col = unname(cores), 
       lwd = 2, bty = "n", title = "", cex = 1.2)  # Aumentar tamanho da legenda

anos <- seq(1750, 2025, by = 50)
abline(v = anos, col = "gray70", lty = 2, lwd = 1)

# Terminal event
plot(amphibia_terminal_100, se = TRUE, col = cores["Amphibia"], 
     legend = FALSE,
     ylim = c(0, 0.25),
     xlim = c(1750, 2025),
     xlab = "",  # Omitir título do eixo x
     cex.main = 1.5)
plot(reptilia_terminal_100, se = TRUE, add = TRUE, col = cores["Reptilia"])
plot(aves_terminal_100, se = TRUE, add = TRUE, col = cores["Aves"])
plot(mammalia_terminal_100, se = TRUE, add = TRUE, col = cores["Mammalia"])
title(main = "Completeness", line = 0.5, cex.main = 1.5)

anos <- seq(1750, 2025, by = 50)
abline(v = anos, col = "gray70", lty = 2, lwd = 1)

# Mean events
amphibia_rm <- recurrentMarginalPhreg(amphibia_event_100,
                                      amphibia_terminal_100)
reptilia_rm <- recurrentMarginalPhreg(reptilia_event_100,
                                      reptilia_terminal_100)
aves_rm <- recurrentMarginalPhreg(aves_event_100,
                                  aves_terminal_100)
mammalia_rm <- recurrentMarginalPhreg(mammalia_event_100,
                                      mammalia_terminal_100)

# Mean between event + terminal
plot(amphibia_rm, se = TRUE, col = cores["Amphibia"], 
     legend = FALSE,
     ylim = c(0, 35),
     xlim = c(1750, 2025),
     xlab = "",  # Omitir título do eixo x
     cex.main = 1.5)
plot(reptilia_rm, se = TRUE, add = TRUE, col = cores["Reptilia"])
plot(aves_rm, se = FALSE, add = TRUE, col = cores["Aves"])
plot(mammalia_rm, se = TRUE, add = TRUE, col = cores["Mammalia"])
title(main = "Geodiscovery + Completeness", line = 0.5, cex.main = 1.5)

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

# Definir a ordem desejada das variáveis
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

# Definir a ordem desejada dos grupos (classes)
ordem_grupos <- c("Amphibia", "Reptilia", "Aves", "Mammalia")

# Criar background_data com uma cor por grupo
background_data <- expand.grid(
  Class = factor(ordem_grupos, levels = ordem_grupos),
  variable = factor(ordem_variaveis, levels = ordem_variaveis)
) %>%
  mutate(
    # Mapear a cor baseada no grupo
    fill = cores[as.character(Class)]
  )

# Converter ambas as colunas para fatores com a ordem especificada
coef_plot_ordenado <- coef_plot %>%
  mutate(
    variable = factor(variable, levels = ordem_variaveis),
    Class = factor(Class, levels = ordem_grupos)
  )

# Reordenar também o background_data para corresponder
background_data_ordenado <- background_data %>%
  mutate(Class = factor(Class, levels = ordem_grupos))

# Agora usar o dataframe ordenado no gráfico
coefplot <- ggplot(coef_plot_ordenado,
                   aes(x = hr, y = variable, color = modelo)) +
  # Background colorido (usar o background_data ordenado)
  geom_rect(
    data = background_data_ordenado,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.05,
    color = NA
  ) +
  # Linha vertical de referência
  geom_vline(xintercept = 1, linetype = "dashed", 
             color = "black", linewidth = 0.6, alpha = 0.8) +
  # Barras de erro
  geom_errorbarh(
    aes(xmin = lower, xmax = upper), 
    position = position_dodge(width = 0.6),
    height = 0.25,
    linewidth = 1.2,
    alpha = 1
  ) +
  # Pontos
  geom_point(
    position = position_dodge(width = 0.6), 
    size = 3.5,
    shape = 16
  ) +
  # Labels
  labs(
    x = "Hazard Ratio (HR)", 
    y = "",
    color = "Model"
  ) +
  # Escalas
  scale_color_manual(
    values = c("Completeness" = "#A23B72", "Geodiscovery" = "#2E86AB"),
    labels = c("Completeness", "Geodiscovery")
  ) +
  scale_fill_identity() +
  # Facets - agora os grupos aparecerão na ordem especificada
  facet_wrap(~ Class, ncol = 4) +
  # Reverter a ordem do eixo Y para que BodyLength_mm fique no topo
  scale_y_discrete(limits = rev) +
  # Tema
  theme_minimal() +
  theme(
    # Texto
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 15, color = "black"),
    
    # Legendas
    legend.position = "top",  # MUDADO: de "bottom" para "top"
    legend.title = element_blank(),  # MUDADO: remove o título "Model"
    legend.text = element_text(size = 18),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(b = 0, t = 0),  # ADICIONADO: ajusta margens
    
    # Strips (títulos dos facets) - agora na ordem correta
    strip.text = element_text(
      face = "bold", 
      size = 15,
      color = "black",
      margin = margin(b = 8, t = 8)
    ),
    strip.background = element_rect(
      fill = "white",
      color = "grey80",
      linewidth = 0.5
    ),
    
    # Painéis - REDUZIR espaçamento
    panel.spacing.x = unit(0.3, "lines"),  # Espaço horizontal entre grupos
    panel.spacing.y = unit(0.1, "lines"),  # Espaço vertical (se tiver múltiplas linhas)
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      fill = NA, 
      color = "grey80", 
      linewidth = 0.4
    ),
    
    # Plot geral
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  # Adicionar linhas divisórias horizontais entre variáveis
  geom_hline(
    data = data.frame(
      yintercept = seq(1.5, length(ordem_variaveis) - 0.5, by = 1)
    ),
    aes(yintercept = yintercept),
    color = "grey70",
    linewidth = 0.2,
    linetype = "solid"
  ); coefplot

ggsave(file.path("G:", "Meu Drive", "Artigos", "wallace_time", "Figures",
                 "Fig2.CoefPlot.pdf"),
       plot=coefplot, width=15, height=8, units="in",
       dpi = 'print', cairo_pdf)

## GOF plot ----
# Event
dev.off()
par(mfrow=c(4,8),
    mar=c(2.5, 3, 1.5, 0.5),
    oma=c(0,0,0,0),
    mgp=c(1.7,0.6,0),
    xaxs="i",
    yaxs="i",
    cex.axis=1,
    cex.lab=5)
par(las=2)

plot_gof_custom(gof.amphibia, 
                titles = TRUE,
                bg.col = "#75C6FF",
                year_seq = seq(1750, 2025, by = 25),
                show.xaxis = TRUE,
                show.yaxis = TRUE,
                show.xticks = TRUE,
                show.yticks = TRUE)
plot_gof_custom(gof.reptilia,
                titles = FALSE,
                bg.col = "#86A94D",
                year_seq = seq(1750, 2025, by = 25),
                show.xaxis = TRUE,
                show.yaxis = TRUE,
                show.xticks = TRUE,
                show.yticks = TRUE)
plot_gof_custom(gof.aves,
                titles = FALSE,
                bg.col = "#7F5C3B",
                year_seq = seq(1750, 2025, by = 25),
                show.xaxis = TRUE,
                show.yaxis = TRUE,
                show.xticks = TRUE,
                show.yticks = TRUE)
plot_gof_custom(gof.mammalia,
                titles = FALSE,
                bg.col = "#D66FFF",
                year_seq = seq(1750, 2025, by = 25),
                show.xaxis = TRUE,
                show.yaxis = TRUE,
                show.xticks = TRUE,
                show.yticks = TRUE)

# Terminal 
dev.off()
par(mfrow=c(4,8),
    mar=c(2.5,2,1.5,0.5),
    oma=c(0,0,0,0),
    mgp=c(1.7,0.6,0),
    xaxs="i",
    yaxs="i",
    cex.axis=1,
    cex.lab=5)
par(las=2)

plot_gof_custom(gof.terminal.amphibia,
                titles = TRUE,
                year_seq = seq(1900, 2025, by = 25),
                bg.col = "#75C6FF")
plot_gof_custom(gof.terminal.reptilia,
                titles = FALSE,
                year_seq = seq(1900, 2025, by = 25),
                bg.col = "#86A94D")
plot_gof_custom(gof.terminal.aves,
                titles = FALSE,
                year_seq = seq(1850, 2025, by = 25),
                bg.col = "#7F5C3B")
plot_gof_custom(gof.terminal.mammalia, 
                titles = FALSE,
                year_seq = seq(1900, 2025, by = 25),
                bg.col = "#D66FFF")

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
par(mfrow=c(4,8),
    mar=c(2.5,2,1.5,0.5),
    oma=c(0,0,0,0),
    mgp=c(1.7,0.6,0),
    xaxs="i",
    yaxs="i",
    cex.axis=1,
    cex.lab=5)

plot_surv_panel(fit.amphibia.surv,
                show_main = TRUE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#2DA5E8",
                shade_alpha = 0.2)

plot_surv_panel(fit.reptilia.surv,
                show_main = FALSE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#86A94D",
                shade_alpha = 0.2)

plot_surv_panel(fit.aves.surv,
                show_main = FALSE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#7F5C3B",
                shade_alpha = 0.2)

plot_surv_panel(fit.mammalia.surv,
                show_main = FALSE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#D66FFF",
                shade_alpha = 0.2) 

# Completeness
dev.off()
par(mfrow=c(4,8),
    mar=c(2.5,2,1.5,0.5),
    oma=c(0,0,0,0),
    mgp=c(1.7,0.6,0),
    xaxs="i",
    yaxs="i",
    cex.axis=1,
    cex.lab=5)      

plot_surv_panel(fit.amphibia.completeness,
                show_main = TRUE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#2DA5E8",
                shade_alpha = 0.2)

plot_surv_panel(fit.reptilia.completeness,
                show_main = FALSE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#86A94D",
                shade_alpha = 0.2)

plot_surv_panel(fit.aves.completeness,
                show_main = FALSE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#7F5C3B",
                shade_alpha = 0.2)

plot_surv_panel(fit.mammalia.completeness,
                show_main = FALSE,
                show_xlab = FALSE,
                show_ylab = FALSE,
                show_years = TRUE,
                coef_col = "#D66FFF",
                shade_alpha = 0.2) 

# Hist plot events
df_wallacean_100 %>%
  filter(status == 1) %>%  # só eventos
  ggplot(aes(x = year(date), fill = factor(GlobalNorth))) +
  geom_histogram(binwidth = 1, color = "white",
                 alpha = 0.6) +
  # Linhas históricas
  geom_vline(xintercept = 1914, linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = 1939, linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = 2002, linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = 2008, linetype = "dashed", size = 0.7) +
  scale_fill_manual(
    values = c("0" = "#c40e3e", "1" = "#2DA5E8"),
    labels = c("0" = "Global South", "1" = "Global North"),
    name = ""
  ) +
  labs(x = "Year", y = "Event frequency") +
  theme_bw()

# Define a ordem desejada para sobreposição
class_order <- c("Aves", "Mammalia", "Reptilia", "Amphibia")

df_wallacean_100 %>% 
  filter(status == 1) %>%  # só eventos
  mutate(Class = factor(Class, levels = class_order)) %>%
  ggplot(aes(x = year(date), fill = factor(Class))) +
  geom_histogram(binwidth = 1, color = "white", alpha = 1) +
  scale_fill_manual(
    values = c("Amphibia" = "#75C6FF", "Reptilia" = "#86A94D",
               "Aves" = "#7F5C3B", "Mammalia" = "#D66FFF"),
    name = ""
  ) +
  labs(x = "Year", y = "Geodiscovery frequency") +
  theme_bw()

df_wallacean_100 %>% 
  filter(status == 2) %>%  # só eventos
  mutate(Class = factor(Class, levels = class_order)) %>%
  ggplot(aes(x = year(date), fill = factor(Class))) +
  geom_histogram(binwidth = 1, color = "white", alpha = 1) +
  scale_fill_manual(
    values = c("Amphibia" = "#75C6FF", "Reptilia" = "#86A94D",
               "Aves" = "#7F5C3B", "Mammalia" = "#D66FFF"),
    name = ""
  ) +
  labs(x = "Year", y = "Completeness frequency") +
  theme_bw()

