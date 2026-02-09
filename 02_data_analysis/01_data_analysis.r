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
                            xlab = "",          # label do eixo X
                            ylab = "",          # label do eixo Y
                            cex.lab = 1,        # tamanho dos labels
                            cex.axis = 1,       # tamanho dos números dos eixos
                            show.xlab = TRUE,   # mostrar label do eixo X?
                            show.ylab = TRUE,   # mostrar label do eixo Y?
                            show.xaxis = TRUE,  # mostrar números do eixo X?
                            show.yaxis = TRUE,  # mostrar números do eixo Y?
                            show.xticks = TRUE, # mostrar ticks do eixo X?
                            show.yticks = TRUE, # mostrar ticks do eixo Y?
                            ...) {
  if (is.null(type)) {
    if (x$type == "prop") type <- "time"
    if (x$type == "modelmatrix") type <- "modelmatrix"
    if (x$type == "Zmodelmatrix") type <- "z"
  }
  
  # Preparar parâmetros dos eixos
  xlab_final <- if (show.xlab) xlab else ""
  ylab_final <- if (show.ylab) ylab else ""
  
  if (type == "time" || type == "modelmatrix") {
    p <- ncol(x$score)
    for (i in 1:p) {
      simU <- x$simUt[, (0:49) * p + i]
      rsU <- max(abs(simU))
      rsU <- max(rsU, abs(x$score[, i]))
      
      # Plotar sem números/ticks primeiro
      plot(x$jumptimes, x$score[, i], type = "n", ylim = c(-rsU, rsU),
           xlab = xlab_final, ylab = ylab_final,
           cex.lab = cex.lab, cex.axis = cex.axis,
           xaxt = "n", yaxt = "n",  # Sem eixos inicialmente
           ...)
      
      # Adicionar eixos customizados
      if (show.xaxis || show.xticks) {
        axis(1, labels = show.xaxis, tick = show.xticks, cex.axis = cex.axis)
      }
      
      if (show.yaxis || show.yticks) {
        axis(2, labels = show.yaxis, tick = show.yticks, cex.axis = cex.axis)
      }
      
      # Adicionar grade (opcional)
      box()
      
      if (titles) title(main = rownames(x$res)[i])
      
      # Adicionar as linhas
      matlines(x$jumptimes, simU, type = "s", lwd = 0.3,
               col = if (!is.null(bg.col)) bg.col else col)
      lines(x$jumptimes, x$score[, i], type = "s", lwd = 1.5)
    }
  } else {
    for (i in 1:length(x$Zres)) {
      xr <- x$Zres[[i]]
      obsz <- c(tail(xr$score, 1))
      times <- xr$xaxs
      rsU <- max(max(abs(obsz)), max(abs(xr$simUtlast[1:50, ])))
      
      # Plotar sem eixos inicialmente
      plot(times, obsz, type = "n", ylim = c(-rsU, rsU),
           xlab = xlab_final, ylab = ylab_final,
           cex.lab = cex.lab,
           xaxt = "n", yaxt = "n",  # Sem eixos inicialmente
           ...)
      
      # Adicionar eixos customizados
      if (show.xaxis || show.xticks) {
        axis(1, labels = show.xaxis, tick = show.xticks, cex.axis = cex.axis)
      }
      
      if (show.yaxis || show.yticks) {
        axis(2, labels = show.yaxis, tick = show.yticks, cex.axis = cex.axis)
      }
      
      # Adicionar grade
      box()
      
      if (titles) title(main = rownames(x$res)[i])
      
      # Adicionar as linhas
      matlines(times, t(xr$simUtlast[1:50, ]), type = "l", lwd = 0.3,
               col = if (!is.null(bg.col)) bg.col else col)
      lines(times, obsz, lwd = 2, col = 1)
    }
  }
}

summary_aalen <- function(fit, Class, Type) {
  # Captura a saída do summary
  summary_lines <- capture.output(summary(fit))
  
  # Encontra o início da seção de termos paramétricos
  start_idx <- grep("Parametric terms", summary_lines)
  if (length(start_idx) == 0) stop("Seção 'Parametric terms' não encontrada.")
  
  # Cabeçalho da tabela
  header_line <- summary_lines[start_idx + 1]
  
  # Encontra o fim da tabela (linha em branco ou início do "Call")
  table_start <- start_idx + 2
  call_idx <- grep("^\\s*Call", summary_lines)
  table_end <- if (length(call_idx) > 0) call_idx[1] - 2 else which(summary_lines == "")[1] - 2
  table_lines <- summary_lines[table_start:table_end]
  
  # Junta cabeçalho e dados
  table_text <- c(header_line, table_lines)
  
  # Separa os termos (primeira palavra de cada linha, exceto o cabeçalho)
  term_names <- sub("^\\s*(\\S+).*", "\\1", table_text[-1])  # ignora a primeira linha (header)
  # Remove o termo da linha antes de ler a tabela
  clean_lines <- gsub("^\\s*\\S+\\s+", "", table_text[-1])
  # Junta novamente: header + linhas limpas
  table_clean <- c(table_text[1], clean_lines)
  
  # Lê a tabela
  param_df <- read.table(text = table_text,
                         header = FALSE,
                         fill = TRUE,
                         check.names = FALSE)
  
  param_df <- param_df %>%
    rename(term = V1,
           coef = V2,
           SE = V3,
           'Robust SE' = V4,
           z = V5,
           'P-val' = V6) %>%
           #'lower2.5%' = V7
           #'upper97.5%' = V8
    filter(term != 'Coef.') %>%
    mutate(Class = Class,
           Type = Type)
  
  return(param_df)
}

plot_cumulative_events <- function(fits, 
                                   cores = NULL,
                                   title = NULL,
                                   legend = TRUE,
                                   xlab = NULL, 
                                   ylab = NULL) {
  # 2. Extrai times, mu e se.mu
  times_list <- lapply(fits, function(f) f$times)
  mu_list    <- lapply(fits, function(f) f$mu)
  se_list    <- lapply(fits, function(f) f$se.mu)
  
  # 3. xlim unificado
  xrange <- range(unlist(times_list))
  
  # 4. ylim unificado
  all_lower <- unlist(mapply(function(mu, se) mu - 1.96*se, mu_list, se_list))
  all_upper <- unlist(mapply(function(mu, se) mu + 1.96*se, mu_list, se_list))
  yrange    <- range(all_lower, all_upper, na.rm = TRUE)
  
  # 5. Define cores padrão se não fornecidas
  if (is.null(cores)) {
    cores <- c("#4FCFF5", "#A8896A", "#F58C25", "#F54952")
  }
  
  # 6. Abre o gráfico vazio
  par(cex.axis = 1.5, cex.lab = 1.5)
  
  first_grp <- names(fits)[1]
  plot(times_list[[first_grp]], mu_list[[first_grp]], type = "n",
       xlim = xrange, ylim = yrange,
       xlab = xlab, ylab = ylab, main = title)

  # 7. Desenha cada curva com sombra de IC
  i <- 1
  for(g in names(fits)) {
    col     <- cores[(i-1) %% length(cores) + 1]  # Recycle colors if needed
    t       <- times_list[[g]]
    mu      <- mu_list[[g]]
    se      <- se_list[[g]]
    lower   <- mu - 1.96*se
    upper   <- mu + 1.96*se
    
    # desenha a sombra (polígono) do IC
    polygon(
      x = c(t, rev(t)),
      y = c(upper, rev(lower)),
      col = adjustcolor(col, alpha.f = 0.2),
      border = NA
    )
    
    # desenha a linha central
    lines(t, mu, type="s", col=col, lty=1, lwd=2)
    i <- i + 1
  }
  
  # 8. Adiciona legenda se solicitado
  if (legend) {
    legend("topleft",
           legend = names(fits),
           col    = cores[1:length(fits)],
           lty    = 1,
           lwd    = 2,
           bty    = "n")
  }
}

# Load data ----
load(file = file.path(
  local_directory,
  "02_data_analysis",
  "data_wallacean_time.RData")
)

load(file = file.path(
  "00_raw_data",
  "adm_unit_global.RData")
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
      st_drop_geometry(), by = "name_en"
  )

df_amphibia_100 <- df_wallacean_100 %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Amphibia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensity+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation))) 

df_reptilia_100 <- df_wallacean_100 %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Reptilia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensity+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation))) 

df_aves_100 <- df_wallacean_100 %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Aves") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensity+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize)),
         Elevation = scale(log(Elevation)))  

df_mammalia_100 <- df_wallacean_100 %>%
  mutate(status = event + WallaceCompletude) %>%
  filter(Class == "Mammalia") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         HumanDensity = scale(log(HumanDensity+1)),
         Latitude = scale(abs(Latitude)),
         RangeSize = scale(log(RangeSize))) 

# Recurrent event model ----
tetrapods_event_100 <- phreg(Surv(t.start.year, t.stop.year, status==1)~
                              cluster(speciesKey)+
                              strata(GlobalNorth),
                            data= df_wallacean_100)
plot(tetrapods_event_100, se = TRUE)

tetrapods_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~
                                 cluster(speciesKey)+
                                 strata(GlobalNorth),
                               data=df_wallacean_100)
plot(tetrapods_terminal_100, se = TRUE)

## Amphibia ----
# Cox model
amphibia_event_100 <- phreg(Surv(t.start.year, t.stop.year, status==1)~
                              cluster(speciesKey)+
                              Verticality +
                              Nocturnality +
                              BodyLength_mm +
                              Latitude +
                              HumanDensity +
                              Elevation +
                              RangeSize,
                            data=df_amphibia_100)

summary(amphibia_event_100)
plot(amphibia_event_100, se = TRUE)

gof.amphibia <- gof(amphibia_event_100)

amphibia_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~
                                 cluster(speciesKey)+
                                 Verticality +
                                 Nocturnality +
                                 BodyLength_mm +
                                 Latitude +
                                 HumanDensity +
                                 Elevation +
                                 RangeSize,
                               data=df_amphibia_100)
summary(amphibia_terminal_100)
gof.terminal.amphibia <- gof(amphibia_terminal_100)

# Aalen model
fit.amphibia.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
    cluster(speciesKey),
  data= df_amphibia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.amphibia.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
    cluster(speciesKey),
  data= df_amphibia_100,
  start.time = 1850,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.amphibia <- recurrent.marginal.mean(
  fit.amphibia.surv, 
  fit.amphibia.completeness
)

## Reptilia ----
# Cox model
reptilia_event_100 <- phreg(Surv(t.start.year,
                                 t.stop.year,
                                 status==1)~
                              cluster(speciesKey)+
                              Verticality +
                              Nocturnality +
                              BodyLength_mm +
                              Latitude +
                              HumanDensity +
                              Elevation +
                              RangeSize,
                            data=df_reptilia_100)
summary(reptilia_event_100)
gof.reptilia <- gof(reptilia_event_100)

reptilia_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~
                                 cluster(speciesKey)+
                                 Verticality +
                                 Nocturnality +
                                 BodyLength_mm +
                                 Latitude +
                                 HumanDensity +
                                 Elevation +
                                 RangeSize,
                               data=df_reptilia_100)
summary(reptilia_terminal_100)
gof.terminal.reptilia <- gof(reptilia_terminal_100)

# Aalen model
fit.reptilia.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
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
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
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
  fit.reptilia.death
)

## Aves ----
aves_event_100 <- phreg(Surv(t.start.year,
                             t.stop.year,
                             status==1)~
                          cluster(speciesKey)+
                          Verticality +
                          Nocturnality +
                          BodyLength_mm +
                          Latitude +
                          HumanDensity +
                          Elevation +
                          RangeSize,
                        data=df_aves_100)
summary(aves_event_100)
gof.aves <- gof(aves_event_100)

aves_terminal_100 <- phreg(Surv(t.start.year,
                                t.stop.year,
                                status==2)~
                             cluster(speciesKey)+
                             Verticality +
                             Nocturnality +
                             BodyLength_mm +
                             Latitude +
                             HumanDensity +
                             Elevation +
                             RangeSize+
                             strata(GlobalNorth),
                           data=df_aves_100)
summary(aves_terminal_100)
gof.terminal.aves <- gof(aves_terminal_100)

# Aalen model
fit.aves.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
    cluster(speciesKey),
  data= df_aves_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)
par(mfrow=c(2,4))
plot(fit.aves.surv)

fit.aves.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
    cluster(speciesKey),
  data= df_aves_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)
par(mfrow=c(2,4))
plot(fit.aves.completeness)

fit.aves <- recurrent.marginal.mean(
  fit.aves.surv,
  fit.aves.completeness
)

## Mammals ----
mammalia_event_100 <- phreg(Surv(t.start.year,
                                 t.stop.year,
                                 status==1)~
                              cluster(speciesKey)+
                              Verticality +
                              Nocturnality +
                              BodyLength_mm +
                              Latitude +
                              HumanDensity +
                              Elevation +
                              RangeSize+
                              strata(GlobalNorth),
                            data=df_mammalia_100)
summary(mammalia_event_100)
gof.mammalia <- gof(mammalia_event_100)

mammalia_terminal_100 <- phreg(Surv(t.start.year,
                                    t.stop.year,
                                    status==2)~                              
                                 cluster(speciesKey)+
                                 Verticality +
                                 Nocturnality +
                                 BodyLength_mm +
                                 Latitude +
                                 HumanDensity +
                                 Elevation +
                                 RangeSize+
                                 strata(GlobalNorth),
                               data=df_mammalia_100)
summary(mammalia_terminal_100)
gof.terminal.mammalia <- gof(mammalia_terminal_100)

# Aalen model
fit.mammalia.surv <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 1) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
    cluster(speciesKey),
  data= df_mammalia_100,
  start.time = 1760,
  max.time     = 2025,
  #residuals = 1,
  robust       = 1,       # variância robusta
  n.sim        = 100,     # número de simulações para CI
  resample.iid = 1       # para funções de residuais/CI i.i.d.
)

fit.mammalia.completeness <- aalen(
  formula = Surv(t.start.year, t.stop.year, status == 2) ~
    Verticality +
    Nocturnality +
    BodyLength_mm +
    HumanDensity +
    Latitude +
    Elevation +
    RangeSize +
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
    "aalen_completeness_models.RData"))

# Results plots ----
## Cumulative hazard plots ----
tetrapods_100 <- recurrentMarginal(Event(t.start.year, t.stop.year, status==1)~
                           strata(GlobalNorth)+
                           cluster(speciesKey),
                         data=df_wallacean_100,
                         cause=1,
                         death.code=2)

plot(tetrapods_100,
     se = TRUE,
     legend = FALSE,
     ylim = c(0, 30),
     xlim = c(1750, 2025),
     xlab = "Year",
     ylab = "Marginal mean",
     cex.main = 1.5,
     col = c("#c40e3e", "#2DA5E8"
)) 
# insert lines every 25 years
axis(1, at = seq(1750, 2025, by = 25), las = 1)  
abline(v = seq(1750, 2025, by = 25), 
       col = "gray70",  # cor cinza clara
       lty = 3,         # linha tracejada
       lwd = 0.5)       # linha bem fina

# TODO: Does the difference persist?
base <- mets:::basecumhaz(
  tetrapods_100,
  joint = 1,
  robust = FALSE,
  cumhaz = "cumhaz",
  se.cumhaz = "se.cumhaz"
)

H0 <- as.data.frame(base[[1]]$cumhaz)
H1 <- as.data.frame(base[[2]]$cumhaz)

H0 <- H0[, c("time", "cumhaz")]
H1 <- H1[, c("time", "cumhaz")]

H <- full_join(H0, H1, by = "time",
               suffix = c("_0", "_1")) |>
  arrange(time) |>
  fill(cumhaz_0, cumhaz_1, .direction = "down") 

# absolute difference
H_abs <- H |>
  mutate(delta_cumhaz = cumhaz_1 - cumhaz_0) |>
  filter(!is.na(delta_cumhaz))

plot(H_abs$time, H_abs$delta_cumhaz,
     type = "s",
     lwd = 2,
     ylim = c(0, 0.2),
     xlim = c(1875, 2025),
     col = "purple",
     xlab = "Ano",
     ylab = expression(Delta~(t)))
abline(h = 0, lty = 2)
axis(1, at = seq(1850, 2025, by = 25), las = 1)  
abline(v = seq(1850, 2025, by = 25), 
       col = "gray70",  # cor cinza clara
       lty = 3,         # linha tracejada
       lwd = 0.5)       # linha bem fina

# ratio difference
H_raz <- H |> mutate(ratio_cumhaz = cumhaz_1 / cumhaz_0)

plot(H_raz$time, H_raz$ratio_cumhaz, 
     type = "l",
     col = "darkgreen",
     ylim = c(0, 15),
     xlim = c(1875, 2025),
     lwd = 2,
     xlab = "Ano",
     ylab = "North/South Ratio")
abline(h = 1, lty = 2)
axis(1, at = seq(1850, 2025, by = 25), las = 1)  
abline(v = seq(1850, 2025, by = 25), 
       col = "gray70",  # cor cinza clara
       lty = 3,         # linha tracejada
       lwd = 0.5)       # linha bem fina

# Each group
cores <- c("Amphibia" = "#75C6FF",
           "Reptilia" = "#86A94D",
           "Aves" = "#7F5C3B",
           "Mammalia" = "#D66FFF")

par(mfrow=c(1,3),
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
  "BodyLength_mm",
  "Nocturnality", 
  "Verticality",
  "RangeSize",
  "Elevation",
  "HumanDensity",
  "Latitude"
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
    axis.text.x = element_text(size = 15),
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

ggsave("Figures/ModelPlot.pdf",
       plot=coefplot, width=15, height=8, units="in",
       dpi = 'print', cairo_pdf)

## GOF plot ----
dev.off()
par(mfrow=c(4,7),
    mar=c(2.5,2.5,2,1),
    mgp=c(1.7,0.7,0),
    cex.axis=1.0,     # números menores
    cex.lab=1.2)      # labels menores
par(las=2)  # coloca labels na vertical

plot_gof_custom(gof.terminal.amphibia, titles = TRUE, bg.col = "#75C6FF")
plot_gof_custom(gof.terminal.reptilia, titles = FALSE, bg.col = "#86A94D")
plot_gof_custom(gof.terminal.aves, titles = FALSE, bg.col = "#7F5C3B")
plot_gof_custom(gof.terminal.mammalia, titles = FALSE, bg.col = "#D66FFF")

plot_gof_custom(gof.amphibia, titles = TRUE, bg.col = "#75C6FF",
                show.xaxis = TRUE, show.yaxis = TRUE,
                show.xticks = TRUE, show.yticks = TRUE)
plot_gof_custom(gof.reptilia, titles = TRUE, bg.col = "#86A94D",
                show.xaxis = TRUE, show.yaxis = TRUE,
                show.xticks = TRUE, show.yticks = TRUE)
plot_gof_custom(gof.aves, titles = TRUE, bg.col = "#7F5C3B",
                show.xaxis = TRUE, show.yaxis = TRUE,
                show.xticks = TRUE, show.yticks = TRUE)
plot_gof_custom(gof.mammalia, titles = TRUE, bg.col = "#D66FFF",
                show.xaxis = TRUE, show.yaxis = TRUE,
                show.xticks = TRUE, show.yticks = TRUE)

## Time varying effects ----
dev.off()

par(mfrow=c(2,4))
plot(fit.amphibia.completeness)
for(i in 1:8) {
  par(mfg = c(ceiling(i/4), (i-1)%%4 + 1))  # posição na matriz 2x4
  
  # Adicionar as linhas
  abline(v = seq(1750, 2025, by = 25), 
         col = "gray70", 
         lty = 3, 
         lwd = 0.5)
}
