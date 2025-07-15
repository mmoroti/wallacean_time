# Packages and functions ----
library(tidyverse)
library(timereg)
library(reda)
library(reReg)
library(recforest)
library(survival)
library(coxme)
library(xgboost)
library(FSA)         # para o teste de Dunn
library(ggpubr)      # para anotar os p-valores (se quiser)
library(rstatix)     # para ajudar na formatação dos resultados

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
           'P-val' = V6,
           'lower2.5%' = V7,
           'upper97.5%' = V8) %>%
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

plot_completeness_boxplot <- function(df,
                                      eixo_y = NULL, 
                                      title_plot = "Completeness") {
  
  # Define cores fixas para as classes
  cores <- c(
    "Amphibia" = "#4FCFF5",
    "Reptilia" = "#A8896A",
    "Aves"     = "#F58C25",
    "Mammalia" = "#F54952"
  )
  
  # Força a ordem dos fatores
  df <- df %>%
    mutate(Class = factor(Class, levels = c("Reptilia", "Mammalia", "Amphibia", "Aves")))
  
  # Teste Kruskal-Wallis
  kruskal <- kruskal.test(tempo_total ~ Class, data = df)
  
  # Se significativo, faz post hoc de Dunn
  if (kruskal$p.value < 0.05) {
    dunn <- dunnTest(tempo_total ~ Class, data = df, method = "bh")
    print("Teste de Dunn (p-valor ajustado):")
    print(dunn$res)
  } else {
    print("Kruskal-Wallis não significativo, p-valor:")
    print(kruskal$p.value)
  }
  
  # Plot
  plot <- ggplot(df, aes(x = Class, y = tempo_total, fill = Class, color = Class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    scale_fill_manual(values = cores) +
    scale_color_manual(values = cores) +
    labs(
      x = "",
      y = eixo_y,
      title = title_plot
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  return(plot)
}

# Load data ----
load(file = file.path(
  "02_data_analysis",
  "data_wallacean_time_edit.RData")
)

## Amphibia ----
df_amphibia_100 <- df_wallacean_100 %>%
  filter(Class == "Amphibia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_amphibia_75 <- df_wallacean_75 %>%
  filter(Class == "Amphibia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_amphibia_50 <- df_wallacean_50 %>%
  filter(Class == "Amphibia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_amphibia_100$WallaceCompletude) # 46 100% 
table(df_amphibia_75$Wallace75) # 135 75%
table(df_amphibia_50$Wallace50) # 474 50%
length(unique(df_amphibia$speciesKey)) # 753 spp

## Reptilia ----
df_reptilia_100 <- df_wallacean_100 %>%
  filter(Class == "Reptilia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_reptilia_75 <- df_wallacean_75 %>%
  filter(Class == "Reptilia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_reptilia_50 <- df_wallacean_50 %>%
  filter(Class == "Reptilia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_reptilia_100$event) # 3741 eventos
table(df_reptilia_100$WallaceCompletude) # 14 100% 

table(df_reptilia_75$event) # 3725 eventos
table(df_reptilia_75$Wallace75) # 59 chegam 75% 
table(df_reptilia_50$Wallace50) # 315 chegam a 50%

length(unique(df_reptilia_75$speciesKey)) # 754 spp

## Aves ----
df_aves_100 <- df_wallacean_100 %>%
  filter(Class == "Aves") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_aves_75 <- df_wallacean_75 %>%
  filter(Class == "Aves") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_aves_50 <- df_wallacean_50 %>%
  filter(Class == "Aves") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_aves_100$event) # 14560 eventos
table(df_aves_100$WallaceCompletude) # 202 100%
length(unique(df_aves_100$speciesKey)) # 2185 spp

table(df_aves_75$event) # 14345 eventos
table(df_aves_75$Wallace75) # 521 75 %
length(unique(df_aves_75$speciesKey)) # 2185 spp

table(df_aves_50$event) # 12250 eventos
table(df_aves_50$Wallace50) # 1510 50%
length(unique(df_aves_50$speciesKey)) # 2185 spp

## Mammals ----
df_mammals_100 <- df_wallacean_100 %>%
  filter(Class == "Mammalia") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_mammals_75 <- df_wallacean_75 %>%
  filter(Class == "Mammalia") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

df_mammals_50 <- df_wallacean_50 %>%
  filter(Class == "Mammalia") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_mammals_100$event) # 3903 eventos
table(df_mammals_100$WallaceCompletude) # 27 100%
length(unique(df_mammals_100$speciesKey)) # 701 spp

table(df_mammals_75$event) # 3875 eventos
table(df_mammals_75$Wallace75) # 48 75 %
length(unique(df_mammals_75$speciesKey)) # 701 spp

table(df_mammals_50$event) # 3553 eventos
table(df_mammals_50$Wallace50) # 311 50%
length(unique(df_mammals_50$speciesKey)) # 2185 spp

# Aalen model (timereg package) ----
# Recurrent event models with terminal event
## Amphibia ----
#start_time <- Sys.time()
# 100 %
fit.amphibia.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_amphibia_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)

amp.surv.100 <- summary_aalen(fit.amphibia.surv,
              Class = "Amphibia",
              Type = "surv.100")

fit.amphibia.death <- aalen(
  formula      = Surv(t.start, t.stop, WallaceCompletude) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_amphibia_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
amp.death.100 <- summary_aalen(fit.amphibia.death,
                              Class = "Amphibia",
                              Type = "death.100")

bind_rows(amp.surv.100,
          amp.death.100)

fit.amphibia <- recurrent.marginal.mean(fit.amphibia.surv, 
                                        fit.amphibia.death)

# 75%
fit.amphibia.surv.75 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_amphibia_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)

fit.amphibia.death.75 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace75) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_amphibia_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)

amp.surv.75 <- summary_aalen(fit.amphibia.surv.75,
                             Class = "Amphibia",
                             Type = "surv.75")

amp.death.75 <- summary_aalen(fit.amphibia.death.75,
                              Class = "Amphibia",
                              Type = "death.75")

fit.amphibia.75 <- recurrent.marginal.mean(fit.amphibia.surv.75, 
                                        fit.amphibia.death.75)

# 50%
fit.amphibia.surv.50 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_amphibia_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.amphibia.surv.50)

fit.amphibia.death.50 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace50) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_amphibia_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.amphibia.death.50)

amp.surv.50 <- summary_aalen(fit.amphibia.surv.50,
                             Class = "Amphibia",
                             Type = "surv.50")

amp.death.50 <- summary_aalen(fit.amphibia.death.50,
                              Class = "Amphibia",
                              Type = "death.50")

fit.amphibia.50 <- recurrent.marginal.mean(fit.amphibia.surv.50, 
                                           fit.amphibia.death.50)

amphibia_model <- bind_rows(
  amp.surv.100,
  amp.surv.75,
  amp.surv.50,
  amp.death.100,
  amp.death.75,
  amp.death.50)

## Reptilia ----
# 100%
fit.reptilia.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_reptilia_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.surv)

fit.reptilia.death <- aalen(
  formula      = Surv(t.start, t.stop, WallaceCompletude) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_reptilia_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.death)

rep.surv.100 <- summary_aalen(fit.reptilia.surv,
                             Class = "Reptilia",
                             Type = "surv.100")

rep.death.100 <- summary_aalen(fit.reptilia.death,
                              Class = "Reptilia",
                              Type = "death.100")

fit.reptilia <- recurrent.marginal.mean(fit.reptilia.surv,
                                           fit.reptilia.death)

# 75%
fit.reptilia.surv.75 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_reptilia_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.surv.75)

fit.reptilia.death.75 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace75) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_reptilia_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.death.75)

rep.surv.75 <- summary_aalen(fit.reptilia.surv.75,
                              Class = "Reptilia",
                              Type = "surv.75")

rep.death.75 <- summary_aalen(fit.reptilia.death.75,
                               Class = "Reptilia",
                               Type = "death.75")

fit.reptilia.75 <- recurrent.marginal.mean(fit.reptilia.surv.75,
                                           fit.reptilia.death.75)

# 50%
fit.reptilia.surv.50 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_reptilia_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.surv.50)

fit.reptilia.death.50 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace50) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_reptilia_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.death.50)

fit.reptilia.50 <- recurrent.marginal.mean(fit.reptilia.surv.50,
                                           fit.reptilia.death.50)

rep.surv.50 <- summary_aalen(fit.reptilia.surv.50,
                             Class = "Reptilia",
                             Type = "surv.50")

rep.death.50 <- summary_aalen(fit.reptilia.death.50,
                              Class = "Reptilia",
                              Type = "death.50")

reptilia_model <- bind_rows(
  rep.surv.100,
  rep.surv.75,
  rep.surv.50,
  rep.death.100,
  rep.death.75,
  rep.death.50)

## Aves ----
fit.aves.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_aves_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.surv)
fit.aves.death <- aalen(
  formula      = Surv(t.start, t.stop, WallaceCompletude) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_aves_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.death)

aves.surv.100 <- summary_aalen(fit.aves.surv,
                             Class = "Aves",
                             Type = "surv.100")

aves.death.100 <- summary_aalen(fit.aves.death,
                              Class = "Aves",
                              Type = "death.100")

fit.aves <- recurrent.marginal.mean(fit.aves.surv, fit.aves.death)

# 75%
fit.aves.surv.75 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_aves_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.surv.75)

fit.aves.death.75 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace75) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_aves_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.death.75)

aves.surv.75 <- summary_aalen(fit.aves.surv.75,
                               Class = "Aves",
                               Type = "surv.75")

aves.death.75 <- summary_aalen(fit.aves.death.75,
                                Class = "Aves",
                                Type = "death.75")

fit.aves.75 <- recurrent.marginal.mean(fit.aves.surv.75, fit.aves.death.75)

# 50%
fit.aves.surv.50 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_aves_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.surv.50)

fit.aves.death.50 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace50) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_aves_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.death.50)

aves.surv.50 <- summary_aalen(fit.aves.surv.50,
                              Class = "Aves",
                              Type = "surv.50")

aves.death.50 <- summary_aalen(fit.aves.death.50,
                               Class = "Aves",
                               Type = "death.50")

aves_model <- bind_rows(
  aves.surv.100,
  aves.surv.75,
  aves.surv.50,
  aves.death.100,
  aves.death.75,
  aves.death.50)

fit.aves.50 <- recurrent.marginal.mean(fit.aves.surv.50, fit.aves.death.50)

## Mammals ----
fit.mammalia.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
    data         = df_mammals_100,
    robust       = 1,       # variância robusta
    n.sim        = 1000,     # número de simulações para CI
    resample.iid = 1,       # para funções de residuais/CI i.i.d.
    start.time   = 0,       # opcional, default = 0
    max.time     = NULL     # default = máximo observado
)
summary(fit.mammalia.surv)

fit.mammalia.death <- aalen(
  formula      = Surv(t.start, t.stop, WallaceCompletude) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_mammals_100,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.mammalia.death)

mammalia.surv.100 <- summary_aalen(fit.mammalia.surv,
                              Class = "Mammalia",
                              Type = "surv.100")

mammalia.death.100 <- summary_aalen(fit.mammalia.death,
                              Class = "Mammalia",
                              Type = "death.100")

fit.mammalia <- recurrent.marginal.mean(fit.mammalia.surv,
                                        fit.mammalia.death)

# 75%
fit.mammalia.surv.75 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_mammals_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.mammalia.surv.75)

fit.mammalia.death.75 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace75) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_mammals_75,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.mammalia.death.75)

mammalia.surv.75 <- summary_aalen(fit.mammalia.surv.75,
                                   Class = "Mammalia",
                                   Type = "surv.75")

mammalia.death.75 <- summary_aalen(fit.mammalia.death.75,
                                    Class = "Mammalia",
                                    Type = "death.75")

fit.mammalia.75 <- recurrent.marginal.mean(fit.mammalia.surv.75,
                                        fit.mammalia.death.75)

# 50%
fit.mammalia.surv.50 <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_mammals_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.mammalia.surv.50)

fit.mammalia.death.50 <- aalen(
  formula      = Surv(t.start, t.stop, Wallace50) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    #const(RangeSize) +
    cluster(speciesKey),
  data         = df_mammals_50,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.mammalia.death.50)

mammalia.surv.50 <- summary_aalen(fit.mammalia.surv.50,
                                  Class = "Mammalia",
                                  Type = "surv.50")

mammalia.death.50 <- summary_aalen(fit.mammalia.death.50,
                                   Class = "Mammalia",
                                   Type = "death.50")

mammals_model <- bind_rows(
  mammalia.surv.100,
  mammalia.surv.75,
  mammalia.surv.50,
  mammalia.death.100,
  mammalia.death.75,
  mammalia.death.50)

fit.mammalia.50 <- recurrent.marginal.mean(fit.mammalia.surv.50,
                                           fit.mammalia.death.50)
end_time <- Sys.time()
print(end_time - start_time)
# Results plots ----
par(mfrow = c(2, 2))
with(fit.amphibia, plot(times, mu, type = "s", main = "Amphibia"))
with(fit.amphibia, lines(times, mu + 1.96 * se.mu, type = "s", lty = 2))
with(fit.amphibia, lines(times, mu - 1.96 * se.mu, type = "s", lty = 2))

with(fit.reptilia, plot(times, mu, type = "s", main = "Reptilia"))
with(fit.reptilia, lines(times, mu + 1.96 * se.mu, type = "s", lty = 2))
with(fit.reptilia, lines(times, mu - 1.96 * se.mu, type = "s", lty = 2))

with(fit.aves, plot(times, mu, type = "s", main = "Aves"))
with(fit.aves, lines(times, mu + 1.96 * se.mu, type = "s", lty = 2))
with(fit.aves, lines(times, mu - 1.96 * se.mu, type = "s", lty = 2))

with(fit.mammalia, plot(times, mu, type = "s", main = "Mammals"))
with(fit.mammalia, lines(times, mu + 1.96 * se.mu, type = "s", lty = 2))
with(fit.mammalia, lines(times, mu - 1.96 * se.mu, type = "s", lty = 2))
dev.off()  

## Completeness scenarios ----
fits.50 <- list(
  Amphibia = fit.amphibia.50,
  Reptilia = fit.reptilia.50,
  Aves     = fit.aves.50,
  Mammalia = fit.mammalia.50
)

fits.75 <- list(
  Amphibia = fit.amphibia.75,
  Reptilia = fit.reptilia.75,
  Aves     = fit.aves.75,
  Mammalia = fit.mammalia.75
)

fits <- list(
  Amphibia = fit.amphibia,
  Reptilia = fit.reptilia,
  Aves     = fit.aves,
  Mammalia = fit.mammalia
)

cores <- c("#4FCFF5",
           "#A8896A",
           "#F58C25",
           "#F54952")

par(mfrow = c(3, 1))
plot_cumulative_events(fits.50,
                       cores = cores,
                       legend = TRUE,
                       title =  "50% Completeness",
                       xlab = "",
                       ylab = "") 

plot_cumulative_events(fits.75,
                       cores = cores,
                       legend = FALSE,
                       title =  "75% Completeness",
                       xlab = "",
                       ylab = "Cumulative number of events μ(t)") 

plot_cumulative_events(fits,
                       cores = cores,
                       legend = FALSE,
                       title =  "100% Completeness",
                       xlab = "Time (Years)",
                       ylab = "") 
dev.off()

## Time-to-completeness ----
temp.100 <- df_wallacean_100 %>%
  group_by(speciesKey, scientificName, Class) %>%
  filter(WallaceCompletude == 1) %>%
  summarise(
    tempo_total = max(t.stop),
  )
temp.75 <- df_wallacean_75 %>%
  group_by(speciesKey, scientificName, Class) %>%
  filter(Wallace75 == 1) %>%
  summarise(
    tempo_total = max(t.stop),
  )
temp.50 <- df_wallacean_50 %>%
  group_by(speciesKey, scientificName, Class) %>%
  filter(Wallace50 == 1) %>%
  summarise(
    tempo_total = max(t.stop),
  )

# 50%
letras <- data.frame(
  Class = c("Reptilia", "Mammalia", "Amphibia", "Aves"),
  letra = c("a", "b", "c", "c")
)
plot.50 <- plot_completeness_boxplot(temp.50, title_plot = "50% completeness",
                                     eixo_y = "") +
  geom_text(data = letras, aes(x = Class, y = 125, label = letra), 
            inherit.aes = FALSE, size = 5, fontface = "bold") 

# 75%
letras <- data.frame(
  Class = c("Reptilia", "Mammalia", "Amphibia", "Aves"),
  letra = c("a", "ab", "ab", "b")
)
plot.75 <- plot_completeness_boxplot(temp.75,
                                     eixo_y = "Time-to-completeness (Years)",
                                     title_plot = "75% completeness") +
  geom_text(data = letras, aes(x = Class, y = 125, label = letra), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

# 100%
letras <- data.frame(
  Class = c("Reptilia", "Mammalia", "Amphibia", "Aves"),
  letra = c("a", "a", "a", "a")
)
plot.100 <- plot_completeness_boxplot(temp.100,
                                      eixo_y = "",
                                      title_plot = "100% completeness") +
  geom_text(data = letras, aes(x = Class, y = 125, label = letra), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

plot <- cowplot::plot_grid(plot.50, plot.75, plot.100,
                   ncol = 1, nrow = 3)

## Coefficient plot ----
aalen_models <- bind_rows(
  amphibia_model,
  reptilia_model,
  aves_model,
  mammals_model
)

head(aalen_models)

## Derivada ---- 
d_mu <- diff(fit.amphibia$mu) / diff(fit.amphibia$times)
t_deriv <- fit.amphibia$times[-1]

plot(t_deriv, log(d_mu), type = "l", ylab = "Taxa de descoberta (novos municípios por ano)",
     xlab = "Tempo (anos)", main = "Taxa de descoberta ao longo do tempo (Amphibia)")

# Suavização com loess
smoothed <- loess(log(d_mu) ~ t_deriv, span = 0.1)

plot(t_deriv, log(d_mu), type = "l", col = "gray", 
     ylab = "Taxa de descoberta", xlab = "Tempo", 
     main = "Taxa de descoberta suavizada (Amphibia)")

lines(t_deriv, predict(smoothed), col = "blue", lwd = 2)

# Cox-Aalen models ----
## Amphibia ----
amphibia.model.surv <- cox.aalen(Surv(t.start,t.stop,event) ~ 1+ 
                                   prop(Verticality) +
                                   prop(Nocturnality) +
                                   prop(log(BodyLength_mm)) +
                                   cluster(speciesKey),
                                 robust = 1,
                                 data= df_amphibia,
                                 resample.iid=1, max.clust=NULL,max.timepoint.sim=NULL)

summary(amphibia.model.surv)

amphibia.model.death <- cox.aalen(Surv(t.start,t.stop,event) ~ 1 +
                                    prop(Verticality^2) +
                                    prop(Nocturnality^2) +
                                    prop(log(BodyLength_mm)) +
                                    cluster(speciesKey),
                                  data= df_amphibia,
                                  robust = 1,
                                  resample.iid=1, max.clust=NULL, max.timepoint.sim=NULL)

amphibia.model <- recurrent.marginal.coxmean(
  amphibia.model.surv,
  amphibia.model.death)

summary(amphibia.model)
with(amphibia.model,plot(times,mu,type="s"))
with(amphibia.model,lines(times,mu+1.96*se.mu,type="s",lty=2))

## Reptilia ----
reptilia.model.surv <- cox.aalen(Surv(t.start,t.stop,event) ~ 1+ 
                                   prop(Verticality) +
                                   prop(Nocturnality) +
                                   prop(log(BodyLength_mm)) +
                                   cluster(speciesKey),
                                 robust = 1,
                                 data= df_reptilia,
                                 resample.iid=1, max.clust=NULL,max.timepoint.sim=NULL)

summary(reptilia.model.surv)

reptilia.model.death <- cox.aalen(Surv(t.start,t.stop,event) ~ 1 +
                                    prop(Verticality^2) +
                                    prop(Nocturnality^2) +
                                    prop(log(BodyLength_mm)) +
                                    cluster(speciesKey),
                                  data= df_reptilia,
                                  robust = 1,
                                  resample.iid=1, max.clust=NULL, max.timepoint.sim=NULL)

reptilia.model <- recurrent.marginal.coxmean(
  reptilia.model.surv,
  reptilia.model.death)

summary(reptilia.model)
with(reptilia.model,plot(times,mu,type="s"))
with(reptilia.model,lines(times,mu+1.96*se.mu,type="s",lty=2))

## Aves ----
aves.model.surv <- cox.aalen(Surv(t.start,t.stop,event) ~ 1+ 
                               prop(Verticality) +
                               prop(Nocturnality) +
                               prop(log(BodyLength_mm)) +
                               cluster(speciesKey),
                             robust = 1,
                             data= df_aves,
                             resample.iid=1, max.clust=NULL,max.timepoint.sim=NULL)

summary(aves.model.surv)
plot(aves.model.surv, cum = TRUE, se = TRUE)
aves.model.death <- cox.aalen(Surv(t.start,t.stop,event) ~ 1 +
                                prop(Verticality^2) +
                                prop(Nocturnality^2) +
                                prop(log(BodyLength_mm)) +
                                cluster(speciesKey),
                              data= df_aves,
                              robust = 1,
                              resample.iid=1, max.clust=NULL, max.timepoint.sim=NULL)
summary(aves.model.death)

aves.model <- recurrent.marginal.coxmean(
  aves.model.surv,
  aves.model.death)

summary(aves.model)
with(aves.model,plot(times,mu,type="s"))
with(aves.model,lines(times,mu+1.96*se.mu,type="s",lty=2))

## Mammalia ----
mammals.model.surv <- cox.aalen(Surv(t.start,t.stop,event) ~ 1+ 
                                  prop(Verticality) +
                                  prop(Nocturnality) +
                                  prop(log(BodyLength_mm)) +
                                  cluster(speciesKey),
                                robust = 1,
                                data= df_mammals,
                                resample.iid=1, max.clust=NULL,max.timepoint.sim=NULL)
summary(mammals.model.surv)

mammals.model.death <- cox.aalen(Surv(t.start,t.stop,event) ~ 1 +
                                   prop(Verticality^2) +
                                   prop(Nocturnality^2) +
                                   prop(log(BodyLength_mm)) +
                                   cluster(speciesKey),
                                 data= df_mammals,
                                 robust = 1,
                                 resample.iid=1, max.clust=NULL, max.timepoint.sim=NULL)
summary(mammals.model.death)

mammals.model <- recurrent.marginal.coxmean(
  mammals.model.surv,
  mammals.model.death)

summary(mammals.model)
with(mammals.model,plot(times,mu,type="s"))
with(mammals.model,lines(times,mu+1.96*se.mu,type="s",lty=2))

# reReg package ----
# PRECISA ENTENDER PQ EM ALGUNS CASOS A DATA EH IGUAL
df_wallacean_time %>%
  filter(t.start >= t.stop) %>%
  View()

df_teste <- df_wallacean_time %>%
  #filter(speciesKey %in% c("2422148","2428622")) %>%
  filter(speciesKey %in% sample(c(df_wallacean_time$speciesKey), 5)) #%>%
  #mutate(t.start = if_else(t.start > 0, t.start + 0.00000001, 0))

View(df_wallacean_time)
 
try(
  with(df_teste,
    reda::Recur(t.start %to% t.stop,
          id = speciesKey,
          event = event,
          terminal = WallaceCompletude))
)

fm <- with(df_teste,
           reda::Recur(t.start %to% t.stop,
                       id = speciesKey,
                       event = event,
                       terminal = WallaceCompletude)~1) 

reReg(fm, data = df_teste, model = "cox", B = 50)

unique(df_wallacean_time$speciesKey)[86]
unique(df_teste$speciesKey)[11]
View(df_wallacean_time)

# recforest package ----
trained_forest <- train_forest(
  data = df_amphibia,
  id_var = "speciesKey",
  covariates = c("Nocturnality", "Verticality", "BodyLength_mm"),
  time_vars = c("t.start", "t.stop"),
  death_var = "WallaceCompletude",
  event = "event",
  n_trees = 3,
  n_bootstrap = 70,
  mtry = 2,
  minsplit = 3,
  nodesize = 15,
  method = "GL",
  min_score = 5,
  max_nodes = 20,
  seed = 111,
  parallel = TRUE,
  verbose = FALSE
)

summary(train_forest)

predictions <- predict(
  trained_forest,
  newdata = df_eventos,
  id_var = "speciesKey",
  covariates = c("Diu", "Noc", "BodyLength_mm"),
  time_vars = c("t.start", "t.stop"),
  death_var = "WallaceCompletude"
) 

# Lista das variáveis binárias em probabilidade
variaveis_binarias <- c("Diu", "Noc", "Fos", "Ter", "Aqu", "Arb")

# Aplicar a transformação para cada variável
for (variavel in variaveis_binarias) {
  data_wallacean_time[[variavel]] <- ifelse(data_wallacean_time[[variavel]] > 0.5, 1, 0)
}

table(data_wallacean_time$Diu)
table(data_wallacean_time$Noc)


View(data_wallacean_time[,c("scientificName", "event", "year_modified", "sigla_admu")])
# Ajustando o modelo de Cox
modelo_cox <- coxme(Surv(time_to_event, event) ~ BodyLength_mm + Diu + Noc +
 Fos + Ter + Aqu + Arb + (1 | speciesKey), data = data_wallacean_time)

# Resumo do modelo
summary(modelo_cox)

# Calcular o C-index
concordancia <- concordance(Surv(data_wallacean_time$time_to_event, data_wallacean_time$event) ~ predicoes)
print(concordancia$concordance)

### Modelos
data_xgb <- data_wallacean_time %>%
  select(speciesKey, BodyLength_mm, Diu, Noc, year_modified, time_to_event, event) %>%
  mutate(time_to_event = time_to_event + 1,
         BodyLength_mm = BodyLength_mm)
glimpse(data_xgb)

# Converter variáveis categóricas (se houver)
data_xgb$speciesKey <- as.factor(data_xgb$speciesKey)

# Separar covariáveis e variáveis de resposta
X <- data_xgb[, c("BodyLength_mm", "Diu", "Noc")]
y <- data_xgb[, c("time_to_event", "event")]

dtrain <- xgb.DMatrix(data = as.matrix(X), label = y$time_to_event)

# Definir parâmetros do modelo
params <- list(
  objective = "survival:cox",  # Objetivo de sobrevivência
  eval_metric = "cox-nloglik",  # Métrica de avaliação
  max_depth = 6,                # Profundidade máxima da árvore
  eta = 0.1,                    # Taxa de aprendizado
  subsample = 0.8,              # Subamostragem
  colsample_bytree = 0.8        # Subamostragem de colunas
)

# Treinar o modelo
model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,  # Número de iterações
  watchlist = list(train = dtrain)
)

# Fazer previsões
predictions <- predict(model, as.matrix(X))

# Adicionar previsões ao data frame
data_xgb$predicted_risk <- predictions

# Gerar curvas de sobrevivência para cada espécie
surv_curves <- survfit(Surv(time_to_event, event) ~ speciesKey, data = data_xgb)
head(surv_curves)

# Plotar curvas de sobrevivência
survminer::ggsurvplot(
  surv_curves,
  fun = "event", 
  data = data_xgb,
  pval = TRUE,
  conf.int = TRUE,
  palette = "hue",
  xlab = "Tempo (anos)",
  ylab = "Probabilidade de Sobrevivência",
  legend = "none"
  #legend.title = "Espécie",
  #legend.labs = levels(data_xgb$speciesKey)
)
survminer::ggsurvplot(
  surv_curves,
  fun = "event", 
  data = data_xgb,
  pval = TRUE,
  conf.int = FALSE,  # Remove a área de confiança (sombrada)
  palette = "hue",
  xlab = "Tempo (anos)",
  ylab = "Probabilidade de Sobrevivência",
  legend = "none",  # Remove a legenda
  linetype = "solid"  # Garantir que as linhas não estão preenchidas
)

selected_species <- sample(unique(data_xgb$speciesKey), 10)  # Seleciona 10 espécies aleatórias
data_xgb_subset <- data_xgb[data_xgb$speciesKey %in% selected_species, ]

survminer::ggsurvplot(
  surv_curves,
  fun = "event", 
  data = data_xgb_subset,
  pval = TRUE,
  conf.int = FALSE, 
  palette = "hue",
  xlab = "Tempo (anos)",
  ylab = "Probabilidade de Sobrevivência",
  legend = "none",
  linetype = "solid"
)

# IA MODEL USING RANDOM FOREST ----
library(randomForestSRC)

model <- rfsrc(Surv(time_to_event, event) ~ BodyLength_mm + Diu + Noc,
  ntree = 1000,
  nodesize = 15,
  mtry = 2,
  strata = data_xgb$speciesKey,
  data = data_xgb)

summary(model)
?rfsrc
importance <- vimp(model)
print(importance)

# IA MODEL USING XGBOOST
# Lista única de espécies
species <- unique(data_xgb$speciesKey)

# Separar as espécies em treino e teste (80/20)
train_species <- sample(species, size = round(0.8 * length(species)))
test_species <- setdiff(species, train_species)

# Filtrar os dados com base nas espécies selecionadas
train_data <- data_xgb[data_xgb$speciesKey %in% train_species, ]
test_data <- data_xgb[data_xgb$speciesKey %in% test_species, ]

# Verificar a divisão
table(train_data$speciesKey %in% test_data$speciesKey)  # Deve ser sempre 

# Garantir que 'event' (censura) está no formato numérico correto
train_data$event <- as.numeric(train_data$event)
test_data$event <- as.numeric(test_data$event)

# Preparar as variáveis
train_matrix <- model.matrix(~ BodyLength_mm + Diu + Noc + time_to_event - 1, data = train_data)
test_matrix <- model.matrix(~ BodyLength_mm + Diu + Noc + time_to_event - 1, data = test_data)

train_label <- train_data$time_to_event  # Tempo até o evento
train_censorship <- train_data$event    # Censura (1 se o evento ocorreu, 0 se censurado)

test_label <- test_data$time_to_event
test_censorship <- test_data$event

# Parâmetros do modelo
params <- list(
  objective = "survival:aft",  # Modelo AFT
  eval_metric = "aft-nloglik",  # Log-verossimilhança de AFT
  aft_loss_distribution = "normal",  # Distribuição normal para o AFT
  aft_loss_distribution_scale = 1.0,  # Escala da distribuição
  eta = 0.1,  # Taxa de aprendizado
  max_depth = 6,  # Profundidade máxima das árvores
  subsample = 0.8,  # Subamostragem
  colsample_bytree = 0.8  # Subamostragem de colunas
)

# Conversão para DMatrix
dtrain <- xgb.DMatrix(data = train_matrix, label = train_label, missing = NA, base_margin = train_censorship)
dtest <- xgb.DMatrix(data = test_matrix, label = test_label, missing = NA, base_margin = test_censorship)

# Treinamento do modelo
set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 1000,
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 10,
  print_every_n = 10
)

# Verifique se há valores ausentes
summary(train_data$time_to_event)
summary(train_data$event)

# Verifique se as variáveis estão no formato correto
str(train_data$time_to_event)
str(train_data$event)

nrow(train_matrix)
length(train_data$time_to_event)
