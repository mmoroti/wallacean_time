library(tidyverse)
library(timereg)
library(reda)
library(reReg)
library(recforest)
library(survival)
library(coxme)
library(xgboost)

load(file = file.path(
      "02_data_analysis",
      "data_wallacean_time.RData")
)

# separando por classe
df_amphibia <- df_wallacean_time %>%
  filter(Class == "Amphibia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_amphibia$event) # 3055 eventos
table(df_amphibia$WallaceCompletude) # 43 100% 
length(unique(df_amphibia$speciesKey)) # 753 spp

df_reptilia <- df_wallacean_time %>%
  filter(Class == "Reptilia") %>%
  mutate(BodyLength_mm = scale(log(BodyLength_mm)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_reptilia$event) # 3741 eventos
table(df_reptilia$WallaceCompletude) # 13 100% 
length(unique(df_reptilia$speciesKey)) # 754 spp

df_aves <- df_wallacean_time %>%
  filter(Class == "Aves") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_aves$event) # 14560 eventos
table(df_aves$WallaceCompletude) # 200 100%
length(unique(df_aves$speciesKey)) # 2185 spp

df_mammals <- df_wallacean_time %>%
  filter(Class == "Mammalia") %>%
  mutate(BodyMass_g = scale(log(BodyMass_g)),
         Verticality = scale(Verticality^2),
         Nocturnality = scale(Nocturnality^2),
         RangeSize = scale(log(RangeSize)))

table(df_mammals$event) # 3903 eventos
table(df_mammals$WallaceCompletude) # 27 spp 100%
length(unique(df_mammals$speciesKey)) # 701 spp

# timereg package ----
# Aalen models ----
## Amphibia ----
fit.amphibia.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_amphibia,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.amphibia.surv)

fit.amphibia.death <- aalen(
  formula      = Surv(t.start, t.stop, WallaceCompletude) ~
    const(Verticality) +
    const(Nocturnality) +
    const(log(BodyLength_mm)) +
    const(RangeSize) +
    cluster(speciesKey),
  data         = df_amphibia,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.amphibia.death)
fit.amphibia <- recurrent.marginal.mean(fit.amphibia.surv, fit.amphibia.death)

## Reptilia ----
fit.reptilia.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyLength_mm) +
    const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_reptilia,
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
    const(RangeSize) +
    cluster(speciesKey),
  data         = df_reptilia,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.reptilia.death)
fit.reptilia <- recurrent.marginal.mean(fit.reptilia.surv, fit.reptilia.death)

## Aves ----
fit.aves.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    const(RangeSize) +
    cluster(speciesKey),
  #id           = "speciesKey",
  data         = df_aves,
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
    const(RangeSize) +
    cluster(speciesKey),
  data         = df_aves,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
summary(fit.aves.death)
fit.aves <- recurrent.marginal.mean(fit.aves.surv, fit.aves.death)

## Mammals ----
fit.mammalia.surv <- aalen(
  formula      = Surv(t.start, t.stop, event) ~
    const(Verticality) +
    const(Nocturnality) +
    const(BodyMass_g) +
    const(RangeSize) +
    cluster(speciesKey),
    data         = df_mammals,
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
    const(RangeSize) +
    cluster(speciesKey),
  data         = df_mammals,
  robust       = 1,       # variância robusta
  n.sim        = 1000,     # número de simulações para CI
  resample.iid = 1,       # para funções de residuais/CI i.i.d.
  start.time   = 0,       # opcional, default = 0
  max.time     = NULL     # default = máximo observado
)
fit.mammalia <- recurrent.marginal.mean(fit.mammalia.surv,
                                        fit.mammalia.death)

## Plot ----
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
