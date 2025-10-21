library(cowplot)
library(INLA)
library(data.table)
library(tidyverse)
library(spdep)
library(dlnm)
library(tsModel)
library(sf)
library(sp)
library(RColorBrewer)
library(geofacet)
library(ggpubr)
library(ggthemes)
library(splines)
library(raster)
library(grid)
library(ggplot2) 
library(ggsci)

#================== Data import =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  mutate(Pref = Pref + 16)%>%
  dplyr::select(year, month, X, case, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop, regnames, MT, RH) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(year, month, X, case, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop, regnames, MT, RH)%>%
  mutate(prov = "AH")

data <- rbind(AH.data, SC.data) %>%
  mutate(Pref.year = Pref,
         Pref.cos1 = Pref,
         Pref.cos2 = Pref,
         Pref.sin1 = Pref,
         Pref.sin2 = Pref,
         cos_i = cos(2*pi/12*X),
         sin_i = sin(2*pi/12*X),
         cos_i2 = cos(4*pi/12*X),
         sin_i2 = sin(4*pi/12*X),
         E = pop*10000)

#================== Baseline model selection =====================
fitmodel <- function(formula, data = data, family = "nbinomial")
{
  model <- inla(formula = formula, data = data, family = family, 
                control.compute = list(dic = TRUE),
                control.fixed = list(correlation.matrix = TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                verbose = FALSE,
                safe = TRUE)
  model <- inla.rerun(model)
  return(model)
}

all.options <- expand.grid(intercept = c("same", "iid", "spatial"), 
                           longterm.spat = c("same", "iid", "spatial"), 
                           longterm.temporal = c("ar1", "rw1"),
                           seasonal.spat = c("same", "iid", "spatial"),
                           seasonal.temp = c("1Fourier", "2Fourier"))

all.options$DIC = NA

for(i in 97:nrow(all.options))
{
  print(i)
  part1 <- case_when(
    all.options$intercept[i] == "same" ~ "1",
    all.options$intercept[i] == "iid" ~ "1 + f(Pref, model= \"iid\", constr = TRUE)",
    all.options$intercept[i] == "spatial" ~ "1 + f(Pref, model = \"bym\", graph = Both.adj)"
  )
  part2 <- case_when(
    all.options$longterm.spat[i] == "same" & all.options$longterm.temporal[i] == "ar1" ~ "f(year, model = \"ar1\", constr = TRUE)",
    all.options$longterm.spat[i] == "same" & all.options$longterm.temporal[i] == "rw1" ~ "f(year, model = \"rw1\", constr = TRUE)",
    all.options$longterm.spat[i] == "iid" & all.options$longterm.temporal[i] == "ar1" ~ "f(year, model = \"ar1\", constr = TRUE,  group = Pref.year, control.group = list(model = \"iid\"))",
    all.options$longterm.spat[i] == "iid" & all.options$longterm.temporal[i] == "rw1" ~ "f(year, model = \"rw1\", constr = TRUE, group = Pref.year, control.group = list(model = \"iid\"))",
    all.options$longterm.spat[i] == "spatial" & all.options$longterm.temporal[i] == "ar1" ~ "f(year, model = \"ar1\", constr = TRUE, group = Pref.year, control.group = list(model = \"besag\", graph = Both.adj))",
    all.options$longterm.spat[i] == "spatial" & all.options$longterm.temporal[i] == "rw1" ~ "f(year, model = \"rw1\", constr = TRUE, group = Pref.year, control.group = list(model = \"besag\", graph = Both.adj))"
  )
  part3 <- case_when(
    all.options$seasonal.spat[i] == "same" & all.options$seasonal.temp[i] == "1Fourier" ~ "sin_i + cos_i",
    all.options$seasonal.spat[i] == "same" & all.options$seasonal.temp[i] == "2Fourier" ~ "sin_i + cos_i + sin_i2 + cos_i2",
    all.options$seasonal.spat[i] == "iid" & all.options$seasonal.temp[i] == "1Fourier" ~ "f(Pref.cos1, cos_i, model = \"iid\") + f(Pref.sin1, sin_i, model = \"iid\")",
    all.options$seasonal.spat[i] == "iid" & all.options$seasonal.temp[i] == "2Fourier" ~ "f(Pref.cos1, cos_i, model = \"iid\") + f(Pref.sin1, sin_i, model = \"iid\") + f(Pref.cos2, cos_i2, model = \"iid\") + f(Pref.sin2, sin_i2, model = \"iid\")",
    all.options$seasonal.spat[i] == "spatial" & all.options$seasonal.temp[i] == "1Fourier" ~ "f(Pref.cos1, cos_i, model = \"besag\", graph = Both.adj) + f(Pref.sin1, sin_i, model = \"besag\", graph = Both.adj)",
    all.options$seasonal.spat[i] == "spatial" & all.options$seasonal.temp[i] == "2Fourier" ~ "f(Pref.cos1, cos_i, model = \"besag\", graph = Both.adj) + f(Pref.sin1, sin_i, model = \"besag\", graph = Both.adj) + f(Pref.cos2, cos_i2, model = \"besag\", graph = Both.adj) + f(Pref.sin2, sin_i2, model = \"besag\", graph = Both.adj)"
  )
  current.formula <- paste("case ~  offset(log(E)) + CNY +", part1, " + ", part2, " + ", part3, sep = "")
  
  for (k in 1:20) {
    message(k)
    result <- try({
      current.model.fit0 <- fitmodel(as.formula(current.formula), data)   
    })
    if (class(result) != "try-error") {
      print("Success!")
      break
    }
  }
  all.options$DIC[i] <- ifelse(class(result) != "try-error", current.model.fit0$dic$dic, NA)
  
  write.csv(all.options, "Model_DIC.csv", row.names = FALSE)
}

#================== Climate variable selection =====================
ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

ns.RH <- ns(data$RH, 3)
colnames(ns.RH) <- paste("RH.", colnames(ns.RH), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.RH

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  ns.RH

current.model <- fitmodel(model1, data)
dic_1 <- current.model$dic$dic

current.model <- fitmodel(model2, data)
dic_2 <- current.model$dic$dic

current.model <- fitmodel(model3, data)
dic_3 <- current.model$dic$dic

#======================== DLNM - Lag ==============================
slag = -3
nlag = 12

ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)

#basis
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(-3:nlag, 3)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_OM

#================== PM2.5 =====================
model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, 
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                        Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                        RR = as.numeric(result.PM2.5$matRRfit), 
                        RR.low = as.numeric(result.PM2.5$matRRlow), 
                        RR.high = as.numeric(result.PM2.5$matRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, 
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_mat <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                      RR = as.numeric(result.SO4$matRRfit), 
                      RR.low = as.numeric(result.SO4$matRRlow), 
                      RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, by= 5.666344,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                      RR = as.numeric(result.NO3$matRRfit), 
                      RR.low = as.numeric(result.NO3$matRRlow), 
                      RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, by= 3.564469,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                      RR = as.numeric(result.NH4$matRRfit), 
                      RR.low = as.numeric(result.NH4$matRRlow), 
                      RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, by= 1.037183,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                     Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                     RR = as.numeric(result.BC$matRRfit), 
                     RR.low = as.numeric(result.BC$matRRlow), 
                     RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, by= 5.179423,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                     Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                     RR = as.numeric(result.OM$matRRfit), 
                     RR.low = as.numeric(result.OM$matRRlow), 
                     RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat$Type <- "PM2.5"
SO4_mat$Type <- "Sulfate"
NO3_mat$Type <- "Nitrate"
NH4_mat$Type <- "Ammonium"
BC_mat$Type <- "BC"
OM_mat$Type <- "OM"

PM2.5_mat$var50 <- "40"
SO4_mat$var50 <- "8.5"
NO3_mat$var50 <- "9"
NH4_mat$var50 <- "7"
BC_mat$var50 <- "2"
OM_mat$var50 <- "10"

PM2.5_mat$var75 <- "60"
SO4_mat$var75 <- "12"
NO3_mat$var75 <- "16"
NH4_mat$var75 <- "10.5"
BC_mat$var75 <- "2.8"
OM_mat$var75 <- "14"

PM2.5_mat$var95 <- "95"
SO4_mat$var95 <- "17"
NO3_mat$var95 <- "25"
NH4_mat$var95 <- "16"
BC_mat$var95 <- "4.8"
OM_mat$var95 <- "24"

count_mat <- bind_rows(PM2.5_mat, SO4_mat, NO3_mat, NH4_mat, BC_mat, OM_mat)

#======================== DLNM - Exposure ==============================
slag = 0
nlag = 12

ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)

#basis
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 3)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, cyclic = TRUE, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_OM

percentiles <- c(0.05, 0.15,  0.25, 0.35,  0.45, 0.55,  0.65, 0.75, 0.85, 0.95)

specific_concentrations_PM2.5 <- quantile(data$PM2.5, percentiles, na.rm = TRUE)
specific_concentrations_SO4 <- quantile(data$SO4, percentiles, na.rm = TRUE)
specific_concentrations_NO3 <- quantile(data$NO3, percentiles, na.rm = TRUE)
specific_concentrations_NH4 <- quantile(data$NH4, percentiles, na.rm = TRUE)
specific_concentrations_BC <- quantile(data$BC, percentiles, na.rm = TRUE)
specific_concentrations_OM <- quantile(data$OM, percentiles, na.rm = TRUE)

#================== PM2.5 =====================
model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                        Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                        RR = as.numeric(result.PM2.5$cumRRfit), 
                        RR.low = as.numeric(result.PM2.5$cumRRlow), 
                        RR.high = as.numeric(result.PM2.5$cumRRhigh))


#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                      RR = as.numeric(result.SO4$cumRRfit), 
                      RR.low = as.numeric(result.SO4$cumRRlow), 
                      RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                      RR = as.numeric(result.NO3$cumRRfit), 
                      RR.low = as.numeric(result.NO3$cumRRlow), 
                      RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                      RR = as.numeric(result.NH4$cumRRfit), 
                      RR.low = as.numeric(result.NH4$cumRRlow), 
                      RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                     Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                     RR = as.numeric(result.BC$cumRRfit), 
                     RR.low = as.numeric(result.BC$cumRRlow), 
                     RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                     Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                     RR = as.numeric(result.OM$cumRRfit), 
                     RR.low = as.numeric(result.OM$cumRRlow), 
                     RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum$Type <- "PM2.5"
SO4_cum$Type <- "Sulfate"
NO3_cum$Type <- "Nitrate"
NH4_cum$Type <- "Ammonium"
BC_cum$Type <- "BC"
OM_cum$Type <- "OM"

exposure <- bind_rows(PM2.5_cum, SO4_cum, NO3_cum, NH4_cum, BC_cum, OM_cum)