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

#======================== Sensitivity Analysis 1 ==============================

#================== DLNM - Lag =====================
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

#================== lag knots 2 =====================
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
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

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$Type <- "PM2.5"
SO4_mat1$Type <- "Sulfate"
NO3_mat1$Type <- "Nitrate"
NH4_mat1$Type <- "Ammonium"
BC_mat1$Type <- "BC"
OM_mat1$Type <- "OM"

PM2.5_mat1$var50 <- "40"
SO4_mat1$var50 <- "8.5"
NO3_mat1$var50 <- "9"
NH4_mat1$var50 <- "7"
BC_mat1$var50 <- "2"
OM_mat1$var50 <- "10"

PM2.5_mat1$var75 <- "60"
SO4_mat1$var75 <- "12"
NO3_mat1$var75 <- "16"
NH4_mat1$var75 <- "10.5"
BC_mat1$var75 <- "2.8"
OM_mat1$var75 <- "14"

PM2.5_mat1$var95 <- "95"
SO4_mat1$var95 <- "17"
NO3_mat1$var95 <- "25"
NH4_mat1$var95 <- "16"
BC_mat1$var95 <- "4.8"
OM_mat1$var95 <- "24"

PM2.5_mat1$knot <- "2"
SO4_mat1$knot <- "2"
NO3_mat1$knot <- "2"
NH4_mat1$knot <- "2"
BC_mat1$knot <- "2"
OM_mat1$knot <- "2"

#================== lag knots 4 =====================
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
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

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat2$Type <- "PM2.5"
SO4_mat2$Type <- "Sulfate"
NO3_mat2$Type <- "Nitrate"
NH4_mat2$Type <- "Ammonium"
BC_mat2$Type <- "BC"
OM_mat2$Type <- "OM"

PM2.5_mat2$var50 <- "40"
SO4_mat2$var50 <- "8.5"
NO3_mat2$var50 <- "9"
NH4_mat2$var50 <- "7"
BC_mat2$var50 <- "2"
OM_mat2$var50 <- "10"

PM2.5_mat2$var75 <- "60"
SO4_mat2$var75 <- "12"
NO3_mat2$var75 <- "16"
NH4_mat2$var75 <- "10.5"
BC_mat2$var75 <- "2.8"
OM_mat2$var75 <- "14"

PM2.5_mat2$var95 <- "95"
SO4_mat2$var95 <- "17"
NO3_mat2$var95 <- "25"
NH4_mat2$var95 <- "16"
BC_mat2$var95 <- "4.8"
OM_mat2$var95 <- "24"

PM2.5_mat2$knot <- "4"
SO4_mat2$knot <- "4"
NO3_mat2$knot <- "4"
NH4_mat2$knot <- "4"
BC_mat2$knot <- "4"
OM_mat2$knot <- "4"

#================== lag knots 5 =====================
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat3 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
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

SO4_mat3<- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                      RR = as.numeric(result.SO4$matRRfit), 
                      RR.low = as.numeric(result.SO4$matRRlow), 
                      RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 


NO3_mat3 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat3 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat3 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat3 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat3$Type <- "PM2.5"
SO4_mat3$Type <- "Sulfate"
NO3_mat3$Type <- "Nitrate"
NH4_mat3$Type <- "Ammonium"
BC_mat3$Type <- "BC"
OM_mat3$Type <- "OM"

PM2.5_mat3$var50 <- "40"
SO4_mat3$var50 <- "8.5"
NO3_mat3$var50 <- "9"
NH4_mat3$var50 <- "7"
BC_mat3$var50 <- "2"
OM_mat3$var50 <- "10"

PM2.5_mat3$var75 <- "60"
SO4_mat3$var75 <- "12"
NO3_mat3$var75 <- "16"
NH4_mat3$var75 <- "10.5"
BC_mat3$var75 <- "2.8"
OM_mat3$var75 <- "14"

PM2.5_mat3$var95 <- "95"
SO4_mat3$var95 <- "17"
NO3_mat3$var95 <- "25"
NH4_mat3$var95 <- "16"
BC_mat3$var95 <- "4.8"
OM_mat3$var95 <- "24"

PM2.5_mat3$knot <- "5"
SO4_mat3$knot <- "5"
NO3_mat3$knot <- "5"
NH4_mat3$knot <- "5"
BC_mat3$knot <- "5"
OM_mat3$knot <- "5"

lagknots_contour <- bind_rows(PM2.5_mat1, PM2.5_mat2, PM2.5_mat3,
                              SO4_mat1, SO4_mat2, SO4_mat3,
                              NO3_mat1, NO3_mat2, NO3_mat3,
                              NH4_mat1, NH4_mat2, NH4_mat3,
                              BC_mat1, BC_mat2, BC_mat3,
                              OM_mat1, OM_mat2, OM_mat3)

#================== DLNM - Exposure =====================
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

#================== lag knots 2 =====================
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 2)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$Type <- "PM2.5"
SO4_cum1$Type <- "Sulfate"
NO3_cum1$Type <- "Nitrate"
NH4_cum1$Type <- "Ammonium"
BC_cum1$Type <- "BC"
OM_cum1$Type <- "OM"

PM2.5_cum1$knot <- "2"
SO4_cum1$knot <- "2"
NO3_cum1$knot <- "2"
NH4_cum1$knot <- "2"
BC_cum1$knot <- "2"
OM_cum1$knot <- "2"

#================== lag knots 4 =====================
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 4)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM ===================== 
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum2$Type <- "PM2.5"
SO4_cum2$Type <- "Sulfate"
NO3_cum2$Type <- "Nitrate"
NH4_cum2$Type <- "Ammonium"
BC_cum2$Type <- "BC"
OM_cum2$Type <- "OM"

PM2.5_cum2$knot <- "4"
SO4_cum2$knot <- "4"
NO3_cum2$knot <- "4"
NH4_cum2$knot <- "4"
BC_cum2$knot <- "4"
OM_cum2$knot <- "4"

#================== lag knots 4 =====================
basis_PM2.5 <- crossbasis(lag_PM2.5, lag = c(slag,nlag),
                          argvar = list(fun = "ns", knots = quantile(data$PM2.5, c(0.25, 0.5, 0.75))),
                          arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_SO4 <- crossbasis(lag_SO4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$SO4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_NO3 <- crossbasis(lag_NO3, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NO3, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_NH4 <- crossbasis(lag_NH4, c(slag,nlag),
                        argvar = list(fun = "ns", knots = quantile(data$NH4, c(0.25, 0.5, 0.75))),
                        arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_BC <- crossbasis(lag_BC, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$BC, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

basis_OM <- crossbasis(lag_OM, c(slag,nlag),
                       argvar = list(fun = "ns", knots = quantile(data$OM, c(0.25, 0.5, 0.75))),
                       arglag = list(fun = "ns", knots = equalknots(slag:nlag, 5)))

colnames(basis_PM2.5) = paste0("basis_PM2.5.", colnames(basis_PM2.5))
colnames(basis_SO4) = paste0("basis_SO4.", colnames(basis_SO4))
colnames(basis_NO3) = paste0("basis_NO3.", colnames(basis_NO3))
colnames(basis_NH4) = paste0("basis_NH4.", colnames(basis_NH4))
colnames(basis_BC) = paste0("basis_BC.", colnames(basis_BC))
colnames(basis_OM) = paste0("basis_OM.", colnames(basis_OM))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum3 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))


#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum3 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum3 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum3 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum3 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum3 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum3$Type <- "PM2.5"
SO4_cum3$Type <- "Sulfate"
NO3_cum3$Type <- "Nitrate"
NH4_cum3$Type <- "Ammonium"
BC_cum3$Type <- "BC"
OM_cum3$Type <- "OM"

PM2.5_cum3$knot <- "5"
SO4_cum3$knot <- "5"
NO3_cum3$knot <- "5"
NH4_cum3$knot <- "5"
BC_cum3$knot <- "5"
OM_cum3$knot <- "5"

lagknots_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, PM2.5_cum3,
                                 SO4_cum1, SO4_cum2, SO4_cum3,
                                 NO3_cum1, NO3_cum2, NO3_cum3,
                                 NH4_cum1, NH4_cum2, NH4_cum3,
                                 BC_cum1, BC_cum2, BC_cum3,
                                 OM_cum1, OM_cum2, OM_cum3)

#======================== Sensitivity Analysis 2 ==============================

#================== DLNM - Lag =====================
slag = -3
nlag = 12

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)

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


#================== Df of MT 2 =====================
ns.MT <- ns(data$MT, 2)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
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

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$Type <- "PM2.5"
SO4_mat1$Type <- "Sulfate"
NO3_mat1$Type <- "Nitrate"
NH4_mat1$Type <- "Ammonium"
BC_mat1$Type <- "BC"
OM_mat1$Type <- "OM"

PM2.5_mat1$var50 <- "40"
SO4_mat1$var50 <- "8.5"
NO3_mat1$var50 <- "9"
NH4_mat1$var50 <- "7"
BC_mat1$var50 <- "2"
OM_mat1$var50 <- "10"

PM2.5_mat1$var75 <- "60"
SO4_mat1$var75 <- "12"
NO3_mat1$var75 <- "16"
NH4_mat1$var75 <- "10.5"
BC_mat1$var75 <- "2.8"
OM_mat1$var75 <- "14"

PM2.5_mat1$var95 <- "95"
SO4_mat1$var95 <- "17"
NO3_mat1$var95 <- "25"
NH4_mat1$var95 <- "16"
BC_mat1$var95 <- "4.8"
OM_mat1$var95 <- "24"

PM2.5_mat1$df <- "2"
SO4_mat1$df <- "2"
NO3_mat1$df <- "2"
NH4_mat1$df <- "2"
BC_mat1$df <- "2"
OM_mat1$df <- "2"

#================== Df of MT 4 =====================
ns.MT <- ns(data$MT, 4)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
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

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat2$Type <- "PM2.5"
SO4_mat2$Type <- "Sulfate"
NO3_mat2$Type <- "Nitrate"
NH4_mat2$Type <- "Ammonium"
BC_mat2$Type <- "BC"
OM_mat2$Type <- "OM"

PM2.5_mat2$var50 <- "40"
SO4_mat2$var50 <- "8.5"
NO3_mat2$var50 <- "9"
NH4_mat2$var50 <- "7"
BC_mat2$var50 <- "2"
OM_mat2$var50 <- "10"

PM2.5_mat2$var75 <- "60"
SO4_mat2$var75 <- "12"
NO3_mat2$var75 <- "16"
NH4_mat2$var75 <- "10.5"
BC_mat2$var75 <- "2.8"
OM_mat2$var75 <- "14"

PM2.5_mat2$var95 <- "95"
SO4_mat2$var95 <- "17"
NO3_mat2$var95 <- "25"
NH4_mat2$var95 <- "16"
BC_mat2$var95 <- "4.8"
OM_mat2$var95 <- "24"

PM2.5_mat2$df <- "4"
SO4_mat2$df <- "4"
NO3_mat2$df <- "4"
NH4_mat2$df <- "4"
BC_mat2$df <- "4"
OM_mat2$df <- "4"

#================== Df of MT 5 =====================
ns.MT <- ns(data$MT, 5)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat3 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
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

###########contour map###########
SO4_mat3<- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                      RR = as.numeric(result.SO4$matRRfit), 
                      RR.low = as.numeric(result.SO4$matRRlow), 
                      RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat3 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat3 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat3 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat3 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat3$Type <- "PM2.5"
SO4_mat3$Type <- "Sulfate"
NO3_mat3$Type <- "Nitrate"
NH4_mat3$Type <- "Ammonium"
BC_mat3$Type <- "BC"
OM_mat3$Type <- "OM"

PM2.5_mat3$var50 <- "40"
SO4_mat3$var50 <- "8.5"
NO3_mat3$var50 <- "9"
NH4_mat3$var50 <- "7"
BC_mat3$var50 <- "2"
OM_mat3$var50 <- "10"

PM2.5_mat3$var75 <- "60"
SO4_mat3$var75 <- "12"
NO3_mat3$var75 <- "16"
NH4_mat3$var75 <- "10.5"
BC_mat3$var75 <- "2.8"
OM_mat3$var75 <- "14"

PM2.5_mat3$var95 <- "95"
SO4_mat3$var95 <- "17"
NO3_mat3$var95 <- "25"
NH4_mat3$var95 <- "16"
BC_mat3$var95 <- "4.8"
OM_mat3$var95 <- "24"

PM2.5_mat3$df <- "5"
SO4_mat3$df <- "5"
NO3_mat3$df <- "5"
NH4_mat3$df <- "5"
BC_mat3$df <- "5"
OM_mat3$df <- "5"

MTdf_contour <- bind_rows(PM2.5_mat1, PM2.5_mat2, PM2.5_mat3,
                          SO4_mat1, SO4_mat2, SO4_mat3,
                          NO3_mat1, NO3_mat2, NO3_mat3,
                          NH4_mat1, NH4_mat2, NH4_mat3,
                          BC_mat1, BC_mat2, BC_mat3,
                          OM_mat1, OM_mat2, OM_mat3)

#================== DLNM - Exposure =====================
slag = 0
nlag = 12

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)

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

#================== Df of MT 2 ====================
ns.MT <- ns(data$MT, 2)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$Type <- "PM2.5"
SO4_cum1$Type <- "Sulfate"
NO3_cum1$Type <- "Nitrate"
NH4_cum1$Type <- "Ammonium"
BC_cum1$Type <- "BC"
OM_cum1$Type <- "OM"

PM2.5_cum1$var50 <- "40"
SO4_cum1$var50 <- "8.5"
NO3_cum1$var50 <- "9"
NH4_cum1$var50 <- "7"
BC_cum1$var50 <- "2"
OM_cum1$var50 <- "10"

PM2.5_cum1$var75 <- "60"
SO4_cum1$var75 <- "12"
NO3_cum1$var75 <- "16"
NH4_cum1$var75 <- "10.5"
BC_cum1$var75 <- "2.8"
OM_cum1$var75 <- "14"

PM2.5_cum1$var95 <- "95"
SO4_cum1$var95 <- "17"
NO3_cum1$var95 <- "25"
NH4_cum1$var95 <- "16"
BC_cum1$var95 <- "4.8"
OM_cum1$var95 <- "24"

PM2.5_cum1$df <- "2"
SO4_cum1$df <- "2"
NO3_cum1$df <- "2"
NH4_cum1$df <- "2"
BC_cum1$df <- "2"
OM_cum1$df <- "2"

#================== Df of MT 4 ====================
ns.MT <- ns(data$MT, 4)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 


NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM ===================== 
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum2$Type <- "PM2.5"
SO4_cum2$Type <- "Sulfate"
NO3_cum2$Type <- "Nitrate"
NH4_cum2$Type <- "Ammonium"
BC_cum2$Type <- "BC"
OM_cum2$Type <- "OM"

PM2.5_cum2$var50 <- "40"
SO4_cum2$var50 <- "8.5"
NO3_cum2$var50 <- "9"
NH4_cum2$var50 <- "7"
BC_cum2$var50 <- "2"
OM_cum2$var50 <- "10"

PM2.5_cum2$var75 <- "60"
SO4_cum2$var75 <- "12"
NO3_cum2$var75 <- "16"
NH4_cum2$var75 <- "10.5"
BC_cum2$var75 <- "2.8"
OM_cum2$var75 <- "14"

PM2.5_cum2$var95 <- "95"
SO4_cum2$var95 <- "17"
NO3_cum2$var95 <- "25"
NH4_cum2$var95 <- "16"
BC_cum2$var95 <- "4.8"
OM_cum2$var95 <- "24"

PM2.5_cum2$df <- "4"
SO4_cum2$df <- "4"
NO3_cum2$df <- "4"
NH4_cum2$df <- "4"
BC_cum2$df <- "4"
OM_cum2$df <- "4"

#================== Df of MT 5 ====================
ns.MT <- ns(data$MT, 5)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum3 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum3<- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                      RR = as.numeric(result.SO4$cumRRfit), 
                      RR.low = as.numeric(result.SO4$cumRRlow), 
                      RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum3 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum3 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum3 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum3 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum3$Type <- "PM2.5"
SO4_cum3$Type <- "Sulfate"
NO3_cum3$Type <- "Nitrate"
NH4_cum3$Type <- "Ammonium"
BC_cum3$Type <- "BC"
OM_cum3$Type <- "OM"

PM2.5_cum3$var50 <- "40"
SO4_cum3$var50 <- "8.5"
NO3_cum3$var50 <- "9"
NH4_cum3$var50 <- "7"
BC_cum3$var50 <- "2"
OM_cum3$var50 <- "10"

PM2.5_cum3$var75 <- "60"
SO4_cum3$var75 <- "12"
NO3_cum3$var75 <- "16"
NH4_cum3$var75 <- "10.5"
BC_cum3$var75 <- "2.8"
OM_cum3$var75 <- "14"

PM2.5_cum3$var95 <- "95"
SO4_cum3$var95 <- "17"
NO3_cum3$var95 <- "25"
NH4_cum3$var95 <- "16"
BC_cum3$var95 <- "4.8"
OM_cum3$var95 <- "24"

PM2.5_cum3$df <- "5"
SO4_cum3$df <- "5"
NO3_cum3$df <- "5"
NH4_cum3$df <- "5"
BC_cum3$df <- "5"
OM_cum3$df <- "5"

MTdf_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, PM2.5_cum3,
                             SO4_cum1, SO4_cum2, SO4_cum3,
                             NO3_cum1, NO3_cum2, NO3_cum3,
                             NH4_cum1, NH4_cum2, NH4_cum3,
                             BC_cum1, BC_cum2, BC_cum3,
                             OM_cum1, OM_cum2, OM_cum3)


#======================== Sensitivity Analysis 3 ==============================

#================== DLNM - Lag =====================
#================== Max Lag 6 =====================
slag = -3
nlag = 6 

ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)


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
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:6, each = nrow(result.PM2.5$matRRfit)), 
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

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:6, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:6, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:6, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:6, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:6, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$Type <- "PM2.5"
SO4_mat1$Type <- "Sulfate"
NO3_mat1$Type <- "Nitrate"
NH4_mat1$Type <- "Ammonium"
BC_mat1$Type <- "BC"
OM_mat1$Type <- "OM"

PM2.5_mat1$var50 <- "40"
SO4_mat1$var50 <- "8.5"
NO3_mat1$var50 <- "9"
NH4_mat1$var50 <- "7"
BC_mat1$var50 <- "2"
OM_mat1$var50 <- "10"

PM2.5_mat1$var75 <- "60"
SO4_mat1$var75 <- "12"
NO3_mat1$var75 <- "16"
NH4_mat1$var75 <- "10.5"
BC_mat1$var75 <- "2.8"
OM_mat1$var75 <- "14"

PM2.5_mat1$var95 <- "95"
SO4_mat1$var95 <- "17"
NO3_mat1$var95 <- "25"
NH4_mat1$var95 <- "16"
BC_mat1$var95 <- "4.8"
OM_mat1$var95 <- "24"

PM2.5_mat1$nlag <- "6"
SO4_mat1$nlag <- "6"
NO3_mat1$nlag <- "6"
NH4_mat1$nlag <- "6"
BC_mat1$nlag <- "6"
OM_mat1$nlag <- "6"

#================== Max Lag 9 =====================
slag = -3
nlag = 9

ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)


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
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:9, each = nrow(result.PM2.5$matRRfit)), 
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

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:9, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:9, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:9, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:9, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:9, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat2$Type <- "PM2.5"
SO4_mat2$Type <- "Sulfate"
NO3_mat2$Type <- "Nitrate"
NH4_mat2$Type <- "Ammonium"
BC_mat2$Type <- "BC"
OM_mat2$Type <- "OM"

PM2.5_mat2$var50 <- "40"
SO4_mat2$var50 <- "8.5"
NO3_mat2$var50 <- "9"
NH4_mat2$var50 <- "7"
BC_mat2$var50 <- "2"
OM_mat2$var50 <- "10"

PM2.5_mat2$var75 <- "60"
SO4_mat2$var75 <- "12"
NO3_mat2$var75 <- "16"
NH4_mat2$var75 <- "10.5"
BC_mat2$var75 <- "2.8"
OM_mat2$var75 <- "14"

PM2.5_mat2$var95 <- "95"
SO4_mat2$var95 <- "17"
NO3_mat2$var95 <- "25"
NH4_mat2$var95 <- "16"
BC_mat2$var95 <- "4.8"
OM_mat2$var95 <- "24"

PM2.5_mat2$nlag <- "9"
SO4_mat2$nlag <- "9"
NO3_mat2$nlag <- "9"
NH4_mat2$nlag <- "9"
BC_mat2$nlag <- "9"
OM_mat2$nlag <- "9"

nlag_contour <- bind_rows(PM2.5_mat1, PM2.5_mat2, 
                          SO4_mat1, SO4_mat2, 
                          NO3_mat1, NO3_mat2,
                          NH4_mat1, NH4_mat2, 
                          BC_mat1, BC_mat2,
                          OM_mat1, OM_mat2)


#================== DLNM - Exposure =====================
#================== Max Lag 6 =====================
slag = 0
nlag = 6

ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)


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
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:6, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:6, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:6, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:6, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:6, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:6, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$Type <- "PM2.5"
SO4_cum1$Type <- "Sulfate"
NO3_cum1$Type <- "Nitrate"
NH4_cum1$Type <- "Ammonium"
BC_cum1$Type <- "BC"
OM_cum1$Type <- "OM"

PM2.5_cum1$var50 <- "40"
SO4_cum1$var50 <- "8.5"
NO3_cum1$var50 <- "9"
NH4_cum1$var50 <- "7"
BC_cum1$var50 <- "2"
OM_cum1$var50 <- "10"

PM2.5_cum1$var75 <- "60"
SO4_cum1$var75 <- "12"
NO3_cum1$var75 <- "16"
NH4_cum1$var75 <- "10.5"
BC_cum1$var75 <- "2.8"
OM_cum1$var75 <- "14"

PM2.5_cum1$var95 <- "95"
SO4_cum1$var95 <- "17"
NO3_cum1$var95 <- "25"
NH4_cum1$var95 <- "16"
BC_cum1$var95 <- "4.8"
OM_cum1$var95 <- "24"

PM2.5_cum1$nlag <- "6"
SO4_cum1$nlag <- "6"
NO3_cum1$nlag <- "6"
NH4_cum1$nlag <- "6"
BC_cum1$nlag <- "6"
OM_cum1$nlag <- "6"

#================== Max Lag 9 =====================
slag = 0
nlag = 9

ns.MT <- ns(data$MT, 3)
colnames(ns.MT) <- paste("MT.", colnames(ns.MT), sep = "")

#lag
lag_PM2.5 <- tsModel::Lag(data$PM2.5, group = data$regnames, k = slag:nlag)
lag_SO4 <- tsModel::Lag(data$SO4, group = data$regnames, k = slag:nlag)
lag_NO3 <- tsModel::Lag(data$NO3, group = data$regnames, k = slag:nlag)
lag_NH4 <- tsModel::Lag(data$NH4, group = data$regnames, k = slag:nlag)
lag_BC <- tsModel::Lag(data$BC, group = data$regnames, k = slag:nlag)
lag_OM <- tsModel::Lag(data$OM, group = data$regnames, k = slag:nlag)


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
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:9, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
model2 <- fitmodel(model2, data)     

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25)) 

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:9, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))


#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:9, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:9, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:9, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM ===================== 
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:9, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum2$Type <- "PM2.5"
SO4_cum2$Type <- "Sulfate"
NO3_cum2$Type <- "Nitrate"
NH4_cum2$Type <- "Ammonium"
BC_cum2$Type <- "BC"
OM_cum2$Type <- "OM"

PM2.5_cum2$var50 <- "40"
SO4_cum2$var50 <- "8.5"
NO3_cum2$var50 <- "9"
NH4_cum2$var50 <- "7"
BC_cum2$var50 <- "2"
OM_cum2$var50 <- "10"

PM2.5_cum2$var75 <- "60"
SO4_cum2$var75 <- "12"
NO3_cum2$var75 <- "16"
NH4_cum2$var75 <- "10.5"
BC_cum2$var75 <- "2.8"
OM_cum2$var75 <- "14"

PM2.5_cum2$var95 <- "95"
SO4_cum2$var95 <- "17"
NO3_cum2$var95 <- "25"
NH4_cum2$var95 <- "16"
BC_cum2$var95 <- "4.8"
OM_cum2$var95 <- "24"

PM2.5_cum2$nlag <- "9"
SO4_cum2$nlag <- "9"
NO3_cum2$nlag <- "9"
NH4_cum2$nlag <- "9"
BC_cum2$nlag <- "9"
OM_cum2$nlag <- "9"

nlag_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, 
                             SO4_cum1, SO4_cum2, 
                             NO3_cum1, NO3_cum2,
                             NH4_cum1, NH4_cum2, 
                             BC_cum1, BC_cum2, 
                             OM_cum1, OM_cum2)

#======================== Sensitivity Analysis 4 ==============================

#================== DLNM - Lag =====================
#================== -lag 3 =====================
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
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_PM2.5

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_SO4

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NO3

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_NH4

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  ns.MT +
  basis_BC

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
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

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-6:12, each = nrow(result.PM2.5$matRRfit)), 
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

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-6:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
model3 <- fitmodel(model3, data)     

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25)) 

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-6:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
model4 <- fitmodel(model4, data)     

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-6:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
model5 <- fitmodel(model5, data)     

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25)) 

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-6:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
model6 <- fitmodel(model6, data)     

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25)) 

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-6:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$Type <- "PM2.5"
SO4_mat1$Type <- "Sulfate"
NO3_mat1$Type <- "Nitrate"
NH4_mat1$Type <- "Ammonium"
BC_mat1$Type <- "BC"
OM_mat1$Type <- "OM"

PM2.5_mat1$var50 <- "40"
SO4_mat1$var50 <- "8.5"
NO3_mat1$var50 <- "9"
NH4_mat1$var50 <- "7"
BC_mat1$var50 <- "2"
OM_mat1$var50 <- "10"

PM2.5_mat1$var75 <- "60"
SO4_mat1$var75 <- "12"
NO3_mat1$var75 <- "16"
NH4_mat1$var75 <- "10.5"
BC_mat1$var75 <- "2.8"
OM_mat1$var75 <- "14"

PM2.5_mat1$var95 <- "95"
SO4_mat1$var95 <- "17"
NO3_mat1$var95 <- "25"
NH4_mat1$var95 <- "16"
BC_mat1$var95 <- "4.8"
OM_mat1$var95 <- "24"

PM2.5_mat1$nlag <- "-6 ~ 12"
SO4_mat1$nlag <- "-6 ~ 12"
NO3_mat1$nlag <- "-6 ~ 12"
NH4_mat1$nlag <- "-6 ~ 12"
BC_mat1$nlag <- "-6 ~ 12"
OM_mat1$nlag <- "-6 ~ 12"

nlag_contour <- bind_rows(PM2.5_mat1,
                          SO4_mat1, 
                          NO3_mat1,
                          NH4_mat1,  
                          BC_mat1, 
                          OM_mat1)


#================== Figures export =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM) 

AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM)

data_pollution <- rbind(AH.data, SC.data)

#================== Sensitivity Analysis 1 =====================
PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- lagknots_contour   %>%
  filter(Type == "PM2.5") %>%
  group_by(knot) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = knot) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.24)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ knot, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab((expression(PM[2.5] ~ "(μg/m³)"))) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

SO4_concentration <- data_pollution$SO4

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4_contour <- lagknots_contour   %>%
  filter(Type == "Sulfate") %>%
  group_by(knot) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = knot) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.24)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ knot, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Sulfate (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NO3_concentration <- data_pollution$NO3

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3_contour <- lagknots_contour   %>%
  filter(Type == "Nitrate") %>%
  group_by(knot) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = knot) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.24)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ knot, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Nitrate (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NH4_concentration <- data_pollution$NH4

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4_contour <- lagknots_contour   %>%
  filter(Type == "Ammonium") %>%
  group_by(knot) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = knot) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.24)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ knot, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Ammonium (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

OM_concentration <- data_pollution$OM

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM_contour <- lagknots_contour   %>%
  filter(Type == "OM") %>%
  group_by(knot) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = knot) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.24)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ knot, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("OM (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

BC_concentration <- data_pollution$BC

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC_contour <- lagknots_contour   %>%
  filter(Type == "BC") %>%
  group_by(knot) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = knot) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.24)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ knot, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("BC (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

All_contour <- plot_grid(PM2.5_contour, SO4_contour, NO3_contour, NH4_contour, OM_contour, BC_contour, nrow = 3)

pdf("./Figures/Figure S22.pdf", height = 8, width = 14)
All_contour 
dev.off()

lagknots_contour <- lagknots_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) # 其他值保持不变
  ))

lagknots_contour <- lagknots_contour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#a50f15", "#2171b5", "#238b45")

Lag_knot <- lagknots_contour %>%
  group_by(knot, Type, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ Type, scales = "fixed") +
  theme(legend.position = "bottom", 
        legend.justification = c(0, 1),
        plot.title = element_text(hjust = 0, size = 10, face = "bold")) +
  labs(x = "Lag (months)", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_blank())

Lag_knot <- ggdraw(Lag_knot) +
  annotate("text", x = 0.125, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.435, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure S23.pdf", height = 7, width = 14)
Lag_knot 
dev.off()

lagknots_cumcontour <- lagknots_cumcontour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c( "#008B45FF", "#3B4992FF","#EE0000FF")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- lagknots_cumcontour %>%
  filter( Type  == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(knot,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.45)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

SO4_concentration <- data_pollution$SO4[data_pollution$SO4 <= 17.5]

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4 <- lagknots_cumcontour %>%
  filter(Type  == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(knot, Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.45)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

NO3_concentration <- data_pollution$NO3[data_pollution$NO3 <= 26]

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3 <- lagknots_cumcontour %>%
  filter(Type  == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(knot,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.45)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

NH4_concentration <- data_pollution$NH4[data_pollution$NH4 <= 16.5]

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4 <- lagknots_cumcontour %>%
  filter( Type  == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(knot,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.45)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

OM_concentration <- data_pollution$OM[data_pollution$OM <= 25]

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM <- lagknots_cumcontour %>%
  filter( Type  == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(knot,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.45)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

BC_concentration <- data_pollution$BC[data_pollution$BC <= 5]

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC <- lagknots_cumcontour %>%
  filter( Type  == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(knot,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = knot, col = as.factor(knot), fill = as.factor(knot))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.45)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

legend <- get_legend(
  PM2.5 + 
    theme(legend.position = "right",  
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0.2, "cm")) + 
    guides(color = guide_legend(ncol = 1)))

exposure <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure <- ggdraw(exposure) +
  annotate("text", x = 0.2, y = 0.985, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 5, hjust = 0.5) +
  annotate("text", x = 0.54, y = 0.985, label = "B.Sulfate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.86, y = 0.985, label = "C.Nitrate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.2, y = 0.5, label = "D.Ammonium (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.53, y = 0.5, label = "E.OM (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.86, y = 0.5, label = "F.BC (μg/m³)", size = 5, hjust = 0.5)

final_plot <- ggdraw() +
  draw_plot(exposure, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.92,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S24.pdf", height = 7, width = 14)
final_plot
dev.off()

#================== Sensitivity Analysis 2 =====================
PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- MTdf_contour   %>%
  filter(Type == "PM2.5") %>%
  group_by(df) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = df) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.19)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ df, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab((expression(PM[2.5] ~ "(μg/m³)"))) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

SO4_concentration <- data_pollution$SO4

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4_contour <- MTdf_contour   %>%
  filter(Type == "Sulfate") %>%
  group_by(df) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = df) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.19)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ df, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Sulfate (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NO3_concentration <- data_pollution$NO3

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3_contour <- MTdf_contour   %>%
  filter(Type == "Nitrate") %>%
  group_by(df) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = df) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.19)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ df, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Nitrate (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NH4_concentration <- data_pollution$NH4

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4_contour <- MTdf_contour   %>%
  filter(Type == "Ammonium") %>%
  group_by(df) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = df) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.19)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ df, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Ammonium (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

OM_concentration <- data_pollution$OM

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM_contour <- MTdf_contour   %>%
  filter(Type == "OM") %>%
  group_by(df) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = df) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.19)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ df, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("OM (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

BC_concentration <- data_pollution$BC

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC_contour <- MTdf_contour   %>%
  filter(Type == "BC") %>%
  group_by(df) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = df) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.19)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ df, scales = "fixed", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("BC (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

All_contour <- plot_grid(PM2.5_contour, SO4_contour, NO3_contour, NH4_contour, OM_contour, BC_contour, nrow = 3)

pdf("./Figures/Figure S25.pdf", height = 8, width = 14)
All_contour 
dev.off()

MTdf_contour <- MTdf_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) # 其他值保持不变
  ))

MTdf_contour <- MTdf_contour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#a50f15", "#2171b5", "#238b45")

Lag_df <- MTdf_contour %>%
  group_by(df, Type, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ Type, scales = "fixed") +
  theme(legend.position = "bottom", 
        legend.justification = c(0, 1),
        plot.title = element_text(hjust = 0, size = 10, face = "bold")) +
  labs(x = "Lag (months)", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_blank())

Lag_df <- ggdraw(Lag_df) +
  annotate("text", x = 0.125, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.435, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure S26.pdf", height = 7, width = 14)
Lag_df
dev.off()

MTdf_cumcontour <- MTdf_cumcontour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c( "#008B45FF", "#3B4992FF","#EE0000FF")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- MTdf_cumcontour %>%
  filter( Type  == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(df,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

SO4_concentration <- data_pollution$SO4[data_pollution$SO4 <= 17.5]

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4 <- MTdf_cumcontour %>%
  filter(Type  == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(df, Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

NO3_concentration <- data_pollution$NO3[data_pollution$NO3 <= 26]

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3 <- MTdf_cumcontour %>%
  filter(Type  == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(df,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

NH4_concentration <- data_pollution$NH4[data_pollution$NH4 <= 16.5]

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4 <- MTdf_cumcontour %>%
  filter(Type  == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(df,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

OM_concentration <- data_pollution$OM[data_pollution$OM <= 25]

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM <- MTdf_cumcontour %>%
  filter(Type  == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(df,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

BC_concentration <- data_pollution$BC[data_pollution$BC <= 5]

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC <- MTdf_cumcontour %>%
  filter(Type  == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(df,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = df, col = as.factor(df), fill = as.factor(df))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

legend <- get_legend(
  PM2.5 + 
    theme(legend.position = "right",  
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0.2, "cm")) + 
    guides(color = guide_legend(ncol = 1)))

exposure <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure <- ggdraw(exposure) +
  annotate("text", x = 0.2, y = 0.985, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 5, hjust = 0.5) +
  annotate("text", x = 0.54, y = 0.985, label = "B.Sulfate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.86, y = 0.985, label = "C.Nitrate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.2, y = 0.5, label = "D.Ammonium (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.53, y = 0.5, label = "E.OM (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.86, y = 0.5, label = "F.BC (μg/m³)", size = 5, hjust = 0.5)

final_plot <- ggdraw() +
  draw_plot(exposure, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.92,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S27.pdf", height = 7, width = 14)
final_plot
dev.off()

#================== Sensitivity Analysis 3 =====================
PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- nlag_contour   %>%
  filter(Type == "PM2.5") %>%
  group_by(nlag) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = nlag) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.16)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ nlag, scales = "free", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab((expression(PM[2.5] ~ "(μg/m³)"))) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

SO4_concentration <- data_pollution$SO4

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4_contour <- nlag_contour   %>%
  filter(Type == "Sulfate") %>%
  group_by(nlag) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = nlag) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.16)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ nlag, scales = "free", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Sulfate (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NO3_concentration <- data_pollution$NO3

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3_contour <- nlag_contour   %>%
  filter(Type == "Nitrate") %>%
  group_by(nlag) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = nlag) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.16)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ nlag, scales = "free", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Nitrate (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NH4_concentration <- data_pollution$NH4

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4_contour <- nlag_contour   %>%
  filter(Type == "Ammonium") %>%
  group_by(nlag) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = nlag) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.16)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ nlag, scales = "free", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("Ammonium (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

OM_concentration <- data_pollution$OM

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM_contour <- nlag_contour   %>%
  filter(Type == "OM") %>%
  group_by(nlag) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = nlag) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.16)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ nlag, scales = "free", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("OM (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

BC_concentration <- data_pollution$BC

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC_contour <- nlag_contour   %>%
  filter(Type == "BC") %>%
  group_by(nlag) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = nlag) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.68, 1.16)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ nlag, scales = "free", nrow = 1) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90,
                                   hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("BC (μg/m³)")  +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

All_contour <- plot_grid(PM2.5_contour, SO4_contour, NO3_contour, NH4_contour, OM_contour, BC_contour, nrow = 3)

pdf("./Figures/Figure S28.pdf", height = 10, width = 14)
All_contour 
dev.off()

nlag_contour <- nlag_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) 
  ))

nlag_contour <- nlag_contour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#000000", pal_lancet()(5))

nlag_df <- nlag_contour %>%
  group_by(nlag, Type, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = Type, col = as.factor(Type), fill = as.factor(Type))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal, labels = c(expression(paste(PM[2.5])), "Sulfate", "Nitrate", "Ammonium", "BC", "OM")) +
  scale_fill_manual(values = col.pal, guide = FALSE) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9)) +
  facet_grid(VariableValue ~ nlag, scales = "free") +
  theme(legend.position = "bottom", 
        legend.justification = c(0, 1),
        plot.title = element_text(hjust = 0, size = 10, face = "bold")) +
  labs(x = "Lag (months)", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1)) +
  guides(col = guide_legend(nrow = 1, byrow = TRUE)) 

pdf("./Figures/Figure S29.pdf", height = 7, width = 14)
nlag_df
dev.off()

nlag_cumcontour <- nlag_cumcontour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#3B4992FF","#EE0000FF")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- nlag_cumcontour %>%
  filter(Type  == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(nlag,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = nlag, col = as.factor(nlag), fill = as.factor(nlag))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

SO4_concentration <- data_pollution$SO4[data_pollution$SO4 <= 17.5]

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4 <- nlag_cumcontour %>%
  filter(Type  == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(nlag, Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = nlag, col = as.factor(nlag), fill = as.factor(nlag))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

NO3_concentration <- data_pollution$NO3[data_pollution$NO3 <= 26]

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3 <- nlag_cumcontour %>%
  filter(Type  == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(nlag,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = nlag, col = as.factor(nlag), fill = as.factor(nlag))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

NH4_concentration <- data_pollution$NH4[data_pollution$NH4 <= 16.5]

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4 <- nlag_cumcontour %>%
  filter(Type  == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(nlag,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = nlag, col = as.factor(nlag), fill = as.factor(nlag))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

OM_concentration <- data_pollution$OM[data_pollution$OM <= 25]

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM <- nlag_cumcontour %>%
  filter(Type  == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(nlag,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = nlag, col = as.factor(nlag), fill = as.factor(nlag))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

BC_concentration <- data_pollution$BC[data_pollution$BC <= 5]

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC <- nlag_cumcontour %>%
  filter(Type  == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(nlag,  Type ) %>%
  ggplot(aes(x = VariableValue, y = RR, group = nlag, col = as.factor(nlag), fill = as.factor(nlag))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~ Type , scales = "free") +
  scale_y_continuous(limits = c(0.7, 1.5)) +
  theme(legend.position = c(0.01, 0.9)) +
  labs(x = "", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "b", data = regular_data, aes(x = VariableValue), inherit.aes = FALSE, color = "black") +
  geom_rug(sides = "b", data = special_data, aes(x = VariableValue), inherit.aes = FALSE,color = "red") +
  geom_rug(sides = "b", data = green_rug_data, aes(x = VariableValue), inherit.aes = FALSE, color = "green")

legend <- get_legend(
  PM2.5 + 
    theme(legend.position = "right",  
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0.2, "cm")) + 
    guides(color = guide_legend(ncol = 1)))

exposure <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure <- ggdraw(exposure) +
  annotate("text", x = 0.2, y = 0.985, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 5, hjust = 0.5) +
  annotate("text", x = 0.54, y = 0.985, label = "B.Sulfate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.86, y = 0.985, label = "C.Nitrate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.2, y = 0.5, label = "D.Ammonium (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.53, y = 0.5, label = "E.OM (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.86, y = 0.5, label = "F.BC (μg/m³)", size = 5, hjust = 0.5)

final_plot <- ggdraw() +
  draw_plot(exposure, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.92,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S30.pdf", height = 7, width = 14)
final_plot
dev.off()

#================== Sensitivity Analysis 4  =====================
PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- nlag_contour %>%
  filter(Type == "PM2.5") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90, hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("") + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

SO4_concentration <- data_pollution$SO4

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4 <- nlag_contour %>%
  filter(Type == "Sulfate") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90, hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("") + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NO3_concentration <- data_pollution$NO3

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3 <- nlag_contour %>%
  filter(Type == "Nitrate") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90, hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("") + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

NH4_concentration <- data_pollution$NH4

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4 <- nlag_contour %>%
  filter(Type == "Ammonium") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90, hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("") + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

OM_concentration <- data_pollution$OM

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM <- nlag_contour %>%
  filter(Type == "OM") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90, hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("") + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

BC_concentration <- data_pollution$BC

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC <- nlag_contour %>%
  filter(Type == "BC") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.line = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(angle = 90, hjust = 0.5)) +
  xlab("Lag (months)") +
  ylab("") + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank()) +
  geom_rug(sides = "r", data = regular_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "black") +
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red") +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green")

count1 <- plot_grid(PM2.5, SO4, NO3, nrow = 1)

count2 <- plot_grid(NH4, OM, BC, nrow = 1)

count3 <- ggdraw() +
  draw_plot(count1, x = 0, y = 0, width = 1, height = 0.95) + 
  annotate("text", x = 0.17, y = 0.955, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 5, hjust = 0.5) +
  annotate("text", x = 0.49, y = 0.955, label = "B.Sulfate (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.83, y = 0.955, label = "C.Nitrate (μg/m³)", size = 5, hjust = 0.5) 

count4 <- ggdraw() +
  draw_plot(count2, x = 0, y = 0, width = 1, height = 0.95) +
  annotate("text", x = 0.16, y = 0.955, label = "D.Ammonium (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.49, y = 0.955, label = "E.OM (μg/m³)", size = 5, hjust = 0.5) +
  annotate("text", x = 0.83, y = 0.955, label = "F.BC (μg/m³)", size = 5, hjust = 0.5)

count <- plot_grid(count3, count4, nrow = 2)

pdf("./Figures/Figure S31.pdf", height = 6, width = 12)
count 
dev.off()

nlag_contour <- nlag_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) 
  ))

nlag_contour <- nlag_contour %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#000000", pal_lancet()(5))

Lag = nlag_contour %>%
  group_by(VariableValue, Type) %>%
  ggplot(aes(x = Lag, y = RR, group = Type, col = as.factor(Type), fill = as.factor(Type) )) +
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal, labels = c(expression(paste(PM[2.5])), "Sulfate", "Nitrate", "Ammonium", "BC", "OM")) +
  scale_fill_manual(values = col.pal, guide = FALSE) +
  theme(legend.position = c(0.01, 0.9)) +
  facet_wrap(~VariableValue, scales = "fixed", nrow = 3) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3, 6, 9, 12)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Lag (months)", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1)) +
  guides(col = guide_legend(nrow = 1, byrow = TRUE))


pdf("./Figures/Figure S32.pdf", height = 8, width = 12)
Lag 
dev.off()
