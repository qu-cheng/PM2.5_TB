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
  dplyr::select(year, month, X, man, woman, age1, age2, age3, age4, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop1, pop2, pop3, pop4, pop5, pop6, regnames, MT) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(year, month, X, man, woman, age1, age2, age3, age4, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop1, pop2, pop3, pop4, pop5, pop6, regnames, MT)%>%
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
         E1 = pop1*10000,
         E2 = pop2*10000,
         E3 = pop3*10000,
         E4 = pop4*10000,
         E5 = pop5*10000,
         E6 = pop6*10000)

#======================== DLNM - Lag ==============================
fitmodel <- function(formula, data = data, family = "nbinomial")  # GZ.daily.subset is a subset of GZ.daily in dengue transmission years and months
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

#================== PM2.5 =====================
PM2.5_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

##man
PM2.5_model1 <- fitmodel(PM2.5_model1, data) 

coef <- PM2.5_model1$summary.fixed$mean

vcov <- PM2.5_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model1$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

##woman
PM2.5_model2 <- fitmodel(PM2.5_model2, data) 

coef <- PM2.5_model2$summary.fixed$mean

vcov <- PM2.5_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model2$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

##[0, 25)
PM2.5_model3 <- fitmodel(PM2.5_model3, data) 

coef <- PM2.5_model3$summary.fixed$mean

vcov <- PM2.5_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model3$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat3 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

##[25, 45)
PM2.5_model4 <- fitmodel(PM2.5_model4, data) 

coef <- PM2.5_model4$summary.fixed$mean

vcov <- PM2.5_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model4$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat4 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

##[45, 65)
PM2.5_model5 <- fitmodel(PM2.5_model5, data) 

coef <- PM2.5_model5$summary.fixed$mean

vcov <- PM2.5_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model5$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat5 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

##[65, ~)
PM2.5_model6 <- fitmodel(PM2.5_model6, data) 

coef <- PM2.5_model6$summary.fixed$mean

vcov <- PM2.5_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model6$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat6 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

#================== SO4 =====================
SO4_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

##man
SO4_model1 <- fitmodel(SO4_model1, data) 

coef <- SO4_model1$summary.fixed$mean

vcov <- SO4_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model1$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

##woman
SO4_model2 <- fitmodel(SO4_model2, data) 

coef <- SO4_model2$summary.fixed$mean

vcov <- SO4_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

##[0, 25)
SO4_model3 <- fitmodel(SO4_model3, data) 

coef <- SO4_model3$summary.fixed$mean

vcov <- SO4_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model3$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat3 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

##[25, 45)
SO4_model4 <- fitmodel(SO4_model4, data) 

coef <- SO4_model4$summary.fixed$mean

vcov <- SO4_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model4$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat4 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

##[45, 65)
SO4_model5 <- fitmodel(SO4_model5, data) 

coef <- SO4_model5$summary.fixed$mean

vcov <- SO4_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model5$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat5 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

##[65, ~)
SO4_model6 <- fitmodel(SO4_model6, data) 

coef <- SO4_model6$summary.fixed$mean

vcov <- SO4_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model6$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat6 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
NO3_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

##man
NO3_model1 <- fitmodel(NO3_model1, data) 

coef <- NO3_model1$summary.fixed$mean

vcov <- NO3_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model1$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

##woman
NO3_model2 <- fitmodel(NO3_model2, data) 

coef <- NO3_model2$summary.fixed$mean

vcov <- NO3_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model2$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

##[0, 25)
NO3_model3 <- fitmodel(NO3_model3, data) 

coef <- NO3_model3$summary.fixed$mean

vcov <- NO3_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat3 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

##[25, 45)
NO3_model4 <- fitmodel(NO3_model4, data) 

coef <- NO3_model4$summary.fixed$mean

vcov <- NO3_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model4$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat4 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

##[45, 65)
NO3_model5 <- fitmodel(NO3_model5, data) 

coef <- NO3_model5$summary.fixed$mean

vcov <- NO3_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model5$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat5 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

##[65, ~)
NO3_model6 <- fitmodel(NO3_model6, data) 

coef <- NO3_model6$summary.fixed$mean

vcov <- NO3_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model6$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat6 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
NH4_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

##man
NH4_model1 <- fitmodel(NH4_model1, data) 

coef <- NH4_model1$summary.fixed$mean

vcov <- NH4_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model1$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

##woman
NH4_model2 <- fitmodel(NH4_model2, data) 

coef <- NH4_model2$summary.fixed$mean

vcov <- NH4_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model2$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

##[0, 25)
NH4_model3 <- fitmodel(NH4_model3, data) 

coef <- NH4_model3$summary.fixed$mean

vcov <- NH4_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model3$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat3 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

##[25, 45)
NH4_model4 <- fitmodel(NH4_model4, data) 

coef <- NH4_model4$summary.fixed$mean

vcov <- NH4_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat4 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

##[45, 65)
NH4_model5 <- fitmodel(NH4_model5, data) 

coef <- NH4_model5$summary.fixed$mean

vcov <- NH4_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model5$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat5 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

##[65, ~)
NH4_model6 <- fitmodel(NH4_model6, data) 

coef <- NH4_model6$summary.fixed$mean

vcov <- NH4_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model6$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat6 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== OM =====================
OM_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

##man
OM_model1 <- fitmodel(OM_model1, data) 

coef <- OM_model1$summary.fixed$mean

vcov <- OM_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model1$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

##woman
OM_model2 <- fitmodel(OM_model2, data) 

coef <- OM_model2$summary.fixed$mean

vcov <- OM_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model2$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

##[0, 25)
OM_model3 <- fitmodel(OM_model3, data) 

coef <- OM_model3$summary.fixed$mean

vcov <- OM_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model3$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat3 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

##[25, 45)
OM_model4 <- fitmodel(OM_model4, data) 

coef <- OM_model4$summary.fixed$mean

vcov <- OM_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model4$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat4 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

##[45, 65)
OM_model5 <- fitmodel(OM_model5, data) 

coef <- OM_model5$summary.fixed$mean

vcov <- OM_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model5$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat5 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

##[65, ~)
OM_model6 <- fitmodel(OM_model6, data) 

coef <- OM_model6$summary.fixed$mean

vcov <- OM_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat6 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== BC =====================
BC_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

##man
BC_model1 <- fitmodel(BC_model1, data) 

coef <- BC_model1$summary.fixed$mean

vcov <- BC_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model1$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

##woman
BC_model2 <- fitmodel(BC_model2, data) 

coef <- BC_model2$summary.fixed$mean

vcov <- BC_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model2$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

##[0, 25)
BC_model3 <- fitmodel(BC_model3, data) 

coef <- BC_model3$summary.fixed$mean

vcov <- BC_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model3$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat3 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

##[25, 45)
BC_model4 <- fitmodel(BC_model4, data) 

coef <- BC_model4$summary.fixed$mean

vcov <- BC_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model4$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat4 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                        Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                        RR = as.numeric(result.BC$matRRfit), 
                        RR.low = as.numeric(result.BC$matRRlow), 
                        RR.high = as.numeric(result.BC$matRRhigh))

##[45, 65)
BC_model5 <- fitmodel(BC_model5, data) 

coef <- BC_model5$summary.fixed$mean

vcov <- BC_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat5 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

##[65, ~)
BC_model6 <- fitmodel(BC_model6, data) 

coef <- BC_model6$summary.fixed$mean

vcov <- BC_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model6$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat6 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== Results export =====================
PM2.5_mat1$group <- "Male"
PM2.5_mat2$group <- "Female"
SO4_mat1$group <- "Male"
SO4_mat2$group <- "Female"
NO3_mat1$group <- "Male"
NO3_mat2$group <- "Female"
NH4_mat1$group <- "Male"
NH4_mat2$group <- "Female"
BC_mat1$group <- "Male"
BC_mat2$group <- "Female"
OM_mat1$group <- "Male"
OM_mat2$group <- "Female"
PM2.5_mat3$group <- "[0, 25)"
PM2.5_mat4$group <- "[25, 45)"
PM2.5_mat5$group <- "[45, 65)"
PM2.5_mat6$group <- "[65, ~)"
SO4_mat3$group <- "[0, 25)"
SO4_mat4$group <- "[25, 45)"
SO4_mat5$group <- "[45, 65)"
SO4_mat6$group <- "[65, ~)"
NO3_mat3$group <- "[0, 25)"
NO3_mat4$group <- "[25, 45)"
NO3_mat5$group <- "[45, 65)"
NO3_mat6$group <- "[65, ~)"
NH4_mat3$group <- "[0, 25)"
NH4_mat4$group <- "[25, 45)"
NH4_mat5$group <- "[45, 65)"
NH4_mat6$group <- "[65, ~)"
BC_mat3$group <- "[0, 25)"
BC_mat4$group <- "[25, 45)"
BC_mat5$group <- "[45, 65)"
BC_mat6$group <- "[65, ~)"
OM_mat3$group <- "[0, 25)"
OM_mat4$group <- "[25, 45)"
OM_mat5$group <- "[45, 65)"
OM_mat6$group <- "[65, ~)"

PM2.5_mat1$pollution <- "PM2.5"
PM2.5_mat2$pollution <- "PM2.5"
SO4_mat1$pollution <- "Sulfate"
SO4_mat2$pollution <- "Sulfate"
NO3_mat1$pollution <- "Nitrate"
NO3_mat2$pollution <- "Nitrate"
NH4_mat1$pollution <- "Ammonium"
NH4_mat2$pollution <- "Ammonium"
BC_mat1$pollution <- "BC"
BC_mat2$pollution <- "BC"
OM_mat1$pollution <- "OM"
OM_mat2$pollution <- "OM"

PM2.5_mat3$pollution <- "PM2.5"
PM2.5_mat4$pollution <- "PM2.5"
PM2.5_mat5$pollution <- "PM2.5"
PM2.5_mat6$pollution <- "PM2.5"

SO4_mat3$pollution <- "Sulfate"
SO4_mat4$pollution <- "Sulfate"
SO4_mat5$pollution <- "Sulfate"
SO4_mat6$pollution <- "Sulfate"

NO3_mat3$pollution <- "Nitrate"
NO3_mat4$pollution <- "Nitrate"
NO3_mat5$pollution <- "Nitrate"
NO3_mat6$pollution <- "Nitrate"

NH4_mat3$pollution <- "Ammonium"
NH4_mat4$pollution <- "Ammonium"
NH4_mat5$pollution <- "Ammonium"
NH4_mat6$pollution <- "Ammonium"

BC_mat3$pollution <- "BC"
BC_mat4$pollution <- "BC"
BC_mat5$pollution <- "BC"
BC_mat6$pollution <- "BC"

OM_mat3$pollution <- "OM"
OM_mat4$pollution <- "OM"
OM_mat5$pollution <- "OM"
OM_mat6$pollution <- "OM"

PM2.5_mat1$var50 <- "40"
PM2.5_mat2$var50 <- "40"
SO4_mat1$var50 <- "8.5"
SO4_mat2$var50 <- "8.5"
NO3_mat1$var50 <- "9"
NO3_mat2$var50 <- "9"
NH4_mat1$var50 <- "7"
NH4_mat2$var50 <- "7"
BC_mat1$var50 <- "2"
BC_mat2$var50 <- "2"
OM_mat1$var50 <- "10"
OM_mat2$var50 <- "10"
PM2.5_mat3$var50 <- "40"
PM2.5_mat4$var50 <- "40"
PM2.5_mat5$var50 <- "40"
PM2.5_mat6$var50 <- "40"

SO4_mat3$var50 <- "8.5"
SO4_mat4$var50 <- "8.5"
SO4_mat5$var50 <- "8.5"
SO4_mat6$var50 <- "8.5"

NO3_mat3$var50 <- "9"
NO3_mat4$var50 <- "9"
NO3_mat5$var50 <- "9"
NO3_mat6$var50 <- "9"

NH4_mat3$var50 <- "7"
NH4_mat4$var50 <- "7"
NH4_mat5$var50 <- "7"
NH4_mat6$var50 <- "7"

BC_mat3$var50 <- "2"
BC_mat4$var50 <- "2"
BC_mat5$var50 <- "2"
BC_mat6$var50 <- "2"

OM_mat3$var50 <- "10"
OM_mat4$var50 <- "10"
OM_mat5$var50 <- "10"
OM_mat6$var50 <- "10"


PM2.5_mat1$var75 <- "60"
PM2.5_mat2$var75 <- "60"
SO4_mat1$var75 <- "12"
SO4_mat2$var75 <- "12"
NO3_mat1$var75 <- "16"
NO3_mat2$var75 <- "16"
NH4_mat1$var75 <- "10.5"
NH4_mat2$var75 <- "10.5"
BC_mat1$var75 <- "2.8"
BC_mat2$var75 <- "2.8"
OM_mat1$var75 <- "14"
OM_mat2$var75 <- "14"
PM2.5_mat3$var75 <- "60"
PM2.5_mat4$var75 <- "60"
PM2.5_mat5$var75 <- "60"
PM2.5_mat6$var75 <- "60"

SO4_mat3$var75 <- "12"
SO4_mat4$var75 <- "12"
SO4_mat5$var75 <- "12"
SO4_mat6$var75 <- "12"

NO3_mat3$var75 <- "16"
NO3_mat4$var75 <- "16"
NO3_mat5$var75 <- "16"
NO3_mat6$var75 <- "16"

NH4_mat3$var75 <- "10.5"
NH4_mat4$var75 <- "10.5"
NH4_mat5$var75 <- "10.5"
NH4_mat6$var75 <- "10.5"

BC_mat3$var75 <- "2.8"
BC_mat4$var75 <- "2.8"
BC_mat5$var75 <- "2.8"
BC_mat6$var75 <- "2.8"

OM_mat3$var75 <- "14"
OM_mat4$var75 <- "14"
OM_mat5$var75 <- "14"
OM_mat6$var75 <- "14"


PM2.5_mat1$var95 <- "95"
PM2.5_mat2$var95 <- "95"
SO4_mat1$var95 <- "17"
SO4_mat2$var95 <- "17"
NO3_mat1$var95 <- "25"
NO3_mat2$var95 <- "25"
NH4_mat1$var95 <- "16"
NH4_mat2$var95 <- "16"
BC_mat1$var95 <- "4.8"
BC_mat2$var95 <- "4.8"
OM_mat1$var95 <- "24"
OM_mat2$var95 <- "24"
PM2.5_mat3$var95 <- "95"
PM2.5_mat4$var95 <- "95"
PM2.5_mat5$var95 <- "95"
PM2.5_mat6$var95 <- "95"

SO4_mat3$var95 <- "17"
SO4_mat4$var95 <- "17"
SO4_mat5$var95 <- "17"
SO4_mat6$var95 <- "17"

NO3_mat3$var95 <- "25"
NO3_mat4$var95 <- "25"
NO3_mat5$var95 <- "25"
NO3_mat6$var95 <- "25"

NH4_mat3$var95 <- "16"
NH4_mat4$var95 <- "16"
NH4_mat5$var95 <- "16"
NH4_mat6$var95 <- "16"

BC_mat3$var95 <- "4.8"
BC_mat4$var95 <- "4.8"
BC_mat5$var95 <- "4.8"
BC_mat6$var95 <- "4.8"

OM_mat3$var95 <- "24"
OM_mat4$var95 <- "24"
OM_mat5$var95 <- "24"
OM_mat6$var95 <- "24"

sex_contour <- bind_rows(PM2.5_mat1, PM2.5_mat2, 
                         SO4_mat1, SO4_mat2, 
                         NO3_mat1, NO3_mat2, 
                         NH4_mat1, NH4_mat2,
                         BC_mat1, BC_mat2, 
                         OM_mat1, OM_mat2)

age_contour <- bind_rows(PM2.5_mat3, PM2.5_mat4, PM2.5_mat5, PM2.5_mat6, 
                         SO4_mat3, SO4_mat4, SO4_mat5, SO4_mat6, 
                         NO3_mat3, NO3_mat4, NO3_mat5, NO3_mat6,
                         NH4_mat3, NH4_mat4, NH4_mat5, NH4_mat6, 
                         BC_mat3, BC_mat4, BC_mat5, BC_mat6,
                         OM_mat3, OM_mat4, OM_mat5, OM_mat6)

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

percentiles <- c(0.05, 0.15,  0.25, 0.35,  0.45, 0.55,  0.65, 0.75, 0.85, 0.95)

specific_concentrations_PM2.5 <- quantile(data$PM2.5, percentiles, na.rm = TRUE)
specific_concentrations_SO4 <- quantile(data$SO4, percentiles, na.rm = TRUE)
specific_concentrations_NO3 <- quantile(data$NO3, percentiles, na.rm = TRUE)
specific_concentrations_NH4 <- quantile(data$NH4, percentiles, na.rm = TRUE)
specific_concentrations_BC <- quantile(data$BC, percentiles, na.rm = TRUE)
specific_concentrations_OM <- quantile(data$OM, percentiles, na.rm = TRUE)

#================== PM2.5 =====================
PM2.5_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

PM2.5_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT

##man
PM2.5_model1 <- fitmodel(PM2.5_model1, data) 

coef <- PM2.5_model1$summary.fixed$mean

vcov <- PM2.5_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model1$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

##woman
PM2.5_model2 <- fitmodel(PM2.5_model2, data) 

coef <- PM2.5_model2$summary.fixed$mean

vcov <- PM2.5_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model2$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

##[0, 25)
PM2.5_model3 <- fitmodel(PM2.5_model3, data) 

coef <- PM2.5_model3$summary.fixed$mean

vcov <- PM2.5_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model3$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum3 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

##[25, 45)
PM2.5_model4 <- fitmodel(PM2.5_model4, data) 

coef <- PM2.5_model4$summary.fixed$mean

vcov <- PM2.5_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model4$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum4 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

##[45, 65)
PM2.5_model5 <- fitmodel(PM2.5_model5, data) 

coef <- PM2.5_model5$summary.fixed$mean

vcov <- PM2.5_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model5$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum5 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

##[65, ~)
PM2.5_model6 <- fitmodel(PM2.5_model6, data) 

coef <- PM2.5_model6$summary.fixed$mean

vcov <- PM2.5_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", PM2.5_model6$names.fixed)

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum6 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
SO4_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

SO4_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT

##man
SO4_model1 <- fitmodel(SO4_model1, data) 

coef <- SO4_model1$summary.fixed$mean

vcov <- SO4_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model1$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

##woman
SO4_model2 <- fitmodel(SO4_model2, data) 

coef <- SO4_model2$summary.fixed$mean

vcov <- SO4_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model2$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

##[0, 25)
SO4_model3 <- fitmodel(SO4_model3, data) 

coef <- SO4_model3$summary.fixed$mean

vcov <- SO4_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model3$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum3 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

##[25, 45)
SO4_model4 <- fitmodel(SO4_model4, data) 

coef <- SO4_model4$summary.fixed$mean

vcov <- SO4_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model4$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum4 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

##[45, 65)
SO4_model5 <- fitmodel(SO4_model5, data) 

coef <- SO4_model5$summary.fixed$mean

vcov <- SO4_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model5$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum5 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

##[65, ~)
SO4_model6 <- fitmodel(SO4_model6, data) 

coef <- SO4_model6$summary.fixed$mean

vcov <- SO4_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", SO4_model6$names.fixed)

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum6 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
NO3_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

NO3_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT

##man
NO3_model1 <- fitmodel(NO3_model1, data) 

coef <- NO3_model1$summary.fixed$mean

vcov <- NO3_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model1$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

##woman
NO3_model2 <- fitmodel(NO3_model2, data) 

coef <- NO3_model2$summary.fixed$mean

vcov <- NO3_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model2$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

##[0, 25)
NO3_model3 <- fitmodel(NO3_model3, data) 

coef <- NO3_model3$summary.fixed$mean

vcov <- NO3_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model3$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum3 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

##[25, 45)
NO3_model4 <- fitmodel(NO3_model4, data) 

coef <- NO3_model4$summary.fixed$mean

vcov <- NO3_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model4$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum4 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

##[45, 65)
NO3_model5 <- fitmodel(NO3_model5, data) 

coef <- NO3_model5$summary.fixed$mean

vcov <- NO3_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model5$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum5 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

##[65, ~)
NO3_model6 <- fitmodel(NO3_model6, data) 

coef <- NO3_model6$summary.fixed$mean

vcov <- NO3_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", NO3_model6$names.fixed)

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum6 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
NH4_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

NH4_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT

##man
NH4_model1 <- fitmodel(NH4_model1, data) 

coef <- NH4_model1$summary.fixed$mean

vcov <- NH4_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model1$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

##woman
NH4_model2 <- fitmodel(NH4_model2, data) 

coef <- NH4_model2$summary.fixed$mean

vcov <- NH4_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model2$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

##[0, 25)
NH4_model3 <- fitmodel(NH4_model3, data) 

coef <- NH4_model3$summary.fixed$mean

vcov <- NH4_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model3$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum3 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

##[25, 45)
NH4_model4 <- fitmodel(NH4_model4, data) 

coef <- NH4_model4$summary.fixed$mean

vcov <- NH4_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model4$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum4 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

##[45, 65)
NH4_model5 <- fitmodel(NH4_model5, data) 

coef <- NH4_model5$summary.fixed$mean

vcov <- NH4_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model5$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum5 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

##[65, ~)
NH4_model6 <- fitmodel(NH4_model6, data) 

coef <- NH4_model6$summary.fixed$mean

vcov <- NH4_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", NH4_model6$names.fixed)

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum6 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== OM =====================
OM_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

OM_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT

##man
OM_model1 <- fitmodel(OM_model1, data) 

coef <- OM_model1$summary.fixed$mean

vcov <- OM_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model1$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

##woman
OM_model2 <- fitmodel(OM_model2, data) 

coef <- OM_model2$summary.fixed$mean

vcov <- OM_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model2$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

##[0, 25)
OM_model3 <- fitmodel(OM_model3, data) 

coef <- OM_model3$summary.fixed$mean

vcov <- OM_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model3$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum3 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

##[25, 45)
OM_model4 <- fitmodel(OM_model4, data) 

coef <- OM_model4$summary.fixed$mean

vcov <- OM_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model4$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum4 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

##[45, 65)
OM_model5 <- fitmodel(OM_model5, data) 

coef <- OM_model5$summary.fixed$mean

vcov <- OM_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model5$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum5 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

##[65, ~)
OM_model6 <- fitmodel(OM_model6, data) 

coef <- OM_model6$summary.fixed$mean

vcov <- OM_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", OM_model6$names.fixed)

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum6 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== BC =====================
BC_model1 <- man ~ offset(log(E1)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model2 <- woman ~ offset(log(E2)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model3 <- age1 ~ offset(log(E3)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model4 <- age2 ~ offset(log(E4)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model5 <- age3 ~ offset(log(E5)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

BC_model6 <- age4 ~ offset(log(E6)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT

##man
BC_model1 <- fitmodel(BC_model1, data) 

coef <- BC_model1$summary.fixed$mean

vcov <- BC_model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model1$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

##woman
BC_model2 <- fitmodel(BC_model2, data) 

coef <- BC_model2$summary.fixed$mean

vcov <- BC_model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model2$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

##[0, 25)
BC_model3 <- fitmodel(BC_model3, data) 

coef <- BC_model3$summary.fixed$mean

vcov <- BC_model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model3$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum3 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

##[25, 45)
BC_model4 <- fitmodel(BC_model4, data) 

coef <- BC_model4$summary.fixed$mean

vcov <- BC_model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model4$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum4 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

##[45, 65)
BC_model5 <- fitmodel(BC_model5, data) 

coef <- BC_model5$summary.fixed$mean

vcov <- BC_model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model5$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum5 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

##[65, ~)
BC_model6 <- fitmodel(BC_model6, data) 

coef <- BC_model6$summary.fixed$mean

vcov <- BC_model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", BC_model6$names.fixed)

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum6 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$group <- "Male"
PM2.5_cum2$group <- "Female"
SO4_cum1$group <- "Male"
SO4_cum2$group <- "Female"
NO3_cum1$group <- "Male"
NO3_cum2$group <- "Female"
NH4_cum1$group <- "Male"
NH4_cum2$group <- "Female"
BC_cum1$group <- "Male"
BC_cum2$group <- "Female"
OM_cum1$group <- "Male"
OM_cum2$group <- "Female"

PM2.5_cum3$group <- "[0, 25)"
PM2.5_cum4$group <- "[25, 45)"
PM2.5_cum5$group <- "[45, 65)"
PM2.5_cum6$group <- "[65, ~)"

SO4_cum3$group <- "[0, 25)"
SO4_cum4$group <- "[25, 45)"
SO4_cum5$group <- "[45, 65)"
SO4_cum6$group <- "[65, ~)"

NO3_cum3$group <- "[0, 25)"
NO3_cum4$group <- "[25, 45)"
NO3_cum5$group <- "[45, 65)"
NO3_cum6$group <- "[65, ~)"

NH4_cum3$group <- "[0, 25)"
NH4_cum4$group <- "[25, 45)"
NH4_cum5$group <- "[45, 65)"
NH4_cum6$group <- "[65, ~)"

BC_cum3$group <- "[0, 25)"
BC_cum4$group <- "[25, 45)"
BC_cum5$group <- "[45, 65)"
BC_cum6$group <- "[65, ~)"

OM_cum3$group <- "[0, 25)"
OM_cum4$group <- "[25, 45)"
OM_cum5$group <- "[45, 65)"
OM_cum6$group <- "[65, ~)"

PM2.5_cum1$pollution <- "PM2.5"
PM2.5_cum2$pollution <- "PM2.5"
SO4_cum1$pollution <- "Sulfate"
SO4_cum2$pollution <- "Sulfate"
NO3_cum1$pollution <- "Nitrate"
NO3_cum2$pollution <- "Nitrate"
NH4_cum1$pollution <- "Ammonium"
NH4_cum2$pollution <- "Ammonium"
BC_cum1$pollution <- "BC"
BC_cum2$pollution <- "BC"
OM_cum1$pollution <- "OM"
OM_cum2$pollution <- "OM"
PM2.5_cum3$pollution <- "PM2.5"
PM2.5_cum4$pollution <- "PM2.5"
PM2.5_cum5$pollution <- "PM2.5"
PM2.5_cum6$pollution <- "PM2.5"

SO4_cum3$pollution <- "Sulfate"
SO4_cum4$pollution <- "Sulfate"
SO4_cum5$pollution <- "Sulfate"
SO4_cum6$pollution <- "Sulfate"

NO3_cum3$pollution <- "Nitrate"
NO3_cum4$pollution <- "Nitrate"
NO3_cum5$pollution <- "Nitrate"
NO3_cum6$pollution <- "Nitrate"

NH4_cum3$pollution <- "Ammonium"
NH4_cum4$pollution <- "Ammonium"
NH4_cum5$pollution <- "Ammonium"
NH4_cum6$pollution <- "Ammonium"

BC_cum3$pollution <- "BC"
BC_cum4$pollution <- "BC"
BC_cum5$pollution <- "BC"
BC_cum6$pollution <- "BC"

OM_cum3$pollution <- "OM"
OM_cum4$pollution <- "OM"
OM_cum5$pollution <- "OM"
OM_cum6$pollution <- "OM"


sex_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, 
                            SO4_cum1, SO4_cum2, 
                            NO3_cum1, NO3_cum2, 
                            NH4_cum1, NH4_cum2,
                            BC_cum1, BC_cum2, 
                            OM_cum1, OM_cum2)

age_cumcontour <- bind_rows(PM2.5_cum3, PM2.5_cum4, PM2.5_cum5, PM2.5_cum6,
                            SO4_cum3, SO4_cum4, SO4_cum5, SO4_cum6,
                            NO3_cum3, NO3_cum4, NO3_cum5, NO3_cum6, 
                            NH4_cum3, NH4_cum4, NH4_cum5, NH4_cum6,
                            BC_cum3, BC_cum4, BC_cum5, BC_cum6,  
                            OM_cum3, OM_cum4, OM_cum5, OM_cum6)

#================== Figures export =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM) 

AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM)

data_pollution <- rbind(AH.data, SC.data)

#================== Sex contour =====================
sex_contour <- sex_contour %>%
  mutate(group = factor(group, levels = c("Male", "Female")))

PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- sex_contour   %>%
  filter(pollution == "PM2.5") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.65, 1.2)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

SO4_contour <- sex_contour   %>%
  filter(pollution == "Sulfate") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.65, 1.2)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

NO3_contour <- sex_contour   %>%
  filter(pollution == "Nitrate") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.65, 1.2)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

NH4_contour <- sex_contour   %>%
  filter(pollution == "Ammonium") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.65, 1.2)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

OM_contour <- sex_contour   %>%
  filter(pollution == "OM") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.65, 1.2)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

BC_contour <- sex_contour   %>%
  filter(pollution == "BC") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.65, 1.2)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

pdf("./Figures/Figure S8.pdf", height = 9, width = 12)
All_contour 
dev.off()

#================== Sex lag =====================
sex_contour <- sex_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) 
  ))

sex_contour <- sex_contour %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

sex_contour <- sex_contour %>%
  mutate(group = factor(group, levels = c("Male", "Female")))

col.pal <- c("#2171b5", "#a50f15")

Lag_sex <- sex_contour %>%
  group_by(group, pollution, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.7, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ pollution, scales = "fixed") +
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

Lag_sex <- ggdraw(Lag_sex) +
  annotate("text", x = 0.13, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.44, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure S9.pdf", height = 7, width = 14)
Lag_sex 
dev.off()

#================== Age contour =====================
PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- age_contour   %>%
  filter(pollution == "PM2.5") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.60, 1.45)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

SO4_contour <- age_contour   %>%
  filter(pollution == "Sulfate") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.60, 1.45)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

NO3_contour <- age_contour   %>%
  filter(pollution == "Nitrate") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.60, 1.45)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

NH4_contour <- age_contour   %>%
  filter(pollution == "Ammonium") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.60, 1.45)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

OM_contour <- age_contour   %>%
  filter(pollution == "OM") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.60, 1.45)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

BC_contour <- age_contour   %>%
  filter(pollution == "BC") %>%
  group_by(group) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = group) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.60, 1.45)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
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

pdf("./Figures/Figure S11.pdf", height = 8, width = 14)
All_contour 
dev.off()

#================== Age lag =====================
age_contour <- age_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) # 其他值保持不变
  ))

age_contour <- age_contour %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#238b45","#2171b5","#810f7c", "#a50f15")

Lag_age <- age_contour %>%
  group_by(group, pollution, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.4, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ pollution, scales = "fixed") +
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

Lag_age <- ggdraw(Lag_age) +
  annotate("text", x = 0.125, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.435, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure S12.pdf", height = 7, width = 14)
Lag_age 
dev.off()

#================== Sex exposure =====================
sex_cumcontour <- sex_cumcontour %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

sex_cumcontour <- sex_cumcontour %>%
  mutate(group = factor(group, levels = c("Male", "Female")))

col.pal <- c("#EE0000FF", "#3B4992FF")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- sex_cumcontour %>%
  filter(pollution == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.75, 1.5)) +
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

SO4 <- sex_cumcontour %>%
  filter(pollution == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.75, 1.5)) +
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

NO3 <- sex_cumcontour %>%
  filter(pollution == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.75, 1.5)) +
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

NH4 <- sex_cumcontour %>%
  filter(pollution == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.75, 1.5)) +
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

OM <- sex_cumcontour %>%
  filter(pollution == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.75, 1.5)) +
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

BC <- sex_cumcontour %>%
  filter(pollution == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.75, 1.5)) +
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
          legend.text = element_text(size = 15),
          legend.spacing.y = unit(0.2, "cm")) + 
    guides(color = guide_legend(ncol = 1)))

exposure_sex <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure_sex <- ggdraw(exposure_sex) +
  annotate("text", x = 0.142, y = 0.975, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 6, hjust = 0.5) +
  annotate("text", x = 0.482, y = 0.975, label = "B.Sulfate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.814, y = 0.975, label = "C.Nitrate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.171, y = 0.475, label = "D.Ammonium (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.466, y = 0.475, label = "E.OM (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.795, y = 0.475, label = "F.BC (μg/m³)", size = 6, hjust = 0.5)

exposure_sex <- ggdraw() +
  draw_plot(exposure_sex, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.91,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S7.pdf", height = 7, width = 12)
exposure_sex
dev.off()

#================== Age exposure =====================
age_cumcontour <- age_cumcontour %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#238b45","#2171b5","#810f7c","#a50f15")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- age_cumcontour %>%
  filter(pollution == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.7, 2)) +
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

SO4 <- age_cumcontour %>%
  filter(pollution == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.7, 2)) +
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

NO3 <- age_cumcontour %>%
  filter(pollution == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.7, 2)) +
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

NH4 <- age_cumcontour %>%
  filter(pollution == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.7, 2)) +
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

OM <- age_cumcontour %>%
  filter(pollution == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.7, 2)) +
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

BC <- age_cumcontour %>%
  filter(pollution == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~pollution, scales = "free") +
  scale_y_continuous(limits = c(0.7, 2)) +
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
          legend.text = element_text(size = 13),
          legend.spacing.y = unit(0.2, "cm")) + 
    guides(color = guide_legend(ncol = 1)))

exposure_age <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure_age <- ggdraw(exposure_age) +
  annotate("text", x = 0.14, y = 0.975, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 6, hjust = 0.5) +
  annotate("text", x = 0.481, y = 0.975, label = "B.Sulfate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.813, y = 0.975, label = "C.Nitrate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.168, y = 0.475, label = "D.Ammonium (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.465, y = 0.475, label = "E.OM (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.795, y = 0.475, label = "F.BC (μg/m³)", size = 6, hjust = 0.5)

exposure_age <- ggdraw() +
  draw_plot(exposure_age, 
            x = 0, 
            y = 0, 
            width = 0.93, 
            height = 1) +
  draw_plot(legend, 
            x = 0.923,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S10.pdf", height = 7, width = 12)
exposure_age
dev.off()