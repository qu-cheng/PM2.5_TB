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

#================== Incidence rate =====================

#================== Data import =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  mutate(Pref = Pref + 16)%>%
  dplyr::select(year, month, X, age1, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop3, regnames, MT) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(year, month, X, age1, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop3, regnames, MT)%>%
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
         E = pop3*10000)

#================== Calculate the annual average incidence rate =====================
annual_cases <- data %>%
  group_by(year, Pref) %>%
  summarise(total_cases = sum(age1, na.rm = TRUE))

SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  mutate(Pref = Pref + 16)%>%
  dplyr::select(year, Pref, pop3) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(year, Pref, pop3)%>%
  mutate(prov = "AH")

data_pop <- rbind(AH.data, SC.data) %>%
  mutate(E = pop3*10000)

unique_data <- data_pop %>%
  distinct()

merged_data <- left_join(annual_cases, unique_data, by = c("year", "Pref"))

merged_data$IR <- (merged_data$total_cases/merged_data$E)*100000

annual_IR <- merged_data %>%
  group_by(Pref) %>%
  summarise(annual_IR = mean(IR, na.rm = TRUE))

data <- left_join(data, annual_IR, by = c("Pref"))

data <- data %>%
  rename(IR = annual_IR)

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

#======================== Linear interaction term ==============================
IR_ind1 <- data$IR - quantile(data$IR, 0.90)
IR_ind2 <- data$IR - quantile(data$IR, 0.5)
IR_ind3 <- data$IR - quantile(data$IR, 0.10)

#================== PM2.5 =====================
IR_basis1_PM2.5 <- basis_PM2.5*IR_ind1
IR_basis2_PM2.5 <- basis_PM2.5*IR_ind2
IR_basis3_PM2.5 <- basis_PM2.5*IR_ind3

colnames(IR_basis1_PM2.5) = paste0("IR_basis1_PM2.5.", colnames(IR_basis1_PM2.5))
colnames(IR_basis2_PM2.5) = paste0("IR_basis2_PM2.5.", colnames(IR_basis2_PM2.5))
colnames(IR_basis3_PM2.5) = paste0("IR_basis3_PM2.5.", colnames(IR_basis3_PM2.5))

model1 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  IR_basis1_PM2.5 +
  IR

model2 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  IR_basis2_PM2.5 +
  IR

model3 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  IR_basis3_PM2.5 +
  IR

model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

model2 <- fitmodel(model2, data) 

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

model3 <- fitmodel(model3, data) 

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat3 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

#================== SO4 =====================
IR_basis1_SO4 <- basis_SO4*IR_ind1
IR_basis2_SO4 <- basis_SO4*IR_ind2
IR_basis3_SO4 <- basis_SO4*IR_ind3

colnames(IR_basis1_SO4) = paste0("IR_basis1_SO4.", colnames(IR_basis1_SO4))
colnames(IR_basis2_SO4) = paste0("IR_basis2_SO4.", colnames(IR_basis2_SO4))
colnames(IR_basis3_SO4) = paste0("IR_basis3_SO4.", colnames(IR_basis3_SO4))

model4 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  IR_basis1_SO4 +
  IR

model5 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  IR_basis2_SO4 +
  IR

model6 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT + 
  IR_basis3_SO4 +
  IR

model4 <- fitmodel(model4, data) 

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model4$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

model5 <- fitmodel(model5, data) 

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model5$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

model6 <- fitmodel(model6, data) 

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model6$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat3 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
IR_basis1_NO3 <- basis_NO3*IR_ind1
IR_basis2_NO3 <- basis_NO3*IR_ind2
IR_basis3_NO3 <- basis_NO3*IR_ind3

colnames(IR_basis1_NO3) = paste0("IR_basis1_NO3.", colnames(IR_basis1_NO3))
colnames(IR_basis2_NO3) = paste0("IR_basis2_NO3.", colnames(IR_basis2_NO3))
colnames(IR_basis3_NO3) = paste0("IR_basis3_NO3.", colnames(IR_basis3_NO3))

model7 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  IR_basis1_NO3 +
  IR

model8 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  IR_basis2_NO3 +
  IR

model9 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  IR_basis3_NO3 +
  IR

model7 <- fitmodel(model7, data) 

coef <- model7$summary.fixed$mean

vcov <- model7$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model7$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

model8 <- fitmodel(model8, data) 

coef <- model8$summary.fixed$mean

vcov <- model8$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model8$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

model9 <- fitmodel(model9, data) 

coef <- model9$summary.fixed$mean

vcov <- model9$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model9$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat3 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
IR_basis1_NH4 <- basis_NH4*IR_ind1
IR_basis2_NH4 <- basis_NH4*IR_ind2
IR_basis3_NH4 <- basis_NH4*IR_ind3

colnames(IR_basis1_NH4) = paste0("IR_basis1_NH4.", colnames(IR_basis1_NH4))
colnames(IR_basis2_NH4) = paste0("IR_basis2_NH4.", colnames(IR_basis2_NH4))
colnames(IR_basis3_NH4) = paste0("IR_basis3_NH4.", colnames(IR_basis3_NH4))

model10 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  IR_basis1_NH4 +
  IR

model11 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  IR_basis2_NH4 +
  IR

model12 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  IR_basis3_NH4 +
  IR

model10 <- fitmodel(model10, data) 

coef <- model10$summary.fixed$mean

vcov <- model10$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model10$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

model11 <- fitmodel(model11, data) 

coef <- model11$summary.fixed$mean

vcov <- model11$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model11$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

model12 <- fitmodel(model12, data) 

coef <- model12$summary.fixed$mean

vcov <- model12$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model12$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat3 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
IR_basis1_BC <- basis_BC*IR_ind1
IR_basis2_BC <- basis_BC*IR_ind2
IR_basis3_BC <- basis_BC*IR_ind3

colnames(IR_basis1_BC) = paste0("IR_basis1_BC.", colnames(IR_basis1_BC))
colnames(IR_basis2_BC) = paste0("IR_basis2_BC.", colnames(IR_basis2_BC))
colnames(IR_basis3_BC) = paste0("IR_basis3_BC.", colnames(IR_basis3_BC))

model13 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  IR_basis1_BC +
  IR

model14 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  IR_basis2_BC +
  IR

model15 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  IR_basis3_BC +
  IR

model13 <- fitmodel(model13, data) 

coef <- model13$summary.fixed$mean

vcov <- model13$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model13$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

model14 <- fitmodel(model14, data) 

coef <- model14$summary.fixed$mean

vcov <- model14$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model14$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

model15 <- fitmodel(model15, data) 

coef <- model15$summary.fixed$mean

vcov <- model15$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model15$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat3 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
IR_basis1_OM <- basis_OM*IR_ind1
IR_basis2_OM <- basis_OM*IR_ind2
IR_basis3_OM <- basis_OM*IR_ind3

colnames(IR_basis1_OM) = paste0("IR_basis1_OM.", colnames(IR_basis1_OM))
colnames(IR_basis2_OM) = paste0("IR_basis2_OM.", colnames(IR_basis2_OM))
colnames(IR_basis3_OM) = paste0("IR_basis3_OM.", colnames(IR_basis3_OM))

model16 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  IR_basis1_OM +
  IR

model17 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  IR_basis2_OM +
  IR

model18 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  IR_basis3_OM +
  IR

model16 <- fitmodel(model16, data) 

coef <- model16$summary.fixed$mean

vcov <- model16$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model16$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))


model17 <- fitmodel(model17, data) 

coef <- model17$summary.fixed$mean

vcov <- model17$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model17$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))


model18 <- fitmodel(model18, data) 

coef <- model18$summary.fixed$mean

vcov <- model18$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model18$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat3 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$IR <- "High IR"
PM2.5_mat2$IR <- "Middle IR"
PM2.5_mat3$IR <- "Low IR"
SO4_mat1$IR <- "High IR"
SO4_mat2$IR <- "Middle IR"
SO4_mat3$IR <- "Low IR"
NO3_mat1$IR <- "High IR"
NO3_mat2$IR <- "Middle IR"
NO3_mat3$IR <- "Low IR"
NH4_mat1$IR <- "High IR"
NH4_mat2$IR <- "Middle IR"
NH4_mat3$IR <- "Low IR"
BC_mat1$IR <- "High IR"
BC_mat2$IR <- "Middle IR"
BC_mat3$IR <- "Low IR"
OM_mat1$IR <- "High IR"
OM_mat2$IR <- "Middle IR"
OM_mat3$IR <- "Low IR"

PM2.5_mat1$Pollution <- "PM2.5"
PM2.5_mat2$Pollution <- "PM2.5"
PM2.5_mat3$Pollution <- "PM2.5"
SO4_mat1$Pollution <- "Sulfate"
SO4_mat2$Pollution <- "Sulfate"
SO4_mat3$Pollution <- "Sulfate"
NO3_mat1$Pollution <- "Nitrate"
NO3_mat2$Pollution <- "Nitrate"
NO3_mat3$Pollution <- "Nitrate"
NH4_mat1$Pollution <- "Ammonium"
NH4_mat2$Pollution <- "Ammonium"
NH4_mat3$Pollution <- "Ammonium"
BC_mat1$Pollution <- "BC"
BC_mat2$Pollution <- "BC"
BC_mat3$Pollution <- "BC"
OM_mat1$Pollution <- "OM"
OM_mat2$Pollution <- "OM"
OM_mat3$Pollution <- "OM"

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
SO4_mat3$var50 <- "8.5"
NO3_mat3$var50 <- "9"
NH4_mat3$var50 <- "7"
BC_mat3$var50 <- "2"
OM_mat3$var50 <- "10"


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
SO4_mat3$var75 <- "12"
NO3_mat3$var75 <- "16"
NH4_mat3$var75 <- "10.5"
BC_mat3$var75 <- "2.8"
OM_mat3$var75 <- "14"


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
SO4_mat3$var95 <- "17"
NO3_mat3$var95 <- "25"
NH4_mat3$var95 <- "16"
BC_mat3$var95 <- "4.8"
OM_mat3$var95 <- "24"


IR_contour <- bind_rows(PM2.5_mat1, PM2.5_mat2, PM2.5_mat3, SO4_mat1, SO4_mat2, SO4_mat3, 
                NO3_mat1, NO3_mat2, NO3_mat3, NH4_mat1, NH4_mat2, NH4_mat3, 
                BC_mat1, BC_mat2, BC_mat3, OM_mat1, OM_mat2, OM_mat3)

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

#======================== Linear interaction term ==============================
IR_ind1 <- data$IR - quantile(data$IR, p = 0.90) 
IR_ind2 <- data$IR - quantile(data$IR, p = 0.5) 
IR_ind3 <- data$IR - quantile(data$IR, p = 0.10) 

#================== PM2.5 =====================
IR_basis1_PM2.5 <- basis_PM2.5*IR_ind1
IR_basis2_PM2.5 <- basis_PM2.5*IR_ind2
IR_basis3_PM2.5 <- basis_PM2.5*IR_ind3

colnames(IR_basis1_PM2.5) = paste0("IR_basis1_PM2.5.", colnames(IR_basis1_PM2.5))
colnames(IR_basis2_PM2.5) = paste0("IR_basis2_PM2.5.", colnames(IR_basis2_PM2.5))
colnames(IR_basis3_PM2.5) = paste0("IR_basis3_PM2.5.", colnames(IR_basis3_PM2.5))

model1 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  IR_basis1_PM2.5 +
  IR

model2 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  IR_basis2_PM2.5 +
  IR

model3 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  IR_basis3_PM2.5 +
  IR

model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

model2 <- fitmodel(model2, data) 

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

model3 <- fitmodel(model3, data) 

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum3 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
IR_basis1_SO4 <- basis_SO4*IR_ind1
IR_basis2_SO4 <- basis_SO4*IR_ind2
IR_basis3_SO4 <- basis_SO4*IR_ind3

colnames(IR_basis1_SO4) = paste0("IR_basis1_SO4.", colnames(IR_basis1_SO4))
colnames(IR_basis2_SO4) = paste0("IR_basis2_SO4.", colnames(IR_basis2_SO4))
colnames(IR_basis3_SO4) = paste0("IR_basis3_SO4.", colnames(IR_basis3_SO4))

model4 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  IR_basis1_SO4 +
  IR

model5 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  IR_basis2_SO4 +
  IR

model6 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT + 
  IR_basis3_SO4 +
  IR

model4 <- fitmodel(model4, data) 

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model4$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

model5 <- fitmodel(model5, data) 

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model5$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

model6 <- fitmodel(model6, data) 

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model6$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum3 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
IR_basis1_NO3 <- basis_NO3*IR_ind1
IR_basis2_NO3 <- basis_NO3*IR_ind2
IR_basis3_NO3 <- basis_NO3*IR_ind3

colnames(IR_basis1_NO3) = paste0("IR_basis1_NO3.", colnames(IR_basis1_NO3))
colnames(IR_basis2_NO3) = paste0("IR_basis2_NO3.", colnames(IR_basis2_NO3))
colnames(IR_basis3_NO3) = paste0("IR_basis3_NO3.", colnames(IR_basis3_NO3))

model7 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  IR_basis1_NO3 +
  IR

model8 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  IR_basis2_NO3 +
  IR

model9 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  IR_basis3_NO3 +
  IR

model7 <- fitmodel(model7, data) 

coef <- model7$summary.fixed$mean

vcov <- model7$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model7$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

model8 <- fitmodel(model8, data) 

coef <- model8$summary.fixed$mean

vcov <- model8$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model8$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

model9 <- fitmodel(model9, data) 

coef <- model9$summary.fixed$mean

vcov <- model9$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model9$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum3 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
IR_basis1_NH4 <- basis_NH4*IR_ind1
IR_basis2_NH4 <- basis_NH4*IR_ind2
IR_basis3_NH4 <- basis_NH4*IR_ind3

colnames(IR_basis1_NH4) = paste0("IR_basis1_NH4.", colnames(IR_basis1_NH4))
colnames(IR_basis2_NH4) = paste0("IR_basis2_NH4.", colnames(IR_basis2_NH4))
colnames(IR_basis3_NH4) = paste0("IR_basis3_NH4.", colnames(IR_basis3_NH4))

model10 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  IR_basis1_NH4 +
  IR

model11 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  IR_basis2_NH4 +
  IR

model12 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  IR_basis3_NH4 +
  IR

model10 <- fitmodel(model10, data) 

coef <- model10$summary.fixed$mean

vcov <- model10$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model10$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

model11 <- fitmodel(model11, data) 

coef <- model11$summary.fixed$mean

vcov <- model11$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model11$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

model12 <- fitmodel(model12, data) 

coef <- model12$summary.fixed$mean

vcov <- model12$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model12$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum3 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
IR_basis1_BC <- basis_BC*IR_ind1
IR_basis2_BC <- basis_BC*IR_ind2
IR_basis3_BC <- basis_BC*IR_ind3

colnames(IR_basis1_BC) = paste0("IR_basis1_BC.", colnames(IR_basis1_BC))
colnames(IR_basis2_BC) = paste0("IR_basis2_BC.", colnames(IR_basis2_BC))
colnames(IR_basis3_BC) = paste0("IR_basis3_BC.", colnames(IR_basis3_BC))

model13 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  IR_basis1_BC +
  IR

model14 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  IR_basis2_BC +
  IR

model15 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  IR_basis3_BC +
  IR

model13 <- fitmodel(model13, data) 

coef <- model13$summary.fixed$mean

vcov <- model13$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model13$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

model14 <- fitmodel(model14, data) 

coef <- model14$summary.fixed$mean

vcov <- model14$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model14$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

model15 <- fitmodel(model15, data) 

coef <- model15$summary.fixed$mean

vcov <- model15$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model15$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum3 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
IR_basis1_OM <- basis_OM*IR_ind1
IR_basis2_OM <- basis_OM*IR_ind2
IR_basis3_OM <- basis_OM*IR_ind3

colnames(IR_basis1_OM) = paste0("IR_basis1_OM.", colnames(IR_basis1_OM))
colnames(IR_basis2_OM) = paste0("IR_basis2_OM.", colnames(IR_basis2_OM))
colnames(IR_basis3_OM) = paste0("IR_basis3_OM.", colnames(IR_basis3_OM))

model16 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  IR_basis1_OM +
  IR

model17 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  IR_basis2_OM +
  IR

model18 <- age1 ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  IR_basis3_OM +
  IR

model16 <- fitmodel(model16, data) 

coef <- model16$summary.fixed$mean

vcov <- model16$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model16$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

model17 <- fitmodel(model17, data) 

coef <- model17$summary.fixed$mean

vcov <- model17$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model17$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

model18 <- fitmodel(model18, data) 

coef <- model18$summary.fixed$mean

vcov <- model18$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model18$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum3 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$IR <- "High IR"
PM2.5_cum2$IR <- "Middle IR"
PM2.5_cum3$IR <- "Low IR"
SO4_cum1$IR <- "High IR"
SO4_cum2$IR <- "Middle IR"
SO4_cum3$IR <- "Low IR"
NO3_cum1$IR <- "High IR"
NO3_cum2$IR <- "Middle IR"
NO3_cum3$IR <- "Low IR"
NH4_cum1$IR <- "High IR"
NH4_cum2$IR <- "Middle IR"
NH4_cum3$IR <- "Low IR"
BC_cum1$IR <- "High IR"
BC_cum2$IR <- "Middle IR"
BC_cum3$IR <- "Low IR"
OM_cum1$IR <- "High IR"
OM_cum2$IR <- "Middle IR"
OM_cum3$IR<- "Low IR"

PM2.5_cum1$pollution <- "PM2.5"
PM2.5_cum2$pollution <- "PM2.5"
PM2.5_cum3$pollution <- "PM2.5"
SO4_cum1$pollution <- "Sulfate"
SO4_cum2$pollution <- "Sulfate"
SO4_cum3$pollution <- "Sulfate"
NO3_cum1$pollution <- "Nitrate"
NO3_cum2$pollution <- "Nitrate"
NO3_cum3$pollution <- "Nitrate"
NH4_cum1$pollution <- "Ammonium"
NH4_cum2$pollution <- "Ammonium"
NH4_cum3$pollution <- "Ammonium"
BC_cum1$pollution <- "BC"
BC_cum2$pollution <- "BC"
BC_cum3$pollution <- "BC"
OM_cum1$pollution <- "OM"
OM_cum2$pollution <- "OM"
OM_cum3$pollution <- "OM"


IR_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, PM2.5_cum3, SO4_cum1, SO4_cum2, SO4_cum3, 
                NO3_cum1, NO3_cum2, NO3_cum3, NH4_cum1, NH4_cum2, NH4_cum3, 
                BC_cum1, BC_cum2, BC_cum3, OM_cum1, OM_cum2, OM_cum3)

#================== Figures export =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM) 

AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM)

data_pollution <- rbind(AH.data, SC.data)

#================== IR contour =====================
IR_contour <- IR_contour %>%
  mutate(IR = factor(IR, levels = c("Low IR", "Middle IR", "High IR")))

PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- IR_contour   %>%
  filter(Pollution == "PM2.5") %>%
  group_by(IR) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = IR) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.59, 2.11)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ IR, scales = "fixed", nrow = 1) +
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

SO4_contour <- IR_contour   %>%
  filter(Pollution == "Sulfate") %>%
  group_by(IR) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = IR) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.59, 2.11)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ IR, scales = "fixed", nrow = 1) +
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

NO3_contour <- IR_contour   %>%
  filter(Pollution == "Nitrate") %>%
  group_by(IR) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = IR) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.59, 2.11)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ IR, scales = "fixed", nrow = 1) +
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

NH4_contour <- IR_contour   %>%
  filter(Pollution == "Ammonium") %>%
  group_by(IR) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = IR) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.59, 2.11)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ IR, scales = "fixed", nrow = 1) +
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

OM_contour <- IR_contour   %>%
  filter(Pollution == "OM") %>%
  group_by(IR) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = IR) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.59, 2.11)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ IR, scales = "fixed", nrow = 1) +
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

BC_contour <- IR_contour   %>%
  filter(Pollution == "BC") %>%
  group_by(IR) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = IR) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.59, 2.11)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ IR, scales = "fixed", nrow = 1) +
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

pdf("./Figures/Figure S14.pdf", height = 8, width = 14)
All_contour 
dev.off()

#================== IR lag =====================
IR_contour <- IR_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue)
  ))

IR_contour <- IR_contour %>%
  mutate(pollution = factor(Pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

IR_contour <- IR_contour %>%
  mutate(IR = factor(IR, levels = c("Low IR", "Middle IR", "High IR")))

col.pal <- c("#238b45", "#2171b5", "#a50f15")

Lag_IR <- IR_contour %>%
  group_by(IR, pollution, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ Pollution, scales = "fixed") +
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

Lag_IR <- ggdraw(Lag_IR) +
  annotate("text", x = 0.125, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.435, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure S15.pdf", height = 7, width = 14)
Lag_IR 
dev.off()

#================== IR exposure =====================
IR_cumcontour <- IR_cumcontour %>%
  mutate(Pollution = factor(Pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

IR_cumcontour <- IR_cumcontour %>%
  mutate(IR = factor(IR, levels = c("Low IR", "Middle IR", "High IR")))

col.pal <- c( "#008B45FF", "#3B4992FF","#EE0000FF")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- data %>%
  filter(Pollution == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(IR, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.5, 2.2)) +
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

SO4 <- data %>%
  filter(Pollution == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(IR, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.5, 2.2)) +
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

NO3 <- data %>%
  filter(Pollution == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(IR, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.5, 2.2)) +
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

NH4 <- data %>%
  filter(Pollution == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(IR, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.5, 2.2)) +
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

OM <- data %>%
  filter(Pollution == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(IR, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.5, 2.2)) +
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

BC <- data %>%
  filter(Pollution == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(IR, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = IR, col = as.factor(IR), fill = as.factor(IR))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.5, 2.2)) +
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

exposure_IR <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure_IRe <- ggdraw(exposure_IR) +
  annotate("text", x = 0.124, y = 0.975, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 6, hjust = 0.5) +
  annotate("text", x = 0.463, y = 0.975, label = "B.Sulfate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.795, y = 0.975, label = "C.Nitrate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.148, y = 0.475, label = "D.Ammonium (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.449, y = 0.475, label = "E.OM (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.779, y = 0.475, label = "F.BC (μg/m³)", size = 6, hjust = 0.5)

exposure_IR <- ggdraw() +
  draw_plot(exposure_IR, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.92,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S13.pdf", height = 7, width = 14)
exposure_IR
dev.off()

#================== NDVI =====================

#================== Data import =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  mutate(Pref = Pref + 16)%>%
  dplyr::select(year, month, X, case, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop, regnames, MT, NDVI) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(year, month, X, case, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop, regnames, MT, NDVI)%>%
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

#======================== Linear interaction term ==============================
NDVI_ind1 <- data$NDVI - quantile(data$NDVI, p = 0.90) 
NDVI_ind2 <- data$NDVI - quantile(data$NDVI, p = 0.5) 
NDVI_ind3 <- data$NDVI - quantile(data$NDVI, p = 0.10) 

#================== PM2.5 =====================
NDVI_basis1_PM2.5 <- basis_PM2.5*NDVI_ind1
NDVI_basis2_PM2.5 <- basis_PM2.5*NDVI_ind2
NDVI_basis3_PM2.5 <- basis_PM2.5*NDVI_ind3

colnames(NDVI_basis1_PM2.5) = paste0("NDVI_basis1_PM2.5.", colnames(NDVI_basis1_PM2.5))
colnames(NDVI_basis2_PM2.5) = paste0("NDVI_basis2_PM2.5.", colnames(NDVI_basis2_PM2.5))
colnames(NDVI_basis3_PM2.5) = paste0("NDVI_basis3_PM2.5.", colnames(NDVI_basis3_PM2.5))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  NDVI_basis1_PM2.5 +
  NDVI

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  NDVI_basis2_PM2.5 +
  NDVI

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  NDVI_basis3_PM2.5 +
  NDVI

model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

model2 <- fitmodel(model2, data) 

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

model3 <- fitmodel(model3, data) 

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat3 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

#================== SO4 =====================
NDVI_basis1_SO4 <- basis_SO4*NDVI_ind1
NDVI_basis2_SO4 <- basis_SO4*NDVI_ind2
NDVI_basis3_SO4 <- basis_SO4*NDVI_ind3

colnames(NDVI_basis1_SO4) = paste0("NDVI_basis1_SO4.", colnames(NDVI_basis1_SO4))
colnames(NDVI_basis2_SO4) = paste0("NDVI_basis2_SO4.", colnames(NDVI_basis2_SO4))
colnames(NDVI_basis3_SO4) = paste0("NDVI_basis3_SO4.", colnames(NDVI_basis3_SO4))

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  NDVI_basis1_SO4 +
  NDVI

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  NDVI_basis2_SO4 +
  NDVI

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT + 
  NDVI_basis3_SO4 +
  NDVI

model4 <- fitmodel(model4, data) 

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model4$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

model5 <- fitmodel(model5, data) 

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model5$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

model6 <- fitmodel(model6, data) 

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model6$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat3 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
NDVI_basis1_NO3 <- basis_NO3*NDVI_ind1
NDVI_basis2_NO3 <- basis_NO3*NDVI_ind2
NDVI_basis3_NO3 <- basis_NO3*NDVI_ind3

colnames(NDVI_basis1_NO3) = paste0("NDVI_basis1_NO3.", colnames(NDVI_basis1_NO3))
colnames(NDVI_basis2_NO3) = paste0("NDVI_basis2_NO3.", colnames(NDVI_basis2_NO3))
colnames(NDVI_basis3_NO3) = paste0("NDVI_basis3_NO3.", colnames(NDVI_basis3_NO3))

model7 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  NDVI_basis1_NO3 +
  NDVI

model8 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  NDVI_basis2_NO3 +
  NDVI

model9 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  NDVI_basis3_NO3 +
  NDVI

model7 <- fitmodel(model7, data) 

coef <- model7$summary.fixed$mean

vcov <- model7$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model7$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

model8 <- fitmodel(model8, data) 

coef <- model8$summary.fixed$mean

vcov <- model8$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model8$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

model9 <- fitmodel(model9, data) 

coef <- model9$summary.fixed$mean

vcov <- model9$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model9$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))


NO3_mat3 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
NDVI_basis1_NH4 <- basis_NH4*NDVI_ind1
NDVI_basis2_NH4 <- basis_NH4*NDVI_ind2
NDVI_basis3_NH4 <- basis_NH4*NDVI_ind3

colnames(NDVI_basis1_NH4) = paste0("NDVI_basis1_NH4.", colnames(NDVI_basis1_NH4))
colnames(NDVI_basis2_NH4) = paste0("NDVI_basis2_NH4.", colnames(NDVI_basis2_NH4))
colnames(NDVI_basis3_NH4) = paste0("NDVI_basis3_NH4.", colnames(NDVI_basis3_NH4))

model10 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  NDVI_basis1_NH4 +
  NDVI

model11 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  NDVI_basis2_NH4 +
  NDVI

model12 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  NDVI_basis3_NH4 +
  NDVI

model10 <- fitmodel(model10, data) 

coef <- model10$summary.fixed$mean

vcov <- model10$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model10$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

model11 <- fitmodel(model11, data) 

coef <- model11$summary.fixed$mean

vcov <- model11$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model11$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

model12 <- fitmodel(model12, data) 

coef <- model12$summary.fixed$mean

vcov <- model12$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model12$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat3 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
NDVI_basis1_BC <- basis_BC*NDVI_ind1
NDVI_basis2_BC <- basis_BC*NDVI_ind2
NDVI_basis3_BC <- basis_BC*NDVI_ind3

colnames(NDVI_basis1_BC) = paste0("NDVI_basis1_BC.", colnames(NDVI_basis1_BC))
colnames(NDVI_basis2_BC) = paste0("NDVI_basis2_BC.", colnames(NDVI_basis2_BC))
colnames(NDVI_basis3_BC) = paste0("NDVI_basis3_BC.", colnames(NDVI_basis3_BC))

model13 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  NDVI_basis1_BC +
  NDVI

model14 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  NDVI_basis2_BC +
  NDVI

model15 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  NDVI_basis3_BC +
  NDVI

model13 <- fitmodel(model13, data) 

coef <- model13$summary.fixed$mean

vcov <- model13$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model13$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

model14 <- fitmodel(model14, data) 

coef <- model14$summary.fixed$mean

vcov <- model14$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model14$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

model15 <- fitmodel(model15, data) 

coef <- model15$summary.fixed$mean

vcov <- model15$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model15$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat3 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
NDVI_basis1_OM <- basis_OM*NDVI_ind1
NDVI_basis2_OM <- basis_OM*NDVI_ind2
NDVI_basis3_OM <- basis_OM*NDVI_ind3

colnames(NDVI_basis1_OM) = paste0("NDVI_basis1_OM.", colnames(NDVI_basis1_OM))
colnames(NDVI_basis2_OM) = paste0("NDVI_basis2_OM.", colnames(NDVI_basis2_OM))
colnames(NDVI_basis3_OM) = paste0("NDVI_basis3_OM.", colnames(NDVI_basis3_OM))

model16 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  NDVI_basis1_OM +
  NDVI

model17 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  NDVI_basis2_OM +
  NDVI

model18 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  NDVI_basis3_OM +
  NDVI

model16 <- fitmodel(model16, data) 

coef <- model16$summary.fixed$mean

vcov <- model16$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model16$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))


model17 <- fitmodel(model17, data) 

coef <- model17$summary.fixed$mean

vcov <- model17$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model17$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))


model18 <- fitmodel(model18, data) 

coef <- model18$summary.fixed$mean

vcov <- model18$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model18$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat3 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$green <- "High NDVI"
PM2.5_mat2$green <- "Middle NDVI"
PM2.5_mat3$green <- "Low NDVI"
SO4_mat1$green <- "High NDVI"
SO4_mat2$green <- "Middle NDVI"
SO4_mat3$green <- "Low NDVI"
NO3_mat1$green <- "High NDVI"
NO3_mat2$green <- "Middle NDVI"
NO3_mat3$green <- "Low NDVI"
NH4_mat1$green <- "High NDVI"
NH4_mat2$green <- "Middle NDVI"
NH4_mat3$green <- "Low NDVI"
BC_mat1$green <- "High NDVI"
BC_mat2$green <- "Middle NDVI"
BC_mat3$green <- "Low NDVI"
OM_mat1$green <- "High NDVI"
OM_mat2$green <- "Middle NDVI"
OM_mat3$green <- "Low NDVI"

PM2.5_mat1$Pollution <- "PM2.5"
PM2.5_mat2$Pollution <- "PM2.5"
PM2.5_mat3$Pollution <- "PM2.5"
SO4_mat1$Pollution <- "Sulfate"
SO4_mat2$Pollution <- "Sulfate"
SO4_mat3$Pollution <- "Sulfate"
NO3_mat1$Pollution <- "Nitrate"
NO3_mat2$Pollution <- "Nitrate"
NO3_mat3$Pollution <- "Nitrate"
NH4_mat1$Pollution <- "Ammonium"
NH4_mat2$Pollution <- "Ammonium"
NH4_mat3$Pollution <- "Ammonium"
BC_mat1$Pollution <- "BC"
BC_mat2$Pollution <- "BC"
BC_mat3$Pollution <- "BC"
OM_mat1$Pollution <- "OM"
OM_mat2$Pollution <- "OM"
OM_mat3$Pollution <- "OM"

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
SO4_mat3$var50 <- "8.5"
NO3_mat3$var50 <- "9"
NH4_mat3$var50 <- "7"
BC_mat3$var50 <- "2"
OM_mat3$var50 <- "10"


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
SO4_mat3$var75 <- "12"
NO3_mat3$var75 <- "16"
NH4_mat3$var75 <- "10.5"
BC_mat3$var75 <- "2.8"
OM_mat3$var75 <- "14"


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
SO4_mat3$var95 <- "17"
NO3_mat3$var95 <- "25"
NH4_mat3$var95 <- "16"
BC_mat3$var95 <- "4.8"
OM_mat3$var95 <- "24"


NDVI_count <- bind_rows(PM2.5_mat1, PM2.5_mat2, PM2.5_mat3, SO4_mat1, SO4_mat2, SO4_mat3, 
                  NO3_mat1, NO3_mat2, NO3_mat3, NH4_mat1, NH4_mat2, NH4_mat3, 
                  BC_mat1, BC_mat2, BC_mat3, OM_mat1, OM_mat2, OM_mat3)

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

#======================== Linear interaction term ==============================
NDVI_ind1 <- data$NDVI - quantile(data$NDVI, p = 0.90) 
NDVI_ind2 <- data$NDVI - quantile(data$NDVI, p = 0.5) 
NDVI_ind3 <- data$NDVI - quantile(data$NDVI, p = 0.10) 

#================== PM2.5 =====================
NDVI_basis1_PM2.5 <- basis_PM2.5*NDVI_ind1
NDVI_basis2_PM2.5 <- basis_PM2.5*NDVI_ind2
NDVI_basis3_PM2.5 <- basis_PM2.5*NDVI_ind3

colnames(NDVI_basis1_PM2.5) = paste0("NDVI_basis1_PM2.5.", colnames(NDVI_basis1_PM2.5))
colnames(NDVI_basis2_PM2.5) = paste0("NDVI_basis2_PM2.5.", colnames(NDVI_basis2_PM2.5))
colnames(NDVI_basis3_PM2.5) = paste0("NDVI_basis3_PM2.5.", colnames(NDVI_basis3_PM2.5))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  NDVI_basis1_PM2.5 +
  NDVI

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  NDVI_basis2_PM2.5 +
  NDVI

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  ns.MT +
  NDVI_basis3_PM2.5 +
  NDVI

model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

model2 <- fitmodel(model2, data) 

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

model3 <- fitmodel(model3, data) 

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum3 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
NDVI_basis1_SO4 <- basis_SO4*NDVI_ind1
NDVI_basis2_SO4 <- basis_SO4*NDVI_ind2
NDVI_basis3_SO4 <- basis_SO4*NDVI_ind3

colnames(NDVI_basis1_SO4) = paste0("NDVI_basis1_SO4.", colnames(NDVI_basis1_SO4))
colnames(NDVI_basis2_SO4) = paste0("NDVI_basis2_SO4.", colnames(NDVI_basis2_SO4))
colnames(NDVI_basis3_SO4) = paste0("NDVI_basis3_SO4.", colnames(NDVI_basis3_SO4))

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  NDVI_basis1_SO4 +
  NDVI

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT +
  NDVI_basis2_SO4 +
  NDVI

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  ns.MT + 
  NDVI_basis3_SO4 +
  NDVI

model4 <- fitmodel(model4, data) 

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model4$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

model5 <- fitmodel(model5, data) 

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model5$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

model6 <- fitmodel(model6, data) 

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model6$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum3 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
NDVI_basis1_NO3 <- basis_NO3*NDVI_ind1
NDVI_basis2_NO3 <- basis_NO3*NDVI_ind2
NDVI_basis3_NO3 <- basis_NO3*NDVI_ind3

colnames(NDVI_basis1_NO3) = paste0("NDVI_basis1_NO3.", colnames(NDVI_basis1_NO3))
colnames(NDVI_basis2_NO3) = paste0("NDVI_basis2_NO3.", colnames(NDVI_basis2_NO3))
colnames(NDVI_basis3_NO3) = paste0("NDVI_basis3_NO3.", colnames(NDVI_basis3_NO3))

model7 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  NDVI_basis1_NO3 +
  NDVI

model8 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  NDVI_basis2_NO3 +
  NDVI

model9 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  ns.MT +
  NDVI_basis3_NO3 +
  NDVI

model7 <- fitmodel(model7, data) 

coef <- model7$summary.fixed$mean

vcov <- model7$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model7$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

model8 <- fitmodel(model8, data) 

coef <- model8$summary.fixed$mean

vcov <- model8$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model8$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

model9 <- fitmodel(model9, data) 

coef <- model9$summary.fixed$mean

vcov <- model9$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model9$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum3 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
NDVI_basis1_NH4 <- basis_NH4*NDVI_ind1
NDVI_basis2_NH4 <- basis_NH4*NDVI_ind2
NDVI_basis3_NH4 <- basis_NH4*NDVI_ind3

colnames(NDVI_basis1_NH4) = paste0("NDVI_basis1_NH4.", colnames(NDVI_basis1_NH4))
colnames(NDVI_basis2_NH4) = paste0("NDVI_basis2_NH4.", colnames(NDVI_basis2_NH4))
colnames(NDVI_basis3_NH4) = paste0("NDVI_basis3_NH4.", colnames(NDVI_basis3_NH4))

model10 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  NDVI_basis1_NH4 +
  NDVI

model11 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  NDVI_basis2_NH4 +
  NDVI

model12 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  ns.MT +
  NDVI_basis3_NH4 +
  NDVI

model10 <- fitmodel(model10, data) 

coef <- model10$summary.fixed$mean

vcov <- model10$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model10$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

model11 <- fitmodel(model11, data) 

coef <- model11$summary.fixed$mean

vcov <- model11$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model11$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

model12 <- fitmodel(model12, data) 

coef <- model12$summary.fixed$mean

vcov <- model12$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model12$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum3 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
NDVI_basis1_BC <- basis_BC*NDVI_ind1
NDVI_basis2_BC <- basis_BC*NDVI_ind2
NDVI_basis3_BC <- basis_BC*NDVI_ind3

colnames(NDVI_basis1_BC) = paste0("NDVI_basis1_BC.", colnames(NDVI_basis1_BC))
colnames(NDVI_basis2_BC) = paste0("NDVI_basis2_BC.", colnames(NDVI_basis2_BC))
colnames(NDVI_basis3_BC) = paste0("NDVI_basis3_BC.", colnames(NDVI_basis3_BC))

model13 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  NDVI_basis1_BC +
  NDVI

model14 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  NDVI_basis2_BC +
  NDVI

model15 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  ns.MT +
  NDVI_basis3_BC +
  NDVI

model13 <- fitmodel(model13, data) 

coef <- model13$summary.fixed$mean

vcov <- model13$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model13$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

model14 <- fitmodel(model14, data) 

coef <- model14$summary.fixed$mean

vcov <- model14$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model14$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

model15 <- fitmodel(model15, data) 

coef <- model15$summary.fixed$mean

vcov <- model15$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model15$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum3 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
NDVI_basis1_OM <- basis_OM*NDVI_ind1
NDVI_basis2_OM <- basis_OM*NDVI_ind2
NDVI_basis3_OM <- basis_OM*NDVI_ind3

colnames(NDVI_basis1_OM) = paste0("NDVI_basis1_OM.", colnames(NDVI_basis1_OM))
colnames(NDVI_basis2_OM) = paste0("NDVI_basis2_OM.", colnames(NDVI_basis2_OM))
colnames(NDVI_basis3_OM) = paste0("NDVI_basis3_OM.", colnames(NDVI_basis3_OM))

model16 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  NDVI_basis1_OM +
  NDVI

model17 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  NDVI_basis2_OM +
  NDVI

model18 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  ns.MT +
  NDVI_basis3_OM +
  NDVI

model16 <- fitmodel(model16, data) 

coef <- model16$summary.fixed$mean

vcov <- model16$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model16$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))


model17 <- fitmodel(model17, data) 

coef <- model17$summary.fixed$mean

vcov <- model17$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model17$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))


model18 <- fitmodel(model18, data) 

coef <- model18$summary.fixed$mean

vcov <- model18$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model18$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum3 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$green <- "High NDVI"
PM2.5_cum2$green <- "Middle NDVI"
PM2.5_cum3$green <- "Low NDVI"
SO4_cum1$green <- "High NDVI"
SO4_cum2$green <- "Middle NDVI"
SO4_cum3$green <- "Low NDVI"
NO3_cum1$green <- "High NDVI"
NO3_cum2$green <- "Middle NDVI"
NO3_cum3$green <- "Low NDVI"
NH4_cum1$green <- "High NDVI"
NH4_cum2$green <- "Middle NDVI"
NH4_cum3$green <- "Low NDVI"
BC_cum1$green <- "High NDVI"
BC_cum2$green <- "Middle NDVI"
BC_cum3$green <- "Low NDVI"
OM_cum1$green <- "High NDVI"
OM_cum2$green <- "Middle NDVI"
OM_cum3$green <- "Low NDVI"

PM2.5_cum1$pollution <- "PM2.5"
PM2.5_cum2$pollution <- "PM2.5"
PM2.5_cum3$pollution <- "PM2.5"
SO4_cum1$pollution <- "Sulfate"
SO4_cum2$pollution <- "Sulfate"
SO4_cum3$pollution <- "Sulfate"
NO3_cum1$pollution <- "Nitrate"
NO3_cum2$pollution <- "Nitrate"
NO3_cum3$pollution <- "Nitrate"
NH4_cum1$pollution <- "Ammonium"
NH4_cum2$pollution <- "Ammonium"
NH4_cum3$pollution <- "Ammonium"
BC_cum1$pollution <- "BC"
BC_cum2$pollution <- "BC"
BC_cum3$pollution <- "BC"
OM_cum1$pollution <- "OM"
OM_cum2$pollution <- "OM"
OM_cum3$pollution <- "OM"


NDVI_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, PM2.5_cum3, SO4_cum1, SO4_cum2, SO4_cum3, 
                  NO3_cum1, NO3_cum2, NO3_cum3, NH4_cum1, NH4_cum2, NH4_cum3, 
                  BC_cum1, BC_cum2, BC_cum3, OM_cum1, OM_cum2, OM_cum3)

#================== Figures export =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM) 

AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM)

data_pollution <- rbind(AH.data, SC.data)

#================== NDVI contour =====================
NDVI_contour <- NDVI_contour %>%
  mutate(green = factor(green, levels = c("Low NDVI", "Middle NDVI", "High NDVI")))

PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- NDVI_contour   %>%
  filter(Pollution == "PM2.5") %>%
  group_by(green) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = green) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.28, 1.44)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ green, scales = "fixed", nrow = 1) +
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

SO4_contour <- NDVI_contour   %>%
  filter(Pollution == "Sulfate") %>%
  group_by(green) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = green) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.28, 1.44)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ green, scales = "fixed", nrow = 1) +
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

NO3_contour <- NDVI_contour   %>%
  filter(Pollution == "Nitrate") %>%
  group_by(green) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = green) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.28, 1.44)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ green, scales = "fixed", nrow = 1) +
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

NH4_contour <- NDVI_contour   %>%
  filter(Pollution == "Ammonium") %>%
  group_by(green) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = green) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.28, 1.44)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ green, scales = "fixed", nrow = 1) +
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

OM_contour <- NDVI_contour   %>%
  filter(Pollution == "OM") %>%
  group_by(green) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = green) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.28, 1.44)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ green, scales = "fixed", nrow = 1) +
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

BC_contour <- NDVI_contour   %>%
  filter(Pollution == "BC") %>%
  group_by(green) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = green) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.28, 1.44)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ green, scales = "fixed", nrow = 1) +
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

pdf("./Figures/Figure S17.pdf", height = 8, width = 14)
All_contour 
dev.off()

#================== NDVI lag =====================
NDVI_contour <- NDVI_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) # 其他值保持不变
  ))

NDVI_contour <- NDVI_contour %>%
  mutate(pollution = factor(Pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

NDVI_contour <- NDVI_contour %>%
  mutate(green = factor(green, levels = c("Low NDVI", "Middle NDVI", "High NDVI")))

col.pal <- c("#a50f15", "#2171b5", "#238b45")

Lag_NDVI <- NDVI_contour %>%
  group_by(green, pollution, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ Pollution, scales = "fixed") +
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

Lag_NDVI <- ggdraw(Lag_NDVI) +
  annotate("text", x = 0.125, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.435, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure S18.pdf", height = 7, width = 14)
Lag_NDVI
dev.off()

#================== NDVI exposure =====================
NDVI_cumcontour <- NDVI_cumcontour %>%
  mutate(Pollution = factor(Pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

NDVI_cumcontour <- NDVI_cumcontour %>%
  mutate(green = factor(green, levels = c("Low NDVI", "Middle NDVI", "High NDVI")))

col.pal <- c("#EE0000FF", "#3B4992FF", "#008B45FF")

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- NDVI_cumcontour %>%
  filter(Pollution == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(green, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.8)) +
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

SO4 <- NDVI_cumcontour %>%
  filter(Pollution == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(green, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.8)) +
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

NO3 <- NDVI_cumcontour %>%
  filter(Pollution == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(green, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.8)) +
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

NH4 <- NDVI_cumcontour %>%
  filter(Pollution == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(green, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.8)) +
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

OM <- NDVI_cumcontour %>%
  filter(Pollution == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(green, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.8)) +
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

BC <- NDVI_cumcontour %>%
  filter(Pollution == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(green, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = green, col = as.factor(green), fill = as.factor(green))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.8)) +
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

exposure_NDVI <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure_NDVI <- ggdraw(exposure_NDVI) +
  annotate("text", x = 0.124, y = 0.975, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 6, hjust = 0.5) +
  annotate("text", x = 0.463, y = 0.975, label = "B.Sulfate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.795, y = 0.975, label = "C.Nitrate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.148, y = 0.475, label = "D.Ammonium (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.449, y = 0.475, label = "E.OM (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.779, y = 0.475, label = "F.BC (μg/m³)", size = 6, hjust = 0.5)

exposure_NDVI<- ggdraw() +
  draw_plot(exposure_NDVI, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.904,      
            y = 0.48,       
            width = 0.04, 
            height = 0.3)


pdf("./Figures/Figure S16.pdf", height = 7, width = 14)
exposure_NDVI
dev.off()

#================== MT =====================

#================== Data import =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  mutate(Pref = Pref + 16)%>%
  dplyr::select(year, month, X, case, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop, regnames, MT) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(year, month, X, case, Pref, CNY, PM2.5, SO4, NO3, NH4, BC, OM, pop, regnames, MT)%>%
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

#======================== DLNM - Lag ==============================
slag = -3
nlag = 12

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

#======================== Linear interaction term ==============================
MT_ind1 <- data$MT - quantile(data$MT, p = 0.90) 
MT_ind2 <- data$MT - quantile(data$MT, p = 0.5) 
MT_ind3 <- data$MT - quantile(data$MT, p = 0.10) 

#================== PM2.5 =====================
MT_basis1_PM2.5 <- basis_PM2.5*MT_ind1
MT_basis2_PM2.5 <- basis_PM2.5*MT_ind2
MT_basis3_PM2.5 <- basis_PM2.5*MT_ind3

colnames(MT_basis1_PM2.5) = paste0("MT_basis1_PM2.5.", colnames(MT_basis1_PM2.5))
colnames(MT_basis2_PM2.5) = paste0("MT_basis2_PM2.5.", colnames(MT_basis2_PM2.5))
colnames(MT_basis3_PM2.5) = paste0("MT_basis3_PM2.5.", colnames(MT_basis3_PM2.5))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  MT_basis1_PM2.5 +
  MT

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  MT_basis2_PM2.5 +
  MT

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  MT_basis3_PM2.5 +
  MT

model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat1 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

model2 <- fitmodel(model2, data) 

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat2 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

model3 <- fitmodel(model3, data) 

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_mat3 <- data.frame(VariableValue = rownames(result.PM2.5$matRRfit),
                         Lag = rep(-3:12, each = nrow(result.PM2.5$matRRfit)), 
                         RR = as.numeric(result.PM2.5$matRRfit), 
                         RR.low = as.numeric(result.PM2.5$matRRlow), 
                         RR.high = as.numeric(result.PM2.5$matRRhigh))

#================== SO4 =====================
MT_basis1_SO4 <- basis_SO4*MT_ind1
MT_basis2_SO4 <- basis_SO4*MT_ind2
MT_basis3_SO4 <- basis_SO4*MT_ind3

colnames(MT_basis1_SO4) = paste0("MT_basis1_SO4.", colnames(MT_basis1_SO4))
colnames(MT_basis2_SO4) = paste0("MT_basis2_SO4.", colnames(MT_basis2_SO4))
colnames(MT_basis3_SO4) = paste0("MT_basis3_SO4.", colnames(MT_basis3_SO4))

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  MT_basis1_SO4 +
  MT

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  MT_basis2_SO4 +
  MT

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  MT_basis3_SO4 +
  MT

model4 <- fitmodel(model4, data) 

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model4$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat1 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

model5 <- fitmodel(model5, data) 

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model5$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat2 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

model6 <- fitmodel(model6, data) 

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model6$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_mat3 <- data.frame(VariableValue = rownames(result.SO4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.SO4$matRRfit)), 
                       RR = as.numeric(result.SO4$matRRfit), 
                       RR.low = as.numeric(result.SO4$matRRlow), 
                       RR.high = as.numeric(result.SO4$matRRhigh))

#================== NO3 =====================
MT_basis1_NO3 <- basis_NO3*MT_ind1
MT_basis2_NO3 <- basis_NO3*MT_ind2
MT_basis3_NO3 <- basis_NO3*MT_ind3

colnames(MT_basis1_NO3) = paste0("MT_basis1_NO3.", colnames(MT_basis1_NO3))
colnames(MT_basis2_NO3) = paste0("MT_basis2_NO3.", colnames(MT_basis2_NO3))
colnames(MT_basis3_NO3) = paste0("MT_basis3_NO3.", colnames(MT_basis3_NO3))

model7 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  MT_basis1_NO3 +
  MT

model8 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  MT_basis2_NO3 +
  MT

model9 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  MT_basis3_NO3 +
  MT

model7 <- fitmodel(model7, data) 

coef <- model7$summary.fixed$mean

vcov <- model7$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model7$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat1 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

model8 <- fitmodel(model8, data) 

coef <- model8$summary.fixed$mean

vcov <- model8$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model8$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat2 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

model9 <- fitmodel(model9, data) 

coef <- model9$summary.fixed$mean

vcov <- model9$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model9$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_mat3 <- data.frame(VariableValue = rownames(result.NO3$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NO3$matRRfit)), 
                       RR = as.numeric(result.NO3$matRRfit), 
                       RR.low = as.numeric(result.NO3$matRRlow), 
                       RR.high = as.numeric(result.NO3$matRRhigh))

#================== NH4 =====================
MT_basis1_NH4 <- basis_NH4*MT_ind1
MT_basis2_NH4 <- basis_NH4*MT_ind2
MT_basis3_NH4 <- basis_NH4*MT_ind3

colnames(MT_basis1_NH4) = paste0("MT_basis1_NH4.", colnames(MT_basis1_NH4))
colnames(MT_basis2_NH4) = paste0("MT_basis2_NH4.", colnames(MT_basis2_NH4))
colnames(MT_basis3_NH4) = paste0("MT_basis3_NH4.", colnames(MT_basis3_NH4))

model10 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  MT_basis1_NH4 +
  MT

model11 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  MT_basis2_NH4 +
  MT

model12 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  MT_basis3_NH4 +
  MT

model10 <- fitmodel(model10, data) 

coef <- model10$summary.fixed$mean

vcov <- model10$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model10$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat1 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

model11 <- fitmodel(model11, data) 

coef <- model11$summary.fixed$mean

vcov <- model11$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model11$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat2 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

model12 <- fitmodel(model12, data) 

coef <- model12$summary.fixed$mean

vcov <- model12$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model12$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_mat3 <- data.frame(VariableValue = rownames(result.NH4$matRRfit),
                       Lag = rep(-3:12, each = nrow(result.NH4$matRRfit)), 
                       RR = as.numeric(result.NH4$matRRfit), 
                       RR.low = as.numeric(result.NH4$matRRlow), 
                       RR.high = as.numeric(result.NH4$matRRhigh))

#================== BC =====================
MT_basis1_BC <- basis_BC*MT_ind1
MT_basis2_BC <- basis_BC*MT_ind2
MT_basis3_BC <- basis_BC*MT_ind3

colnames(MT_basis1_BC) = paste0("MT_basis1_BC.", colnames(MT_basis1_BC))
colnames(MT_basis2_BC) = paste0("MT_basis2_BC.", colnames(MT_basis2_BC))
colnames(MT_basis3_BC) = paste0("MT_basis3_BC.", colnames(MT_basis3_BC))

model13 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  MT_basis1_BC +
  MT

model14 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  MT_basis2_BC +
  MT

model15 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  MT_basis3_BC +
  MT

model13 <- fitmodel(model13, data) 

coef <- model13$summary.fixed$mean

vcov <- model13$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model13$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat1 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

model14 <- fitmodel(model14, data) 

coef <- model14$summary.fixed$mean

vcov <- model14$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model14$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat2 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

model15 <- fitmodel(model15, data) 

coef <- model15$summary.fixed$mean

vcov <- model15$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model15$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_mat3 <- data.frame(VariableValue = rownames(result.BC$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.BC$matRRfit)), 
                      RR = as.numeric(result.BC$matRRfit), 
                      RR.low = as.numeric(result.BC$matRRlow), 
                      RR.high = as.numeric(result.BC$matRRhigh))

#================== OM =====================
MT_basis1_OM <- basis_OM*MT_ind1
MT_basis2_OM <- basis_OM*MT_ind2
MT_basis3_OM <- basis_OM*MT_ind3

colnames(MT_basis1_OM) = paste0("MT_basis1_OM.", colnames(MT_basis1_OM))
colnames(MT_basis2_OM) = paste0("MT_basis2_OM.", colnames(MT_basis2_OM))
colnames(MT_basis3_OM) = paste0("MT_basis3_OM.", colnames(MT_basis3_OM))

model16 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  MT_basis1_OM +
  MT

model17 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  MT_basis2_OM +
  MT

model18 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  MT_basis3_OM +
  MT

model16 <- fitmodel(model16, data) 

coef <- model16$summary.fixed$mean

vcov <- model16$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model16$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))


OM_mat1 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))


model17 <- fitmodel(model17, data) 

coef <- model17$summary.fixed$mean

vcov <- model17$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model17$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat2 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))


model18 <- fitmodel(model18, data) 

coef <- model18$summary.fixed$mean

vcov <- model18$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model18$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_mat3 <- data.frame(VariableValue = rownames(result.OM$matRRfit),
                      Lag = rep(-3:12, each = nrow(result.OM$matRRfit)), 
                      RR = as.numeric(result.OM$matRRfit), 
                      RR.low = as.numeric(result.OM$matRRlow), 
                      RR.high = as.numeric(result.OM$matRRhigh))

#================== Results export =====================
PM2.5_mat1$MT <- "High MT"
PM2.5_mat2$MT <- "Middle MT"
PM2.5_mat3$MT <- "Low MT"
SO4_mat1$MT <- "High MT"
SO4_mat2$MT <- "Middle MT"
SO4_mat3$MT <- "Low MT"
NO3_mat1$MT <- "High MT"
NO3_mat2$MT <- "Middle MT"
NO3_mat3$MT <- "Low MT"
NH4_mat1$MT <- "High MT"
NH4_mat2$MT <- "Middle MT"
NH4_mat3$MT <- "Low MT"
BC_mat1$MT <- "High MT"
BC_mat2$MT <- "Middle MT"
BC_mat3$MT <- "Low MT"
OM_mat1$MT <- "High MT"
OM_mat2$MT <- "Middle MT"
OM_mat3$MT <- "Low MT"

PM2.5_mat1$Pollution <- "PM2.5"
PM2.5_mat2$Pollution <- "PM2.5"
PM2.5_mat3$Pollution <- "PM2.5"
SO4_mat1$Pollution <- "Sulfate"
SO4_mat2$Pollution <- "Sulfate"
SO4_mat3$Pollution <- "Sulfate"
NO3_mat1$Pollution <- "Nitrate"
NO3_mat2$Pollution <- "Nitrate"
NO3_mat3$Pollution <- "Nitrate"
NH4_mat1$Pollution <- "Ammonium"
NH4_mat2$Pollution <- "Ammonium"
NH4_mat3$Pollution <- "Ammonium"
BC_mat1$Pollution <- "BC"
BC_mat2$Pollution <- "BC"
BC_mat3$Pollution <- "BC"
OM_mat1$Pollution <- "OM"
OM_mat2$Pollution <- "OM"
OM_mat3$Pollution <- "OM"

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
SO4_mat3$var50 <- "8.5"
NO3_mat3$var50 <- "9"
NH4_mat3$var50 <- "7"
BC_mat3$var50 <- "2"
OM_mat3$var50 <- "10"


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
SO4_mat3$var75 <- "12"
NO3_mat3$var75 <- "16"
NH4_mat3$var75 <- "10.5"
BC_mat3$var75 <- "2.8"
OM_mat3$var75 <- "14"


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
SO4_mat3$var95 <- "17"
NO3_mat3$var95 <- "25"
NH4_mat3$var95 <- "16"
BC_mat3$var95 <- "4.8"
OM_mat3$var95 <- "24"


MT_count <- bind_rows(PM2.5_mat1, PM2.5_mat2, PM2.5_mat3, SO4_mat1, SO4_mat2, SO4_mat3, 
                      NO3_mat1, NO3_mat2, NO3_mat3, NH4_mat1, NH4_mat2, NH4_mat3, 
                      BC_mat1, BC_mat2, BC_mat3, OM_mat1, OM_mat2, OM_mat3)

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

#======================== Linear interaction term ==============================
MT_ind1 <- data$MT - quantile(data$MT, p = 0.90) 
MT_ind2 <- data$MT - quantile(data$MT, p = 0.5) 
MT_ind3 <- data$MT - quantile(data$MT, p = 0.10) 

#================== PM2.5 =====================
MT_basis1_PM2.5 <- basis_PM2.5*MT_ind1
MT_basis2_PM2.5 <- basis_PM2.5*MT_ind2
MT_basis3_PM2.5 <- basis_PM2.5*MT_ind3

colnames(MT_basis1_PM2.5) = paste0("MT_basis1_PM2.5.", colnames(MT_basis1_PM2.5))
colnames(MT_basis2_PM2.5) = paste0("MT_basis2_PM2.5.", colnames(MT_basis2_PM2.5))
colnames(MT_basis3_PM2.5) = paste0("MT_basis3_PM2.5.", colnames(MT_basis3_PM2.5))

model1 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  MT_basis1_PM2.5 +
  MT

model2 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  MT_basis2_PM2.5 +
  MT

model3 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_PM2.5 +
  MT_basis3_PM2.5 +
  MT

model1 <- fitmodel(model1, data) 

coef <- model1$summary.fixed$mean

vcov <- model1$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model1$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum1 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

model2 <- fitmodel(model2, data) 

coef <- model2$summary.fixed$mean

vcov <- model2$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22]) 

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum2 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

model3 <- fitmodel(model3, data) 

coef <- model3$summary.fixed$mean

vcov <- model3$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_PM2.5", model2$names.fixed[1:22])

result.PM2.5 <- crosspred(basis_PM2.5, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_PM2.5,
                          model.link = "log", bylag = 1, cen = quantile(data$PM2.5, 0.25))

PM2.5_cum3 <- data.frame(VariableValue = rownames(result.PM2.5$cumRRfit),
                         Lag = rep(0:12, each = nrow(result.PM2.5$cumRRfit)), 
                         RR = as.numeric(result.PM2.5$cumRRfit), 
                         RR.low = as.numeric(result.PM2.5$cumRRlow), 
                         RR.high = as.numeric(result.PM2.5$cumRRhigh))

#================== SO4 =====================
MT_basis1_SO4 <- basis_SO4*MT_ind1
MT_basis2_SO4 <- basis_SO4*MT_ind2
MT_basis3_SO4 <- basis_SO4*MT_ind3

colnames(MT_basis1_SO4) = paste0("MT_basis1_SO4.", colnames(MT_basis1_SO4))
colnames(MT_basis2_SO4) = paste0("MT_basis2_SO4.", colnames(MT_basis2_SO4))
colnames(MT_basis3_SO4) = paste0("MT_basis3_SO4.", colnames(MT_basis3_SO4))

model4 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  MT_basis1_SO4 +
  MT

model5 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  MT_basis2_SO4 +
  MT

model6 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_SO4 +
  MT_basis3_SO4 +
  MT

model4 <- fitmodel(model4, data) 

coef <- model4$summary.fixed$mean

vcov <- model4$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model4$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum1 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

model5 <- fitmodel(model5, data) 

coef <- model5$summary.fixed$mean

vcov <- model5$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model5$names.fixed[1:22]) 

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum2 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

model6 <- fitmodel(model6, data) 

coef <- model6$summary.fixed$mean

vcov <- model6$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_SO4", model6$names.fixed[1:22])

result.SO4 <- crosspred(basis_SO4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_SO4,
                        model.link = "log", bylag = 1, cen = quantile(data$SO4, 0.25))

SO4_cum3 <- data.frame(VariableValue = rownames(result.SO4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.SO4$cumRRfit)), 
                       RR = as.numeric(result.SO4$cumRRfit), 
                       RR.low = as.numeric(result.SO4$cumRRlow), 
                       RR.high = as.numeric(result.SO4$cumRRhigh))

#================== NO3 =====================
MT_basis1_NO3 <- basis_NO3*MT_ind1
MT_basis2_NO3 <- basis_NO3*MT_ind2
MT_basis3_NO3 <- basis_NO3*MT_ind3

colnames(MT_basis1_NO3) = paste0("MT_basis1_NO3.", colnames(MT_basis1_NO3))
colnames(MT_basis2_NO3) = paste0("MT_basis2_NO3.", colnames(MT_basis2_NO3))
colnames(MT_basis3_NO3) = paste0("MT_basis3_NO3.", colnames(MT_basis3_NO3))

model7 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  MT_basis1_NO3 +
  MT

model8 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  MT_basis2_NO3 +
  MT

model9 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NO3 +
  MT_basis3_NO3 +
  MT

model7 <- fitmodel(model7, data) 

coef <- model7$summary.fixed$mean

vcov <- model7$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model7$names.fixed[1:22])

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum1 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

model8 <- fitmodel(model8, data) 

coef <- model8$summary.fixed$mean

vcov <- model8$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model8$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum2 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

model9 <- fitmodel(model9, data) 

coef <- model9$summary.fixed$mean

vcov <- model9$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NO3", model9$names.fixed[1:22]) 

result.NO3 <- crosspred(basis_NO3, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NO3,
                        model.link = "log", bylag = 1, cen = quantile(data$NO3, 0.25))

NO3_cum3 <- data.frame(VariableValue = rownames(result.NO3$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NO3$cumRRfit)), 
                       RR = as.numeric(result.NO3$cumRRfit), 
                       RR.low = as.numeric(result.NO3$cumRRlow), 
                       RR.high = as.numeric(result.NO3$cumRRhigh))

#================== NH4 =====================
MT_basis1_NH4 <- basis_NH4*MT_ind1
MT_basis2_NH4 <- basis_NH4*MT_ind2
MT_basis3_NH4 <- basis_NH4*MT_ind3

colnames(MT_basis1_NH4) = paste0("MT_basis1_NH4.", colnames(MT_basis1_NH4))
colnames(MT_basis2_NH4) = paste0("MT_basis2_NH4.", colnames(MT_basis2_NH4))
colnames(MT_basis3_NH4) = paste0("MT_basis3_NH4.", colnames(MT_basis3_NH4))

model10 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  MT_basis1_NH4 +
  MT

model11 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  MT_basis2_NH4 +
  MT

model12 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_NH4 +
  MT_basis3_NH4 +
  MT

model10 <- fitmodel(model10, data) 

coef <- model10$summary.fixed$mean

vcov <- model10$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model10$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum1 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

model11 <- fitmodel(model11, data) 

coef <- model11$summary.fixed$mean

vcov <- model11$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model11$names.fixed[1:22])

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum2 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

model12 <- fitmodel(model12, data) 

coef <- model12$summary.fixed$mean

vcov <- model12$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_NH4", model12$names.fixed[1:22]) 

result.NH4 <- crosspred(basis_NH4, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_NH4,
                        model.link = "log", bylag = 1, cen = quantile(data$NH4, 0.25))

NH4_cum3 <- data.frame(VariableValue = rownames(result.NH4$cumRRfit),
                       Lag = rep(0:12, each = nrow(result.NH4$cumRRfit)), 
                       RR = as.numeric(result.NH4$cumRRfit), 
                       RR.low = as.numeric(result.NH4$cumRRlow), 
                       RR.high = as.numeric(result.NH4$cumRRhigh))

#================== BC =====================
MT_basis1_BC <- basis_BC*MT_ind1
MT_basis2_BC <- basis_BC*MT_ind2
MT_basis3_BC <- basis_BC*MT_ind3

colnames(MT_basis1_BC) = paste0("MT_basis1_BC.", colnames(MT_basis1_BC))
colnames(MT_basis2_BC) = paste0("MT_basis2_BC.", colnames(MT_basis2_BC))
colnames(MT_basis3_BC) = paste0("MT_basis3_BC.", colnames(MT_basis3_BC))

model13 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  MT_basis1_BC +
  MT

model14 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  MT_basis2_BC +
  MT

model15 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_BC +
  MT_basis3_BC +
  MT

model13 <- fitmodel(model13, data) 

coef <- model13$summary.fixed$mean

vcov <- model13$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model13$names.fixed[1:22])

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum1 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

model14 <- fitmodel(model14, data) 

coef <- model14$summary.fixed$mean

vcov <- model14$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model14$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum2 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

model15 <- fitmodel(model15, data) 

coef <- model15$summary.fixed$mean

vcov <- model15$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_BC", model15$names.fixed[1:22]) 

result.BC <- crosspred(basis_BC, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_BC,
                       model.link = "log", bylag = 1, cen = quantile(data$BC, 0.25))

BC_cum3 <- data.frame(VariableValue = rownames(result.BC$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.BC$cumRRfit)), 
                      RR = as.numeric(result.BC$cumRRfit), 
                      RR.low = as.numeric(result.BC$cumRRlow), 
                      RR.high = as.numeric(result.BC$cumRRhigh))

#================== OM =====================
MT_basis1_OM <- basis_OM*MT_ind1
MT_basis2_OM <- basis_OM*MT_ind2
MT_basis3_OM <- basis_OM*MT_ind3

colnames(MT_basis1_OM) = paste0("MT_basis1_OM.", colnames(MT_basis1_OM))
colnames(MT_basis2_OM) = paste0("MT_basis2_OM.", colnames(MT_basis2_OM))
colnames(MT_basis3_OM) = paste0("MT_basis3_OM.", colnames(MT_basis3_OM))

model16 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  MT_basis1_OM +
  MT

model17 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  MT_basis2_OM +
  MT

model18 <- case ~ offset(log(E)) + as.factor(CNY) + 1 + 
  f(Pref, model = "iid", constr = TRUE) + 
  f(year, model = "rw1", constr = TRUE, group = Pref.year, control.group = list(model = "iid")) + 
  f(Pref.cos1, cos_i,model = "iid") + 
  f(Pref.sin1, sin_i, model = "iid") + 
  f(Pref.cos2,cos_i2, model = "iid") + 
  f(Pref.sin2, sin_i2, model = "iid") +
  basis_OM +
  MT_basis3_OM +
  MT

model16 <- fitmodel(model16, data) 

coef <- model16$summary.fixed$mean

vcov <- model16$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model16$names.fixed[1:22])

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum1 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

model17 <- fitmodel(model17, data) 

coef <- model17$summary.fixed$mean

vcov <- model17$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model17$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum2 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

model18 <- fitmodel(model18, data) 

coef <- model18$summary.fixed$mean

vcov <- model18$misc$lincomb.derived.covariance.matrix

indt <- grep("basis_OM", model18$names.fixed[1:22]) 

result.OM <- crosspred(basis_OM, coef = coef[indt], vcov=vcov[indt,indt], cumul=TRUE, at = specific_concentrations_OM,
                       model.link = "log", bylag = 1, cen = quantile(data$OM, 0.25))

OM_cum3 <- data.frame(VariableValue = rownames(result.OM$cumRRfit),
                      Lag = rep(0:12, each = nrow(result.OM$cumRRfit)), 
                      RR = as.numeric(result.OM$cumRRfit), 
                      RR.low = as.numeric(result.OM$cumRRlow), 
                      RR.high = as.numeric(result.OM$cumRRhigh))

#================== Results export =====================
PM2.5_cum1$MT <- "High MT"
PM2.5_cum2$MT <- "Middle MT"
PM2.5_cum3$MT <- "Low MT"
SO4_cum1$MT <- "High MT"
SO4_cum2$MT <- "Middle MT"
SO4_cum3$MT <- "Low MT"
NO3_cum1$MT <- "High MT"
NO3_cum2$MT <- "Middle MT"
NO3_cum3$MT <- "Low MT"
NH4_cum1$MT <- "High MT"
NH4_cum2$MT <- "Middle MT"
NH4_cum3$MT <- "Low MT"
BC_cum1$MT <- "High MT"
BC_cum2$MT <- "Middle MT"
BC_cum3$MT <- "Low MT"
OM_cum1$MT <- "High MT"
OM_cum2$MT <- "Middle MT"
OM_cum3$MT <- "Low MT"

PM2.5_cum1$pollution <- "PM2.5"
PM2.5_cum2$pollution <- "PM2.5"
PM2.5_cum3$pollution <- "PM2.5"
SO4_cum1$pollution <- "Sulfate"
SO4_cum2$pollution <- "Sulfate"
SO4_cum3$pollution <- "Sulfate"
NO3_cum1$pollution <- "Nitrate"
NO3_cum2$pollution <- "Nitrate"
NO3_cum3$pollution <- "Nitrate"
NH4_cum1$pollution <- "Ammonium"
NH4_cum2$pollution <- "Ammonium"
NH4_cum3$pollution <- "Ammonium"
BC_cum1$pollution <- "BC"
BC_cum2$pollution <- "BC"
BC_cum3$pollution <- "BC"
OM_cum1$pollution <- "OM"
OM_cum2$pollution <- "OM"
OM_cum3$pollution <- "OM"

MT_cumcontour <- bind_rows(PM2.5_cum1, PM2.5_cum2, PM2.5_cum3, SO4_cum1, SO4_cum2, SO4_cum3, 
                           NO3_cum1, NO3_cum2, NO3_cum3, NH4_cum1, NH4_cum2, NH4_cum3, 
                           BC_cum1, BC_cum2, BC_cum3, OM_cum1, OM_cum2, OM_cum3)

#================== Figures export =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM) 

AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM)

data_pollution <- rbind(AH.data, SC.data)

#================== NDVI contour =====================
MT_contour <- MT_contour %>%
  mutate(MT = factor(MT, levels = c("Low MT", "Middle MT", "High MT")))

PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5_contour <- MT_contour   %>%
  filter(Pollution == "PM2.5") %>%
  group_by(MT) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = MT) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.55, 1.50)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ MT, scales = "fixed", nrow = 1) +
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

SO4_contour <- MT_contour   %>%
  filter(Pollution == "Sulfate") %>%
  group_by(MT) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = MT) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.55, 1.50)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ MT, scales = "fixed", nrow = 1) +
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

NO3_contour <- MT_contour   %>%
  filter(Pollution == "Nitrate") %>%
  group_by(MT) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = MT) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.55, 1.50)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ MT, scales = "fixed", nrow = 1) +
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

NH4_contour <- MT_contour   %>%
  filter(Pollution == "Ammonium") %>%
  group_by(MT) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = MT) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.55, 1.50)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ MT, scales = "fixed", nrow = 1) +
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

OM_contour <- MT_contour   %>%
  filter(Pollution == "OM") %>%
  group_by(MT) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = MT) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.55, 1.50)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ MT, scales = "fixed", nrow = 1) +
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

BC_contour <- MT_contour   %>%
  filter(Pollution == "BC") %>%
  group_by(MT) %>%
  mutate(VariableValue = as.numeric(VariableValue)) %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = MT) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.55, 1.50)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ MT, scales = "fixed", nrow = 1) +
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

pdf("./Figures/Figure S20.pdf", height = 8, width = 14)
All_contour 
dev.off()

#================== NDVI lag =====================
MT_contour <- MT_contour %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) # 其他值保持不变
  ))

MT_contour <- MT_contour %>%
  mutate(pollution = factor(Pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

MT_contour <- MT_contour %>%
  mutate(MT = factor(MT, levels = c("Low MT", "Middle MT", "High MT")))

col.pal <- c("#238b45", "#2171b5", "#a50f15")

Lag_MT <- MT_contour %>%
  group_by(MT, pollution, VariableValue) %>%
  ggplot(aes(x = Lag, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 0.5, position = position_dodge(width = 0.8))+
  geom_point(size = 1.3, position = position_dodge(width = 0.8))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
  facet_grid(VariableValue ~ Pollution, scales = "fixed") +
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

Lag_MT <- ggdraw(Lag_MT) +
  annotate("text", x = 0.125, y = 0.985, label = expression(paste("A.", PM[2.5])), size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.28, y = 0.985, label = "B.Sulfate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.435, y = 0.985, label = "C.Nitrate", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.589, y = 0.985, label = "D.Ammonium", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.747, y = 0.985, label = "E.OM", size = 4.5, hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.985, label = "F.BC", size = 4.5, hjust = 0.5)

pdf("./Figures/Figure 21.pdf", height = 7, width = 14)
Lag_MT 
dev.off()

#================== MT exposure =====================
MT_cumcontour <- MT_cumcontour %>%
  mutate(Pollution = factor(Pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

MT_cumcontoura <- MT_cumcontour %>%
  mutate(MT = factor(MT, levels = c("Low MT", "Middle MT", "High MT")))

PM2.5_concentration <- data_pollution$PM2.5[data_pollution$PM2.5 <= 95]

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

col.pal <- c("#008B45FF","#3B4992FF", "#EE0000FF")

PM2.5 <- MT_cumcontoura %>%
  filter(Pollution == "PM2.5" & VariableValue <= 95) %>%
  filter(Lag == "3") %>%
  group_by(MT, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.7)) +
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

SO4 <- MT_cumcontoura %>%
  filter(Pollution == "Sulfate" & VariableValue <= 17.5) %>%
  filter(Lag == "3") %>%
  group_by(MT, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.7)) +
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

NO3 <- MT_cumcontoura %>%
  filter(Pollution == "Nitrate" & VariableValue <= 26) %>%
  filter(Lag == "3") %>%
  group_by(MT, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.7)) +
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

NH4 <- MT_cumcontoura %>%
  filter(Pollution == "Ammonium" & VariableValue <= 16.5) %>%
  filter(Lag == "3") %>%
  group_by(MT, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.7)) +
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

OM <- MT_cumcontoura %>%
  filter(Pollution == "OM" & VariableValue <= 25) %>%
  filter(Lag == "3") %>%
  group_by(MT, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.7)) +
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

BC <- MT_cumcontoura %>%
  filter(Pollution == "BC" & VariableValue <= 5) %>%
  filter(Lag == "3") %>%
  group_by(MT, Pollution) %>%
  ggplot(aes(x = VariableValue, y = RR, group = MT, col = as.factor(MT), fill = as.factor(MT))) + 
  geom_ribbon(aes(ymin = RR.low, ymax = RR.high),  alpha = 0.3, color = NA) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_grid(~Pollution, scales = "free") +
  scale_y_continuous(limits = c(0.6, 1.7)) +
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

exposure_MT <- plot_grid(PM2.5, SO4, NO3, NH4, OM, BC, align = "hv", ncol = 3)

exposure_MT <- ggdraw(exposure_MT) +
  annotate("text", x = 0.124, y = 0.975, label = expression(paste("A.", PM[2.5]~ "(μg/m³)")), size = 6, hjust = 0.5) +
  annotate("text", x = 0.463, y = 0.975, label = "B.Sulfate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.795, y = 0.975, label = "C.Nitrate (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.148, y = 0.475, label = "D.Ammonium (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.449, y = 0.475, label = "E.OM (μg/m³)", size = 6, hjust = 0.5) +
  annotate("text", x = 0.779, y = 0.475, label = "F.BC (μg/m³)", size = 6, hjust = 0.5)

exposure_MT <- ggdraw() +
  draw_plot(exposure_MT, 
            x = 0, 
            y = 0, 
            width = 0.92, 
            height = 1) +
  draw_plot(legend, 
            x = 0.915,      
            y = 0.5,       
            width = 0.04, 
            height = 0.3)

pdf("./Figures/Figure S19.pdf", height = 7, width = 14)
final_plot
dev.off()