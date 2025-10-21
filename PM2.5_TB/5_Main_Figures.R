source("1_main model.R")

#================== Main model Figures =====================

SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM) 

AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM)

data_pollution <- rbind(AH.data, SC.data)

#======================== Figure S6 ==============================
PM2.5_concentration <- data_pollution$PM2.5

PM2.5_data <- data.frame(VariableValue = PM2.5_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(41.95169512, 61.53785475, 94.96438359), "Special", "Regular"),
         greenRug = if_else(VariableValue == 29.92511951, "Green", "Other"))

regular_data <- filter(PM2.5_data, Special == "Regular")
special_data <- filter(PM2.5_data, Special == "Special")
green_rug_data <- filter(PM2.5_data, greenRug == "Green")

PM2.5 <- count_mat %>%
  filter(Type == "PM2.5") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
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
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red", size = 1, alpha = 0.8) +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green", size = 1, alpha = 0.8)

SO4_concentration <- data_pollution$SO4

SO4_data <- data.frame(VariableValue = SO4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(8.710487072, 11.92059308, 17.02274153), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047688239, "Green", "Other"))

regular_data <- filter(SO4_data, Special == "Regular")
special_data <- filter(SO4_data, Special == "Special")
green_rug_data <- filter(SO4_data, greenRug == "Green")

SO4 <- count_mat %>%
  filter(Type == "Sulfate") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
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
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red", size = 1, alpha = 0.8) +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green", size = 1, alpha = 0.8)

NO3_concentration <- data_pollution$NO3

NO3_data <- data.frame(VariableValue = NO3_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.238248117, 15.52828443, 25.01219685), "Special", "Regular"),
         greenRug = if_else(VariableValue == 6.047218366, "Green", "Other"))

regular_data <- filter(NO3_data, Special == "Regular")
special_data <- filter(NO3_data, Special == "Special")
green_rug_data <- filter(NO3_data, greenRug == "Green")

NO3 <- count_mat %>%
  filter(Type == "Nitrate") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
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
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red", size = 1, alpha = 0.8) +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green", size = 1, alpha = 0.8)

NH4_concentration <- data_pollution$NH4

NH4_data <- data.frame(VariableValue = NH4_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(6.800077307, 10.34297981, 16.03825827), "Special", "Regular"),
         greenRug = if_else(VariableValue == 4.568814284, "Green", "Other"))

regular_data <- filter(NH4_data, Special == "Regular")
special_data <- filter(NH4_data, Special == "Special")
green_rug_data <- filter(NH4_data, greenRug == "Green")

NH4 <- count_mat %>%
  filter(Type == "Ammonium") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
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
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red", size = 1, alpha = 0.8) +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green", size = 1, alpha = 0.8)

OM_concentration <- data_pollution$OM

OM_data <- data.frame(VariableValue = OM_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(9.855507357, 14.28278867, 24.01885335), "Special", "Regular"),
         greenRug = if_else(VariableValue == 7.449308079, "Green", "Other"))

regular_data <- filter(OM_data, Special == "Regular")
special_data <- filter(OM_data, Special == "Special")
green_rug_data <- filter(OM_data, greenRug == "Green")

OM <- count_mat %>%
  filter(Type == "OM") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
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
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red", size = 1, alpha = 0.8) +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green", size = 1, alpha = 0.8)

BC_concentration <- data_pollution$BC

BC_data <- data.frame(VariableValue = BC_concentration) %>%
  mutate(Special = if_else(VariableValue %in% c(1.925726678, 2.746651586, 4.802384784), "Special", "Regular"),
         greenRug = if_else(VariableValue == 1.44297593, "Green", "Other"))

regular_data <- filter(BC_data, Special == "Regular")
special_data <- filter(BC_data, Special == "Special")
green_rug_data <- filter(BC_data, greenRug == "Green")

BC <- count_mat %>%
  filter(Type == "BC") %>%
  ggplot(aes(x = Lag, y = VariableValue, fill = RR), group = Type) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "#006d2c", high = "#810f7c", limits = c(0.52, 1.33)) +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0), breaks = c(-3, 0, 3, 6, 9, 12)) +
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
  geom_rug(sides = "r", data = special_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "red", size = 1, alpha = 0.8) +
  geom_rug(sides = "r", data = green_rug_data, aes(x = NULL, y = VariableValue, fill = NULL), color = "green", size = 1, alpha = 0.8)

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

pdf("./Figures/Figure S6.pdf", height = 6, width = 12)
count 
dev.off()

#======================== Figure 3 ==============================
count_mat <- count_mat %>%
  filter(VariableValue == var75 | VariableValue == var50 | VariableValue == var95) %>%
  mutate(VariableValue = case_when(
    VariableValue == var75 ~ "75th percentile",
    VariableValue == var50 ~ "50th percentile",
    VariableValue == var95 ~ "95th percentile",
    TRUE ~ as.character(VariableValue) # 其他值保持不变
  ))

count_mat <- count_mat %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#000000", pal_lancet()(5))

Lag = count_mat %>%
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
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 9, 12)) +
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


pdf("./Figures/Figure 3.pdf", height = 8, width = 12)
Lag 
dev.off()

#======================== Figure 4 ==============================
exposure <- exposure %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

# 定义污染物
pollutants <- c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")

combined_data <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(Type == pollutant_name & Lag == "3")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(0, 95, length.out = n_points)
  pollutant_data <- data.frame(
    Type = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

combined_data <- combined_data %>%
  mutate(Type = factor(Type, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#000000", pal_lancet()(5))

combined_plot <- combined_data %>%
  filter(Percentile %in% seq(5, 95, by = 10)) %>%
  group_by(Type) %>%
  ggplot(aes(x = Percentile, y = RR, group = Type, col = as.factor(Type), fill = as.factor(Type))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 2, position = position_dodge(width = 3))+
  geom_point(size = 1.3, position = position_dodge(width = 3))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = col.pal, labels = c(expression(paste(PM[2.5])), "Sulfate", "Nitrate", "Ammonium", "BC", "OM")) +
  scale_fill_manual(values = col.pal, guide = FALSE) +
  scale_y_continuous(limits = c(0.8, 1.6)) +
  scale_x_continuous(
    breaks = seq(5, 95, by = 10),
    labels = function(x) paste0(x)
  ) +
  theme(legend.position = "bottom", 
        legend.justification = c(0, 1),
        plot.title = element_text(hjust = 0, size = 10, face = "bold"))+
  theme_cowplot() +
  labs(
    x = "Percentile", 
    y = "RR",
    color = "",
    fill = "")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 1)) +
  guides(col = guide_legend(nrow = 1, byrow = TRUE)) 


pdf("./Figures/Figure 4.pdf", height = 7, width = 12)
print(combined_plot)
dev.off()

#======================== Figure 5 ==============================

#======================== Sex ==============================
data <- sex_cumcontour

pollutants <- c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")

data_man <- map_dfr(pollutants, function(pollutant_name) {
  
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Male")
  
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  
  
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  
  
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_woman <- map_dfr(pollutants, function(pollutant_name) {
  
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Female")
  
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  
  
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  
  
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_man$group <- "Male"
data_woman$group <- "Female"

data_sex <- bind_rows(data_man, data_woman)

data_sex <- data_sex %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

data_sex <- data_sex %>%
  mutate(group = factor(group, levels = c("Male", "Female")))

col.pal <- c("#EE0000FF", "#3B4992FF")

PM2.5 <- data_sex  %>%
  filter(pollution == "PM2.5") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = Percentile, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 3, position = position_dodge(width = 3))+
  geom_point(size = 1.3, position = position_dodge(width = 3))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  theme(legend.position = c(0.01, 0.9)) +
  scale_x_continuous(
    breaks = seq(5, 95, by = 10),
    labels = function(x) paste0(x)
  ) +
  labs(x = "Percentile", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank())

sex <- ggdraw(PM2.5) +
  annotate("text", x = 0.2, y = 0.96, label = "A. Sex", size = 4, hjust = 0.5) 

#======================== Age ==============================
data <- age_cumcontour

pollutants <- c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")

data_age1 <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "[0, 25)")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_age2 <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "[25, 45)")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_age3 <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "[45, 65)")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_age4 <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "[65, ~)")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_age1$group <- "[0, 25)"
data_age2$group <- "[25, 45)"
data_age3$group <- "[45, 65)"
data_age4$group <- "[65, ~)"

data_age <- bind_rows(data_age1, data_age2, data_age3, data_age4)

data_age <- data_age %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

col.pal <- c("#EE0000FF", "#3B4992FF", "#008B45FF", "#8B008BFF")

PM2.5 <- data_age  %>%
  filter(pollution == "PM2.5") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = Percentile, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 3, position = position_dodge(width = 5))+
  geom_point(size = 1.3, position = position_dodge(width = 5))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  theme(legend.position = c(0.01, 0.9)) +
  scale_x_continuous(
    breaks = seq(5, 95, by = 10),
    labels = function(x) paste0(x)
  ) +
  labs(x = "Percentile", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank())

age <- ggdraw(PM2.5) +
  annotate("text", x = 0.24, y = 0.96, label = "B. Age group", size = 4, hjust = 0.5) 

#======================== IR ==============================
data <- IR_cumcontour

pollutants <- c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")

data_HighIR <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "High IR")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_MiddleIR <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Middle IR")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_LowIR <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Low IR")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_HighIR$group <- "High IR"
data_MiddleIR$group <- "Middle IR"
data_LowIR$group <- "Low IR"

data_IR <- bind_rows(data_HighIR, data_MiddleIR, data_LowIR)

data_IR <- data_IR %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

data_IR <- data_IR %>%
  mutate(group = factor(group, levels = c("Low IR", "Middle IR", "High IR")))

col.pal <- c( "#008B45FF", "#3B4992FF","#EE0000FF")

PM2.5 <- data_IR  %>%
  filter(pollution == "PM2.5") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = Percentile, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 3, position = position_dodge(width = 4))+
  geom_point(size = 1.3, position = position_dodge(width = 4))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  theme(legend.position = c(0.01, 0.9)) +
  scale_x_continuous(
    breaks = seq(5, 95, by = 10),
    labels = function(x) paste0(x)
  ) +
  labs(x = "Percentile", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank())

IR <- ggdraw(PM2.5) +
  annotate("text", x = 0.37, y = 0.96, label = "C. Local transmission intensity", size = 4, hjust = 0.5) 

#======================== NDVI ==============================
data <- NDVI_cumcontour

pollutants <- c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")

data_HighNDVI <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "High NDVI")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_MiddleNDVI <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Middle NDVI")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_LowNDVI <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Low NDVI")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_HighNDVI$group <- "High NDVI"
data_MiddleNDVI$group <- "Middle NDVI"
data_LowNDVI$group <- "Low NDVI"

data_NDVI <- bind_rows(data_HighNDVI, data_MiddleNDVI, data_LowNDVI)

data_NDVI <- data_NDVI %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

data_NDVI <- data_NDVI %>%
  mutate(group = factor(group, levels = c("Low NDVI", "Middle NDVI", "High NDVI")))

col.pal <- c("#EE0000FF", "#3B4992FF", "#008B45FF")

PM2.5 <- data_NDVI  %>%
  filter(pollution == "PM2.5") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = Percentile, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 3, position = position_dodge(width = 4))+
  geom_point(size = 1.3, position = position_dodge(width = 4))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  theme(legend.position = c(0.01, 0.9)) +
  scale_x_continuous(
    breaks = seq(5, 95, by = 10),
    labels = function(x) paste0(x)
  ) +
  labs(x = "Percentile", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank())


NDVI <- ggdraw(PM2.5) +
  annotate("text", x = 0.25, y = 0.96, label = "D. Greenness", size = 4, hjust = 0.5) 

#======================== MT ==============================
data <- MT_cumcontour

pollutants <- c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")

data_HighMT <- map_dfr(pollutants, function(pollutant_name) {
  
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "High MT")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_MiddleMT <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Middle MT")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_LowMT <- map_dfr(pollutants, function(pollutant_name) {
  model_data <- data %>%
    filter(pollution == pollutant_name & Lag == "3" & group == "Low MT")
  if (nrow(model_data) == 0) {
    return(NULL)
  }
  n_points <- nrow(model_data)
  percentiles <- seq(5, 95, length.out = n_points)
  pollutant_data <- data.frame(
    pollution = pollutant_name,
    Percentile = percentiles,
    RR = model_data$RR,
    RR.low = model_data$RR.low,
    RR.high = model_data$RR.high
  )
  return(pollutant_data)
})

data_HighMT$group <- "High MT"
data_MiddleMT$group <- "Middle MT"
data_LowMT$group <- "Low MT"

data_MT <- bind_rows(data_HighMT, data_MiddleMT, data_LowMT)

data_MT <- data_MT %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "OM", "BC")))

data_MT <- data_MT %>%
  mutate(group = factor(group, levels = c("Low MT", "Middle MT", "High MT")))

col.pal <- c("#008B45FF", "#3B4992FF", "#EE0000FF")

PM2.5 <- data_MT  %>%
  filter(pollution == "PM2.5") %>%
  group_by(group, pollution) %>%
  ggplot(aes(x = Percentile, y = RR, group = group, col = as.factor(group), fill = as.factor(group))) + 
  geom_errorbar(aes(ymin=RR.low,ymax=RR.high), width = 3, position = position_dodge(width = 4))+
  geom_point(size = 1.3, position = position_dodge(width = 4))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_cowplot() + 
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  theme(legend.position = c(0.01, 0.9)) +
  scale_x_continuous(
    breaks = seq(5, 95, by = 10),
    labels = function(x) paste0(x)
  ) +
  labs(x = "Percentile", 
       y = "RR",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  theme(strip.text.x = element_blank())

MT <- ggdraw(PM2.5) +
  annotate("text", x = 0.26, y = 0.96, label = "E. Temperature", size = 4, hjust = 0.5) 

Fig_5 <- plot_grid(sex, age, IR, NDVI, MT, align = "hv", ncol = 3)

pdf("./Figures/Figure 5.pdf", height = 7, width = 14)
Fig_5
dev.off()

#================== Map =====================

#================== Shp data import =====================
sheng.shp <- "./Data/Chian.shp"
national_map <- st_read(sheng.shp)

jiuduan_shp <- "./Data/Nine-dash line.shp"
jiuduan_lines <- st_read(jiuduan_shp)

#================== China =====================
national_map <- national_map %>%
  mutate(name = ifelse(pr_name == "四川省", "Sichuan Province",
                       ifelse(pr_name == "安徽省", "Anhui Province", pr_name)))

Chian1 <- ggplot(national_map) +
  geom_sf(aes(fill = name)) +
  geom_sf(data = jiuduan_lines, color = "black") +
  coord_sf(ylim = c(19.5, 55)) + 
  scale_fill_manual(values = c("Anhui Province" = "blue", "Sichuan Province" = "green")) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.45),
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 8)) +
  labs(fill = "")

Chian2 <- Chian1 +
  annotation_north_arrow(location = "tl", 
                         pad_x = unit(0.1, "in"),
                         pad_y = unit(0.1, "in"),
                         style = north_arrow_nautical,
                         which_north = "true") +
  annotation_scale(location = "bl", width_hint = 0.2)

Chian3 <- ggplot(national_map) + 
  geom_sf(data = jiuduan_lines, color = "red") +
  geom_sf() +
  coord_sf(ylim = c(5, 17), xlim = c(108, 120)) +
  theme_void() +
  theme(panel.border = element_rect(fill=NA, color= "grey10", linetype= 2, size= 1),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1.5)

Fig_A <- ggdraw() +
  draw_plot(Chian2) +
  draw_plot(Chian3, x = 0.85, y = 0.2, width = 0.07, height = 0.2)


#================== Anhui =====================
anhui.shp <- "./Data/Anhui.shp"
anhui_map <- st_read(anhui.shp)

anhui_data <- read.csv("./Data/anhui new 13-19.csv")

anhui_data <- anhui_data[c("Pref", "PM2.5", "case", "pop")]

anhui_data$E <- anhui_data$pop*10000

anhui_data$ir <- anhui_data$case / anhui_data$E * 10^5

cityname <- c("合肥市" = "Hefei", "芜湖市" = "Wuhu", "蚌埠市" = "Bengbu",
              "淮南市" = "Huainan", "马鞍山市" = "Maanshan", "淮北市" = "Huaibei",
              "铜陵市" = "Tongling", "安庆市" = "Anqing", "黄山市" = "Huangshan",
              "滁州市" = "Chuzhou", "阜阳市" = "Fuyang", "宿州市" = "Suzhou",
              "六安市" = "Lu'an", "亳州市" = "Bozhou", "池州市" = "Chizhou", "宣城市" = "Xuancheng")

anhui_data_PM2.5 = anhui_data %>%
  group_by(Pref) %>%
  summarise(meanPM2.5 = mean(PM2.5)) 

anhui_data_PM2.5 <- anhui_data_PM2.5 %>%
  mutate(Pref = case_when(
    Pref == "1" ~ "合肥市",
    Pref == "2" ~ "芜湖市",
    Pref == "3" ~ "蚌埠市",
    Pref == "4" ~ "淮南市",
    Pref == "5" ~ "马鞍山市",
    Pref == "6" ~ "淮北市",
    Pref == "7" ~ "铜陵市",
    Pref == "8" ~ "安庆市",
    Pref == "9" ~ "黄山市",
    Pref == "10" ~ "滁州市",
    Pref == "11" ~ "阜阳市",
    Pref == "12" ~ "宿州市",
    Pref == "13" ~ "六安市",
    Pref == "14" ~ "亳州市",
    Pref == "15" ~ "池州市",
    Pref == "16" ~ "宣城市"))  

anhui_data_PM2.5$cityname <- cityname[anhui_data_PM2.5$Pref]

anhui_data_PM2.5 <- rename(anhui_data_PM2.5, ct_name = Pref)

anhui_map_PM2.5 <- left_join(anhui_map, anhui_data_PM2.5, by = "ct_name")

Fig_B <- ggplot(anhui_map_PM2.5) +
  geom_sf() +
  geom_sf(aes(fill = meanPM2.5)) +
  scale_fill_gradient(
    low = "#fff5f0", 
    high = "#67000d", 
    name =  expression(paste(PM[2.5], "(μg/m³)"))) +
  theme_minimal() +
  geom_sf_text(aes(label = cityname), 
               size = 2.5,      
               check_overlap = FALSE,  
               color = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.text.x = element_blank(),      
    axis.text.y = element_blank(),  
    axis.ticks.x = element_blank(),     
    axis.ticks.y = element_blank(),     
    axis.title.x = element_blank(),     
    axis.title.y = element_blank(),     
    legend.position = "bottom")


anhui_data_TB = anhui_data %>%
  group_by(Pref) %>%
  summarise(meanir = mean(ir)) 

anhui_data_TB <- anhui_data_TB %>%
  mutate(Pref = case_when(
    Pref == "1" ~ "合肥市",
    Pref == "2" ~ "芜湖市",
    Pref == "3" ~ "蚌埠市",
    Pref == "4" ~ "淮南市",
    Pref == "5" ~ "马鞍山市",
    Pref == "6" ~ "淮北市",
    Pref == "7" ~ "铜陵市",
    Pref == "8" ~ "安庆市",
    Pref == "9" ~ "黄山市",
    Pref == "10" ~ "滁州市",
    Pref == "11" ~ "阜阳市",
    Pref == "12" ~ "宿州市",
    Pref == "13" ~ "六安市",
    Pref == "14" ~ "亳州市",
    Pref == "15" ~ "池州市",
    Pref == "16" ~ "宣城市")) 

anhui_data_TB$cityname <- cityname[anhui_data_TB$Pref]

anhui_data_TB <- rename(anhui_data_TB, ct_name = Pref)

anhui_map_TB <- left_join(anhui_map, anhui_data_TB, by = "ct_name")

Fig_C <- ggplot(anhui_map_TB) +
  geom_sf() +
  geom_sf(aes(fill = meanir)) +
  scale_fill_gradient(
    low = "#fff5f0", 
    high = "#67000d", 
    name =  "PTB/per 100,000") +
  theme_minimal() +
  geom_sf_text(aes(label = cityname), 
               size = 2.5,      
               check_overlap = FALSE,  
               color = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.text.x = element_blank(),      
    axis.text.y = element_blank(),      
    axis.ticks.x = element_blank(),     
    axis.ticks.y = element_blank(),    
    axis.title.x = element_blank(),     
    axis.title.y = element_blank(),     
    legend.position = "bottom")

#================== Sichuan =====================
sichuan.shp <- "./Data/Sichuan.shp"
sichuan_map <- st_read(sichuan.shp)

sichuan_data <- read.csv("./Data/sichuan new 13-19.csv")

sichuan_data <- sichuan_data[c("Pref", "PM2.5", "case", "pop")]

sichuan_data$E <- sichuan_data$pop*10000

sichuan_data$ir <- sichuan_data$case / sichuan_data$E * 10^5

cityname <- c("成都市" = "Chengdu", "自贡市" = "Zigong", "攀枝花市" = "Panzhihua",
              "泸州市" = "Luzhou", "德阳市" = "Deyang", "绵阳市" = "Mianyang",
              "广元市" = "Guangyuan", "遂宁市" = "Suining", "内江市" = "Neijiang",
              "乐山市" = "Leshan", "南充市" = "Nanchong", "眉山市" = "Meishan",
              "宜宾市" = "Yibin", "广安市" = "guang'an", "达州市" = "Dazhou", "雅安市" = "Yaan",
              "巴中市" = "Bazhong", "资阳市" = "Ziyang", "阿坝藏族羌族自治州" = "Aba", "甘孜藏族自治州" = "Ganzi", "凉山彝族自治州" = "Liangshan")

sichuan_data_PM2.5 = sichuan_data %>%
  group_by(Pref) %>%
  summarise(meanPM2.5 = mean(PM2.5)) 

sichuan_data_PM2.5 <- sichuan_data_PM2.5 %>%
  mutate(Pref = case_when(
    Pref == "1" ~ "成都市",
    Pref == "2" ~ "自贡市",
    Pref == "3" ~ "攀枝花市",
    Pref == "4" ~ "泸州市",
    Pref == "5" ~ "德阳市",
    Pref == "6" ~ "绵阳市",
    Pref == "7" ~ "广元市",
    Pref == "8" ~ "遂宁市",
    Pref == "9" ~ "内江市",
    Pref == "10" ~ "乐山市",
    Pref == "11" ~ "南充市",
    Pref == "12" ~ "眉山市",
    Pref == "13" ~ "宜宾市",
    Pref == "14" ~ "广安市",
    Pref == "15" ~ "达州市",
    Pref == "16" ~ "雅安市",
    Pref == "17" ~ "巴中市",
    Pref == "18" ~ "资阳市",
    Pref == "19" ~ "阿坝藏族羌族自治州",
    Pref == "20" ~ "甘孜藏族自治州",
    Pref == "21" ~ "凉山彝族自治州"))  

sichuan_data_PM2.5$cityname <- cityname[sichuan_data_PM2.5$Pref]

sichuan_data_PM2.5 <- rename(sichuan_data_PM2.5, ct_name = Pref)

sichuan_map_PM2.5 <- left_join(sichuan_map, sichuan_data_PM2.5, by = "ct_name")

sichuan_centroids <- st_centroid(sichuan_map_PM2.5)

sichuan_centroids <- cbind(sichuan_centroids, st_coordinates(sichuan_centroids))

sichuan_centroids_adjusted <- sichuan_centroids %>%
  mutate(
    Y_adjusted = case_when(
      cityname == "Luzhou" ~ Y - 0.3, 
      TRUE ~ Y  
    ),
    X_adjusted = X  
  )

Fig_D <- ggplot(sichuan_map_PM2.5) +
  geom_sf() +
  geom_sf(aes(fill = meanPM2.5)) +
  scale_fill_gradient(
    low = "#fff5f0", 
    high = "#67000d", 
    name = expression(paste(PM[2.5], "(μg/m³)"))) +
  theme_minimal() +
  geom_text(data = sichuan_centroids_adjusted, 
            aes(X_adjusted, Y_adjusted, label = cityname),
            size = 2.5,
            color = "black",
            check_overlap = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom")

sichuan_data_TB = sichuan_data %>%
  group_by(Pref) %>%
  summarise(meanir = mean(ir)) 

sichuan_data_TB <- sichuan_data_TB %>%
  mutate(Pref = case_when(
    Pref == "1" ~ "成都市",
    Pref == "2" ~ "自贡市",
    Pref == "3" ~ "攀枝花市",
    Pref == "4" ~ "泸州市",
    Pref == "5" ~ "德阳市",
    Pref == "6" ~ "绵阳市",
    Pref == "7" ~ "广元市",
    Pref == "8" ~ "遂宁市",
    Pref == "9" ~ "内江市",
    Pref == "10" ~ "乐山市",
    Pref == "11" ~ "南充市",
    Pref == "12" ~ "眉山市",
    Pref == "13" ~ "宜宾市",
    Pref == "14" ~ "广安市",
    Pref == "15" ~ "达州市",
    Pref == "16" ~ "雅安市",
    Pref == "17" ~ "巴中市",
    Pref == "18" ~ "资阳市",
    Pref == "19" ~ "阿坝藏族羌族自治州",
    Pref == "20" ~ "甘孜藏族自治州",
    Pref == "21" ~ "凉山彝族自治州")) 

sichuan_data_TB$cityname <- cityname[sichuan_data_TB$Pref]

sichuan_data_TB <- rename(sichuan_data_TB, ct_name = Pref)

sichuan_map_TB <- left_join(sichuan_map, sichuan_data_TB, by = "ct_name")

sichuan_centroids <- st_centroid(sichuan_map_TB)

sichuan_centroids <- cbind(sichuan_centroids, st_coordinates(sichuan_centroids))

sichuan_centroids_adjusted <- sichuan_centroids %>%
  mutate(
    Y_adjusted = case_when(
      cityname == "Luzhou" ~ Y - 0.3,  
      TRUE ~ Y  
    ),
    X_adjusted = X  
  )

Fig_E <- ggplot(sichuan_map_TB) +
  geom_sf() +
  geom_sf(aes(fill = meanir)) +
  scale_fill_gradient(
    low = "#fff5f0", 
    high = "#67000d", 
    name = "PTB/per 100,000") +
  theme_minimal() +
  geom_text(data = sichuan_centroids_adjusted, 
            aes(X_adjusted, Y_adjusted, label = cityname),
            size = 2.5,
            color = "black",
            check_overlap = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom")

Fig_F1 <- plot_grid(Fig_E, Fig_C, Fig_D, Fig_B,  nrow = 2)

Fig_F2 <- plot_grid(Fig_A, Fig_F1, nrow = 1, labels = c("A", "B"), rel_widths = c(0.5, 0.5))

pdf("./Figures/Figure 1.pdf", height = 8, width = 12)
Fig_F2
dev.off()


#================== Descriptive =====================

#================== Data import =====================
data_an <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM, case, pop, year, month, regnames) %>%
  mutate(province = "Anhui Province")

data_si <- read.csv("./Data/sichuan new 13-19.csv")%>%
  dplyr::select(PM2.5, SO4, NO3, NH4, BC, OM, case, pop, year, month, regnames) %>%
  mutate(province = "Sichuan Province")

data <- bind_rows(data_an, data_si)%>%
  mutate(E = pop*10000)

data$ir <- data$case/data$E*10^5

#================== TB total =====================
data_TB <- data %>%
  group_by(year, month, province) %>%
  summarise(meanir = mean(ir)) 

col.pal <- c("#EE0000FF", "#008B45FF")

data_TB   <- data_TB  %>%
  mutate(province = factor(province, levels = c("Sichuan Province", "Anhui Province")))

Fig_A <- data_TB %>%
  mutate(time = as.Date(paste(year, "-", month,"-1", sep = ""), "%Y-%m-%d")) %>%
  group_by(province) %>%
  ggplot(aes(x = time, y = meanir, group = province, col = as.factor(province))) +
  geom_rect(xmin = ymd("2013-01-01"),xmax = ymd("2014-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2015-01-01"),xmax = ymd("2016-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2017-01-01"),xmax = ymd("2018-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2019-01-01"),xmax = ymd("2019-12-31"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_line(size = 1) +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  theme_cowplot() +
  scale_x_continuous(breaks = seq(as.Date("2013/7/1"), as.Date("2019/7/1"), "years"),
                     labels = seq(2013, 2019, 1)) +
  labs(x = "Time", 
       y = "PTB incidence rate (per 100,000 people)",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = c(0.7, 0.99),    
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

#================== Pollutant total =====================
data_PM2.5 <- data %>%
  group_by(year, month, province) %>%
  summarise(mean = mean(PM2.5))

data_PM2.5$pollution <- "PM2.5"

data_SO4 <- data %>%
  group_by(year, month, province) %>%
  summarise(mean = mean(SO4))

data_SO4$pollution <- "Sulfate"

data_NO3 <- data %>%
  group_by(year, month, province) %>%
  summarise(mean = mean(NO3))

data_NO3$pollution <- "Nitrate"

data_NH4 <- data %>%
  group_by(year, month, province) %>%
  summarise(mean = mean(NH4))

data_NH4$pollution <- "Ammonium"

data_BC <- data %>%
  group_by(year, month, province) %>%
  summarise(mean = mean(BC))

data_BC$pollution <- "BC"

data_OM <- data %>%
  group_by(year, month, province) %>%
  summarise(mean = mean(OM))

data_OM$pollution <- "OM"

data_pollution <- bind_rows(data_PM2.5, data_SO4, data_NO3, data_NH4, data_BC, data_OM)

data_pollution <- data_pollution %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "BC", "OM")))

col.pal <- c("#000000", pal_lancet()(5))

data_pollution <- data_pollution  %>%
  mutate(province = factor(province, levels = c("Sichuan Province", "Anhui Province")))

Fig_B <- data_pollution %>%
  group_by(pollution) %>%
  mutate(time = as.Date(paste(year, "-", month,"-1", sep = ""), "%Y-%m-%d")) %>%
  group_by(province) %>%
  ggplot(aes(x = time, y = mean, group = pollution, col = as.factor(pollution))) +
  geom_rect(xmin = ymd("2013-01-01"),xmax = ymd("2014-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2015-01-01"),xmax = ymd("2016-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2017-01-01"),xmax = ymd("2018-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2019-01-01"),xmax = ymd("2019-12-31"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_line(size = 0.8) +
  scale_color_manual(values = col.pal, labels = c(expression(paste(PM[2.5])), "Sulfate", "Nitrate", "Ammonium", "BC", "OM")) +
  scale_fill_manual(values = col.pal) +
  theme_cowplot() +
  scale_y_sqrt(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160)) +
  facet_wrap(~province, scales = "free", nrow = 2) +
  scale_x_continuous(breaks = seq(as.Date("2013/7/1"), as.Date("2019/7/1"), "years"),
                     labels = seq(2013, 2019, 1)) +
  labs(x = "Time", 
       y = "Average concentration (μg/m³)",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0 , 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  theme(strip.text.x = element_text(hjust = 0, vjust = 1)) +
  guides(col = guide_legend(nrow = 1, byrow = TRUE, 
                            override.aes = list(size = 4)))  

Fig_C <- plot_grid(Fig_A, Fig_B, nrow = 1, labels = c("A", "B"))

pdf("./Figures/Fig2.pdf", height = 6, width = 12)
Fig_C 
dev.off()

#================== TB City=====================
data_TB <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(meanir = mean(ir)) 

Fig_anhui <- data_TB %>%
  filter(province == "Anhui Province") %>%
  mutate(time = as.Date(paste(year, "-", month,"-1", sep = ""), "%Y-%m-%d"),
         regnames = str_to_title(regnames)) %>%
  group_by(regnames) %>%
  ggplot(aes(x = time, y = meanir, group = regnames, col = as.factor(regnames))) +
  geom_rect(xmin = ymd("2013-01-01"),xmax = ymd("2014-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2015-01-01"),xmax = ymd("2016-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2017-01-01"),xmax = ymd("2018-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2019-01-01"),xmax = ymd("2019-12-31"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_line(size = 1) +
  
  theme_cowplot() +
  facet_wrap(~regnames) +
  scale_x_continuous(breaks = seq(as.Date("2013/7/1"), as.Date("2019/7/1"), "2 years"),
                     labels = seq(2013, 2019, 2)) +
  labs(x = "Time", 
       y = "PTB incidence rate (per 100,000 people)",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

pdf("./Figures/FigS1.pdf", height = 6, width = 8)
Fig_anhui
dev.off()

Fig_sichuan <- data_TB %>%
  filter(province == "Sichuan Province") %>%
  mutate(time = as.Date(paste(year, "-", month,"-1", sep = ""), "%Y-%m-%d"),
         regnames = str_to_title(regnames)) %>%
  group_by(regnames) %>%
  ggplot(aes(x = time, y = meanir, group = regnames, col = as.factor(regnames))) +
  geom_rect(xmin = ymd("2013-01-01"),xmax = ymd("2014-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2015-01-01"),xmax = ymd("2016-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2017-01-01"),xmax = ymd("2018-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2019-01-01"),xmax = ymd("2019-12-31"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_line(size = 1) +
  scale_y_sqrt(breaks = c(0, 10, 20, 30)) +
  theme_cowplot() +
  facet_wrap(~regnames, scales = "fixed", nrow = 7) +
  scale_x_continuous(breaks = seq(as.Date("2013/7/1"), as.Date("2019/7/1"), "2 years"),
                     labels = seq(2013, 2019, 2)) +
  labs(x = "Time", 
       y = "PTB incidence rate (per 100,000 people)",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "none",   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

pdf("./Figures/FigS2.pdf", height = 8, width = 8)
Fig_sichuan
dev.off()


#================== Pollutant City=====================
data_PM2.5 <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(mean = mean(PM2.5))

data_PM2.5$pollution <- "PM2.5"

data_SO4 <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(mean = mean(SO4))

data_SO4$pollution <- "Sulfate"

data_NO3 <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(mean = mean(NO3))

data_NO3$pollution <- "Nitrate"

data_NH4 <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(mean = mean(NH4))

data_NH4$pollution <- "Ammonium"

data_BC <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(mean = mean(BC))

data_BC$pollution <- "BC"

data_OM <- data %>%
  group_by(year, month, province, regnames) %>%
  summarise(mean = mean(OM))

data_OM$pollution <- "OM"

data_pollution <- bind_rows(data_PM2.5, data_SO4, data_NO3, data_NH4, data_BC, data_OM)

data_pollution <- data_pollution %>%
  mutate(pollution = factor(pollution, levels = c("PM2.5", "Sulfate", "Nitrate", "Ammonium", "BC", "OM")))

col.pal <- c("#000000", pal_lancet()(5))

Fig_anhui <- data_pollution %>%
  filter(province == "Anhui Province") %>%
  group_by(pollution) %>%
  mutate(time = as.Date(paste(year, "-", month,"-1", sep = ""), "%Y-%m-%d"),
         regnames = str_to_title(regnames)) %>%
  ggplot(aes(x = time, y = mean, group = pollution, col = as.factor(pollution))) +
  geom_rect(xmin = ymd("2013-01-01"),xmax = ymd("2014-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2015-01-01"),xmax = ymd("2016-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2017-01-01"),xmax = ymd("2018-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2019-01-01"),xmax = ymd("2019-12-31"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_line(size = 0.5) +
  scale_color_manual(values = col.pal, labels = c(expression(paste(PM[2.5])), "Sulfate", "Nitrate", "Ammonium", "BC", "OM")) +
  scale_fill_manual(values = col.pal) +
  theme_cowplot() +
  facet_wrap(~regnames) +
  scale_y_sqrt(breaks = c(0, 40, 80, 120, 160)) +
  scale_x_continuous(breaks = seq(as.Date("2013/7/1"), as.Date("2019/7/1"), "2 years"),
                     labels = seq(2013, 2019, 2)) +
  labs(x = "Time", 
       y = "Average concentration (μg/m³)",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0 , 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  guides(col = guide_legend(nrow = 1, byrow = TRUE))

pdf("./Figures/FigS3.pdf", height = 8, width = 8)
Fig_anhui
dev.off()

Fig_sichuan <- data_pollution %>%
  filter(province == "Sichuan Province") %>%
  group_by(pollution) %>%
  mutate(time = as.Date(paste(year, "-", month,"-1", sep = ""), "%Y-%m-%d"),
         regnames = str_to_title(regnames)) %>%
  ggplot(aes(x = time, y = mean, group = pollution, col = as.factor(pollution))) +
  geom_rect(xmin = ymd("2013-01-01"),xmax = ymd("2014-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2015-01-01"),xmax = ymd("2016-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2017-01-01"),xmax = ymd("2018-01-01"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_rect(xmin = ymd("2019-01-01"),xmax = ymd("2019-12-31"),
            ymin = -Inf, ymax = Inf, fill = "gray90", col = "gray90") +
  geom_line(size = 0.5) +
  scale_color_manual(values = col.pal, labels = c(expression(paste(PM[2.5])), "Sulfate", "Nitrate", "Ammonium", "BC", "OM")) +
  scale_fill_manual(values = col.pal) +
  theme_cowplot() +
  facet_wrap(~regnames, nrow = 7) +
  scale_y_sqrt(breaks = c(0, 40, 80, 120, 160)) +
  scale_x_continuous(breaks = seq(as.Date("2013/7/1"), as.Date("2019/7/1"), "2 years"),
                     labels = seq(2013, 2019, 2)) +
  labs(x = "Time", 
       y = "Average concentration (μg/m³)",
       col = "",
       fill = "") +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = c(0 , 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))+
  guides(col = guide_legend(nrow = 1, byrow = TRUE))

pdf("./Figures/FigS4.pdf", height = 10, width = 8)
Fig_sichuan
dev.off()


#================== Spearman correlation =====================

#================== Data import =====================
SC.data <- read.csv("./Data/sichuan new 13-19.csv") %>%
  mutate(Pref = Pref + 16)%>%
  dplyr::select(Pref, PM2.5, SO4, NO3, NH4, BC, OM, MT, RH, NDVI) %>%
  mutate(prov = "SC")
AH.data <- read.csv("./Data/anhui new 13-19.csv") %>%
  dplyr::select(Pref, PM2.5, SO4, NO3, NH4, BC, OM, MT, RH, NDVI)%>%
  mutate(prov = "AH")

data <- rbind(AH.data, SC.data) 

corresults <- data %>%
  group_by(Pref) %>%
  summarize(
    corAA = cor(PM2.5, PM2.5, use = "complete.obs"),
    corAB = cor(PM2.5, SO4, use = "complete.obs"),
    corAC = cor(PM2.5, NO3, use = "complete.obs"),
    corAD = cor(PM2.5, NH4, use = "complete.obs"),
    corAE = cor(PM2.5, BC, use = "complete.obs"),
    corAF = cor(PM2.5, OM, use = "complete.obs"),
    corAG = cor(PM2.5, MT, use = "complete.obs"),
    corAH = cor(PM2.5, RH, use = "complete.obs"),
    corAI = cor(PM2.5, NDVI, use = "complete.obs"),
    corBB = cor(SO4, SO4, use = "complete.obs"),
    corBC = cor(SO4, NO3, use = "complete.obs"),
    corBD = cor(SO4, NH4, use = "complete.obs"),
    corBE = cor(SO4, BC, use = "complete.obs"),
    corBF = cor(SO4, OM, use = "complete.obs"),
    corBG = cor(SO4, MT, use = "complete.obs"),
    corBH = cor(SO4, RH, use = "complete.obs"),
    corBI = cor(SO4, NDVI, use = "complete.obs"),
    corCC = cor(NO3, NO3, use = "complete.obs"),
    corCD = cor(NO3, NH4, use = "complete.obs"),
    corCE = cor(NO3, BC, use = "complete.obs"),
    corCF = cor(NO3, OM, use = "complete.obs"),
    corCG = cor(NO3, MT, use = "complete.obs"),
    corCH = cor(NO3, RH, use = "complete.obs"),
    corCI = cor(NO3, NDVI, use = "complete.obs"),
    corDD = cor(NH4, NH4, use = "complete.obs"),
    corDE = cor(NH4, BC, use = "complete.obs"),
    corDF = cor(NH4, OM, use = "complete.obs"),
    corDG = cor(NH4, MT, use = "complete.obs"),
    corDH = cor(NH4, RH, use = "complete.obs"),
    corDI = cor(NH4, NDVI, use = "complete.obs"),
    corEE = cor(BC, BC, use = "complete.obs"),
    corEF = cor(BC, OM, use = "complete.obs"),
    corEG = cor(BC, MT, use = "complete.obs"),
    corEH = cor(BC, RH, use = "complete.obs"),
    corEI = cor(BC, NDVI, use = "complete.obs"),
    corFF = cor(OM, OM, use = "complete.obs"),
    corFG = cor(OM, MT, use = "complete.obs"),
    corFH = cor(OM, RH, use = "complete.obs"),
    corFI = cor(OM, NDVI, use = "complete.obs"),
    corGG = cor(MT, MT, use = "complete.obs"),
    corGH = cor(MT, RH, use = "complete.obs"),
    corGI = cor(MT, NDVI, use = "complete.obs"),
    corHH = cor(RH, RH, use = "complete.obs"),
    corHI = cor(RH, NDVI, use = "complete.obs"),
    corII = cor(NDVI, NDVI, use = "complete.obs"))

summary_stats <- corresults %>%
  summarize(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),  
    sd = ~ sd(.x, na.rm = TRUE),      
    median = ~ median(.x, na.rm = TRUE), 
    n = ~ sum(!is.na(.x))             
  ))) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(mean|sd|median|n)"
  ) %>%
  mutate(
    lower_ci = mean - 1.96 * (sd / sqrt(n)), 
    upper_ci = mean + 1.96 * (sd / sqrt(n)) 
  )

median_and_variable <- summary_stats %>%
  dplyr::select(variable, median) %>%
  slice(-1)  

median_and_variable <- median_and_variable %>%
  mutate(
    var1 = substr(variable, 4, 4),  
    var2 = substr(variable, 5, 5)  
  )

var_names <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
cor_matrix <- matrix(NA, nrow = 9, ncol = 9, dimnames = list(var_names, var_names))

for (i in 1:nrow(median_and_variable)) {
  row_index <- which(var_names == median_and_variable$var1[i])
  col_index <- which(var_names == median_and_variable$var2[i])
  cor_matrix[row_index, col_index] <- median_and_variable$median[i]
}
name_mapping <- c(
  "A" = "PM2.5",
  "B" = "Sulfate",
  "C" = "Nitrate",
  "D" = "Ammonium",
  "E" = "BC",
  "F" = "OM",
  "G" = "MT",
  "H" = "RH",
  "I" = "NDVI"
)

rownames(cor_matrix) <- name_mapping[rownames(cor_matrix)]
colnames(cor_matrix) <- name_mapping[colnames(cor_matrix)]

A = ggcorrplot(cor_matrix, 
               type = "full", 
               lab = TRUE,
               colors = c("darkred", "white", "steelblue"))   


ggsave("./Figures/FigS5.pdf", plot = A, width = 6, height = 6, dpi = 600)
