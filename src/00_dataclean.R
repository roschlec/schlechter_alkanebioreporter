

##  Clean data
library(tidyverse)
library(here)
library(scales)
library(RColorBrewer)
library(ggpubr)

##  Colour palette
palette_invitro <- c("diesel" = "#1F6683", "no_diesel" = "#DD3C51", "sterile" = "#D1C7B5")
palette_d <- c("1" = "#A6BDDB", "6" = "#02818A", "28" = "#014636")
palette_d_plant <- c("0-inoculum" = "#FED976", "0-PFF2" = "#99D8C9", "2-PFF2" = "#41AE76", "7-PFF2" = "#00441B")

##  Labels
trt_invitro <- c("diesel" = "BHB + 1% v/v Diesel",
                 "no_diesel" = "BHB",
                 "sterile" = "BHB (sterile)")

## Create breaks for x-axis (probability plots)
b <- c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999)

##  Functions
#   Probability distribution function
df_ecdf <- function(x){
    temp_ecdf <- ecdf(x)
    sort_prob <- temp_ecdf(sort(x))
    df <- data.frame(ecdf = sort_prob, x = sort(x))
    return(df)
}

##  Load data
#   In vitro OD measurement - Growth in BHB
od <- read.csv(here('data', 'invitro_od.csv')) %>% 
    na.exclude() %>% 
    group_by(replicate, time_d, treatment) %>% 
    summarise(od = mean(od), .groups = "drop") %>% 
    arrange(treatment, time_d)
#   Plotting
od %>% 
    group_by(time_d, treatment) %>% 
    summarise(meanOD = mean(od),
              sdOD = sd(od)) %>% 
    ggplot(aes(x = time_d, y = meanOD, group = treatment))+
    geom_point(data = od, aes(y = od, color = treatment), alpha = 0.8, size = 2)+
    geom_line(linetype = "dashed")+
    geom_errorbar(aes(ymin = meanOD-sdOD, ymax = meanOD+sdOD), width = 0.3)+
    geom_point(aes(color = treatment), fill = "white", pch = 21, size = 2, alpha = 0.9, stroke = 1)+
    scale_color_manual(name = "Treatment", labels = trt_invitro, values = palette_invitro)+
    scale_y_continuous(name = bquote("Optical density ("~OD[600]~") [a.u.]"), limits = c(0,0.8))+
    scale_x_continuous(name = "Time [d]", limits = c(-0.1, 30), breaks = seq(0,28,7))+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_invitro, stroke = 0)))+
    theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(2, "mm"))

invitro_cell <- read.csv(here('data', 'invitro_microscopy.csv'))

label_trt <- invitro_cell %>% 
    group_by(time_d, treatment) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, sep = "-"))

prob_invitro <- invitro_cell %>% 
    group_by(time_d, treatment) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(label_trt$id) %>% 
    bind_rows(.id = "id") %>% 
    separate(id, into = c('time_d', 'trt'), sep = "-", remove = FALSE)
prob_invitro$time_d <- factor(prob_invitro$time_d, levels = c("1", "6", "28"))
prob_invitro$trt <- factor(prob_invitro$trt, levels = c("no_diesel", "diesel"))

prob_invitro %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, color = time_d))+
    facet_grid(cols = vars(trt), labeller = labeller(.cols = trt_invitro))+
    geom_point(size = 2, alpha = 0.75, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d)+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d, stroke = 0)))+
    theme(
        aspect.ratio = 0.75,
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(2, "mm"))

invitro_lb_cell <- read.csv(here('data', 'invitro_microscopy_LB.csv'))

prob_lb <- invitro_lb_cell %>%
    mutate(Label = str_replace(Label, ".czi:.*", "")) %>% 
    separate(Label, into = c("sampling_day", "condition", "rep"), sep = "-") %>% 
    group_by(rep) %>% 
    group_map(~ df_ecdf(.x$RFU)) %>% 
    setNames(paste0("rep", seq(1,7))) %>% 
    bind_rows(.id = "id")

prob_lb %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, group = id))+
    geom_point(size = 2, alpha = 0.5, stroke = 0)+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    theme(
        aspect.ratio = 0.75,
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(2, "mm"))


#   In planta
cfu <- read.csv(here('data', 'inplanta_cfu.csv'))
cfu_null <- read.csv(here('data', 'inplanta_cfu_null.csv')) %>% 
    select(time_d, replicate, nullCFU = logCFU)

cfu %>% 
    group_by(time_d, replicate) %>% 
    summarise(logCFU = mean(logCFU)) %>% 
    #left_join(., cfu_null, by = c("time_d", "replicate")) %>% 
    pivot_longer(-time_d:-replicate, values_to = "logCFU") %>%
    ggplot(aes(x = time_d, y = logCFU, colour = name))+
    geom_point(alpha = 0.5, position = position_jitter(width = 0.02))+
    stat_summary(fun = mean, geom = "line", linetype = "dashed")+
    stat_summary(fun.y = mean, 
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x),
                 aes(shape = name), geom = "pointrange", fill = "white", color = "black", linewidth = 0.8)+
    scale_y_continuous(name = bquote("Bacterial load ["~log[10]~CFU~gFW^-1~"]"), limits = c(4, 8))+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.1, 8), breaks = c(0, 2, 7))+
    scale_shape_manual(values = c(21, 23))+
    scale_colour_manual(values = c("black", "grey"))+
    guides(shape = "none", color = "none")+
    theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(2, "mm"))

plant_cell <- read.csv(here('data', 'inplanta_microscopy.csv'))

label_trt_plant <- plant_cell %>% 
    group_by(time_d, treatment) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, sep = "-"))

prob_plant <- plant_cell %>% 
    group_by(time_d, treatment) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(label_trt_plant$id) %>% 
    bind_rows(.id = "id") %>% 
    separate(id, into = c('time_d', 'trt'), sep = "-", remove = FALSE)
prob_plant$time_d <- factor(prob_plant$time_d, levels = c("0", "2", "7"))
prob_plant$trt <- factor(prob_plant$trt, levels = c("inoculum", "PFF2"))

prob_plant %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, color = id))+
    geom_point(size = 2, alpha = 0.75, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d_plant, labels = c("0 (Inoculum)", "0", "2", "7"))+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d_plant, stroke = 0)))+
    theme(
        aspect.ratio = 0.75,
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(2, "mm"))


plant_cell %>% 
    filter(treatment == "PFF2") %>% 
    mutate(logrfu = log(rfu)) %>% 
    ggplot(aes(x = time_d, y = logrfu, colour = treatment, group=time_d))+
    geom_point(position = position_jitter(width=0.2), alpha = 0.2)+
    geom_boxplot(color = "black", width = 0.4, outlier.alpha = 0, alpha = 0.5)+
    guides(color = "none")+
    theme(
        aspect.ratio = 0.75,
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(2, "mm"))

plant_cell %>% 
    ggplot(aes(x = log(rfu)))+
    facet_grid(cols = vars(time_d), rows = vars(treatment))+
    geom_histogram()+
    theme_bw()
