#### In planta data - REP 1 ####

##  Libraries
library(tidyverse)
library(here)
library(scales)
library(RColorBrewer)
library(ggpubr)

#   In planta
#   CFU data
cfu <- read.csv(here('data', 'inplanta_cfu.csv'))

plt_cfu <- 
    cfu %>% 
    group_by(time_d, replicate) %>% 
    summarise(logCFU = mean(logCFU), .groups = "drop") %>% 
    pivot_longer(-time_d:-replicate, values_to = "logCFU") %>%
    ggplot(aes(x = time_d, y = logCFU, colour = name))+
    geom_point(alpha = 0.5, position = position_jitter(width = 0.02))+
    stat_summary(fun = mean, geom = "line", linetype = "dashed")+
    stat_summary(fun.y = mean, 
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x),
                 aes(shape = name), geom = "pointrange", fill = "white", color = "black", linewidth = 0.8)+
    scale_y_continuous(name = bquote("Bacterial load ["~log[10]~CFU~gFW^-1~"]"), limits = c(4, 8))+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_shape_manual(values = c(21, 23))+
    scale_colour_manual(values = c("black", "grey"))+
    guides(shape = "none", color = "none")+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

#   Microscopy
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

plt_prob <- 
    prob_plant %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, color = id))+
    geom_point(size = 1, alpha = 0.75, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d_plant, labels = c("0 (Inoculum)", "0", "2", "7"))+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d_plant, stroke = 0)))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

plt_alkb <-
    plant_cell %>% 
    filter(treatment == "PFF2") %>% 
    group_by(time_d, replicate) %>% 
    summarise(logrfu = median(log(rfu)), .groups = "drop") %>% 
    ggplot(aes(x = time_d, y = logrfu, group=time_d))+
    geom_point(position = position_jitter(width = 0.1), alpha = 0.8, size = 1)+
    geom_boxplot(color = "black", width = 0.4, linewidth = 0.2, outlier.alpha = 0, alpha = 0.5)+
    scale_x_continuous(name = "Time [dpi]", 
                       limits = c(-0.2, 8), 
                       breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote(italic("alkB")~"activity ["~log[10]~RFU~"]"),
                       limits = c(-2, 1.5))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))


#   Prob by rep
label_trt_plant2 <- plant_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, replicate, sep = "-"))

prob_plant_replicate <- plant_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(label_trt_plant2$id) %>% 
    bind_rows(.id = "id") %>% 
    separate(id, into = c('time_d', 'trt', 'rep'), sep = "-", remove = FALSE)
prob_plant_replicate$time_d <- factor(prob_plant_replicate$time_d, levels = c("0", "2", "7"))
prob_plant_replicate$trt <- factor(prob_plant_replicate$trt, levels = c("inoculum", "PFF2"))

plt_prob_rep <- prob_plant_replicate %>% 
    filter(ecdf < 1 & trt == "PFF2") %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, fill = time_d))+
    facet_grid(cols = vars(time_d))+
    geom_point(size = 1, alpha = 0.8, stroke = 0.1, pch = 21)+
    scale_fill_manual(name = "Time [dpi]", 
                      values = c("#FED976", "#41AE76", "#00441B"), 
                      labels = c("0", "2", "7"))+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(fill = "none")+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

plt_alkb_rep <- prob_plant_replicate %>% 
    filter(ecdf > 0.9 & trt == "PFF2") %>% 
    group_by(time_d, rep) %>%
    summarise(logrfu = median(log(x)), .groups = "drop") %>% 
    mutate(time_d = case_when(
        time_d == "0" ~ 0,
        time_d == "2" ~ 2,
        time_d == "7" ~ 7)) %>% 
    ggplot(aes(x = time_d, y = logrfu, fill = as.factor(time_d)))+
    geom_point(position = position_jitter(width = 0.1), alpha = 0.8, size = 1)+
    geom_boxplot(color = "black", width = 0.4, linewidth = 0.2, outlier.alpha = 0, alpha = 0.5)+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote(italic("alkB")~"activity ["~log[10]~RFU~"]"))+
    scale_fill_manual(name = "Time [dpi]", 
                      values = c("#FED976", "#41AE76", "#00441B"), 
                      labels = c("0", "2", "7"))+
    guides(fill = "none")+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))


## Wrap plots
wrap_plots(plt_cfu, plt_prob + plt_alkb, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.5,1))

plots <- plt_prob_rep + plt_alkb_rep + plot_layout(widths = c(3,1))

wrap_plots(plt_cfu, plots, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.5,1))
