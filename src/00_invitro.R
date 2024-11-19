#### In vitro data ####

##  Libraries
library(tidyverse)
library(here)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(patchwork)

##  Growth in LB    ####
invitro_lb_cell <- read.csv(here('data', 'invitro_microscopy_LB.csv'))

prob_lb <- invitro_lb_cell %>%
    mutate(Label = str_replace(Label, ".czi:.*", "")) %>% 
    separate(Label, into = c("sampling_day", "condition", "rep"), sep = "-") %>% 
    group_by(rep) %>% 
    group_map(~ df_ecdf(.x$RFU)) %>% 
    setNames(paste0("rep", seq(1,7))) %>% 
    bind_rows(.id = "id")

plot_lb <- prob_lb %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, group = id))+
    geom_point(size = 2, alpha = 0.5, stroke = 0)+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

##  Growth in BHB   ####
#   OD
od <- 
    read.csv(here('data', 'invitro_od.csv')) %>% 
    na.exclude() %>% 
    group_by(replicate, time_d, treatment) %>% 
    summarise(od = mean(od), .groups = "drop") %>% 
    arrange(treatment, time_d)

#   Plotting
plt_od <-
    od %>% 
    group_by(time_d, treatment) %>% 
    summarise(meanOD = mean(od),
              sdOD = sd(od),
              .groups = "drop") %>% 
    ggplot(aes(x = time_d, y = meanOD, group = treatment))+
    geom_point(data = od, aes(y = od, color = treatment), alpha = 0.8, size = 1)+
    geom_line(linetype = "dashed", lineend = "round", alpha = 0.4)+
    geom_errorbar(aes(ymin = meanOD-sdOD, ymax = meanOD+sdOD), width = 0.3)+
    geom_point(aes(color = treatment), fill = "white", pch = 21, size = 2, alpha = 0.9, stroke = 1)+
    scale_color_manual(name = "Treatment", labels = trt_invitro, values = palette_invitro)+
    scale_y_continuous(name = bquote("Optical density ("~OD[600]~") [a.u.]"), limits = c(0,0.8))+
    scale_x_continuous(name = "Time [d]", limits = c(-0.1, 30), breaks = seq(0,28,7))+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_invitro, stroke = 0)))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

#   Microscopy
invitro_cell <- read.csv(here('data', 'invitro_microscopy.csv'))

#   Replicates combined
label_trt <- 
    invitro_cell %>% 
    group_by(time_d, treatment) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, sep = "-"))

#   Probability distribution
prob_invitro <- 
    invitro_cell %>% 
    group_by(time_d, treatment) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(label_trt$id) %>% 
    bind_rows(.id = "id") %>% 
    separate(id, into = c('time_d', 'trt'), sep = "-", remove = FALSE)
prob_invitro$time_d <- factor(prob_invitro$time_d, levels = c("1", "6", "28"))
prob_invitro$trt <- factor(prob_invitro$trt, levels = c("no_diesel", "diesel"))

#   Plotting
plot_microscopy <- 
    prob_invitro %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, color = time_d))+
    facet_grid(cols = vars(trt), labeller = labeller(.cols = trt_invitro))+
    geom_point(size = 1, alpha = 0.75, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d)+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), 
                       labels = scales::percent(b, suffix = ""), 
                       limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d, stroke = 0)))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

#   PER replicates
label_trt_rep <- 
    invitro_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, replicate, sep = "-"))

#   Probability distribution
prob_invitro_rep <- 
    invitro_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(label_trt$id) %>% 
    bind_rows(.id = "id") %>% 
    separate(id, into = c('time_d', 'trt', 'replicate'), sep = "-", remove = FALSE)
prob_invitro_rep$time_d <- factor(prob_invitro_rep$time_d, levels = c("1", "6", "28"))
prob_invitro_rep$trt <- factor(prob_invitro_rep$trt, levels = c("no_diesel", "diesel"))
prob_invitro_rep$replicate <- factor(prob_invitro_rep$replicate)

#   Plotting
plot_microscopy_rep <- 
    prob_invitro_rep %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, color = time_d, group = replicate))+
    facet_grid(cols = vars(trt), labeller = labeller(.cols = trt_invitro))+
    geom_point(size = 1, alpha = 0.5, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d)+
    scale_x_continuous(name = "Cumulative distribution [%]", 
                       breaks = qnorm(b), 
                       labels = scales::percent(b, suffix = ""), 
                       limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d, stroke = 0)))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

#plt_bhb_alkb <-
    invitro_cell %>% 
    group_by(time_d, replicate, treatment) %>% 
    summarise(logrfu = median(log(rfu)), .groups = "drop") %>% 
    ggplot(aes(x = time_d, y = logrfu, group=time_d, color = treatment))+
    geom_point(position = position_jitter(width = 0.1), alpha = 0.8, size = 1)+
    geom_boxplot(color = "black", width = 0.4, linewidth = 0.2, outlier.alpha = 0, alpha = 0.5)+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote(italic("alkB")~"activity ["~log[10]~RFU~"]"))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

## Wrap plots
wrap_plots(plt_od, plot_microscopy, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.5,1))

plot_microscopy_rep

