#### In vitro data ####

##  Libraries
library(tidyverse)
library(here)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(moments)
library(rstatix)
library(ggprism)

source(here('src', '00_dataclean.R'))

##  Growth in BHB   ####
#   OD
od <- 
    read.csv(here('data', 'invitro_od.csv')) %>% 
    na.exclude() %>% 
    group_by(replicate, time_d, treatment) %>% 
    summarise(od = mean(od),
              .groups = "drop") %>% 
    arrange(treatment, time_d)

#  Summary
od %>% 
    group_by(time_d, treatment) %>% 
    summarise(
        mean = mean(od),
        sd = sd(od)) %>% 
    filter(treatment == "diesel") %>% 
    print(n = nrow(.))

#   Plotting
plt_od <-
    od %>% 
    group_by(time_d, treatment) %>% 
    summarise(meanOD = mean(od),
              sdOD = sd(od),
              .groups = "drop") %>% 
    ggplot(aes(x = time_d, y = meanOD, group = treatment))+
    geom_point(data = od, 
               aes(y = od, color = treatment), 
               alpha = 0.8, size = 0.8, stroke = 0)+
    geom_line(linetype = "dashed", lineend = "round", alpha = 0.4)+
    geom_errorbar(aes(ymin = meanOD - sdOD, ymax = meanOD + sdOD), 
                  width = 0.8, linewidth = 0.2)+
    geom_point(aes(color = treatment), 
               fill = "white", pch = 21, size = 0.5, alpha = 0.9, stroke = 1)+
    scale_color_manual(name = "Treatment", labels = trt_invitro, values = palette_invitro)+
    scale_y_continuous(name = bquote("Optical density ("*OD[600]*") [a.u.]"), limits = c(0,0.8))+
    scale_x_continuous(name = "Time [d]", limits = c(-0.1, 30), breaks = seq(0,28,7))+
    guides(color = 
               guide_legend(override.aes = list(alpha = 1, size = 2, fill = palette_invitro, stroke = 0)))+
    theme(
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 4, color = "black"),
        legend.title = element_text(size = 4, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.spacing = unit(0, "pt"),
        legend.key.width = unit(0.5, "lines"),     # Width of the color box
        legend.key.height = unit(0.5, "lines"),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

plt_od
ggsave(here('results', 'figure_1.pdf'), dpi = 600, width = 80, height = 50, units = "mm")
