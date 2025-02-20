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

##  Growth in LB    ####
invitro_lb_cell <- 
    read.csv(here('data', 'invitro_microscopy_LB.csv')) %>% 
    filter(RFU < 3)

prob_lb <- 
    invitro_lb_cell %>%
    mutate(Label = str_replace(Label, ".czi:.*", "")) %>% 
    separate(Label, into = c("sampling_day", "condition", "rep"), sep = "-") %>% 
    group_by(rep) %>% 
    group_map(~ df_ecdf(.x$RFU)) %>% 
    setNames(paste0("rep", seq(1,7))) %>% 
    bind_rows(.id = "id") %>% 
    mutate(medium = "lb")

#   Summarise
summary_lb <-
    invitro_lb_cell %>% 
    summarise(
        n = n(),
        mean_rfu = mean(RFU),
        sd_rfu = sd(RFU),
        iqr_rfu = IQR(RFU),
        q3 = quantile(RFU, probs = 0.75),
        I = q3 + 1.5*iqr_rfu,
        mean_sd = mean_rfu + 2*sd_rfu,
        skewness = skewness(RFU),
        kurtosis = kurtosis(RFU))
summary_lb

threshold <- 
    summary_lb %>% 
    pull(I)

#   Plot
plot_lb <- 
    prob_lb %>% 
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

#   Summary
#   Median and Interquantile range (IQR)
invitro_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    summarise(
        n = n(),
        median = median(rfu),
        q1 = quantile(rfu, probs = 0.25),
        q3 = quantile(rfu, probs = 0.75),
        IQR = paste0(sprintf("%.1f", q1), "--", sprintf("%.1f", q3)))

invitro_cell %>% 
    group_by(treatment, time_d, replicate) %>% 
    mutate(rel_fl = rfu/threshold) %>% 
    summarise(
        n = n(),
        median = median(rfu),
        IQR = IQR(rfu),
        # Proportion of highly fluorescent cells
        prop_high = mean(rfu > threshold) * 100,
        rel_fl = mean(rel_fl),
        # Skewness and Kurtosis
        skewness = skewness(rfu),
        kurtosis = kurtosis(rfu),
        # Coefficient of Variation
        cv = (sd(rfu) / mean(rfu)) * 100)

#   Stats
invitro_cell_sum <- 
    invitro_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    summarise(m_rfu = median(rfu),
              .groups = "drop")

#   Kruskal-Walis on treatments by day  
invitro_cell_sum %>% 
    group_by(time_d) %>% 
    kruskal_test(m_rfu ~ treatment)

invitro_cell_sum %>% 
    group_by(time_d) %>% 
    dunn_test(m_rfu ~ treatment, p.adjust.method = "BH")

#   Kruskal-Walis on day by treatment
invitro_cell_sum %>% 
    group_by(treatment) %>% 
    kruskal_test(m_rfu ~ time_d)

invitro_cell_sum %>% 
    group_by(treatment) %>% 
    dunn_test(m_rfu ~ time_d, p.adjust.method = "BH")

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
plot_microscopy_lb_bhb <-
    prob_invitro %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x))+
    facet_grid(cols = vars(trt), labeller = labeller(.cols = trt_invitro))+
    geom_point(data = prob_lb %>% filter(ecdf < 1), 
               aes(shape = medium),
               size = 1.5, alpha = 0.3, stroke = 0.5)+
    geom_point(aes(color = time_d), size = 1.5, alpha = 0.75, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d)+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), 
                       labels = scales::percent(b, suffix = ""), 
                       limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    scale_shape_manual(name = "Control", values = 18, label = "LB")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d, stroke = 0)),
           shape = guide_legend(override.aes = list(alpha = 1, size = 4, stroke = 0)))+
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

p_values_summary <-
    invitro_cell_sum %>% 
    mutate(rel_rfu = m_rfu/threshold) %>% 
    group_by(time_d) %>% 
    dunn_test(rel_rfu ~ treatment, p.adjust.method = "BH") %>% 
    add_xy_position(x = "time_d", dodge = 3) %>% 
    mutate(x = time_d,
           xmin = x - 0.8,
           xmax = x + 0.8)

plot_microscopy_summary <-
    invitro_cell_sum %>% 
    ggplot(aes(x = time_d, y = m_rfu/threshold))+
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.2)+
    geom_line(aes(group = interaction(time_d, treatment)),
              position = position_dodge(width = 3),
              linewidth = 0.5)+
    geom_point(aes(fill = treatment, group = treatment), 
               position = position_dodge(width = 3),
               pch = 21, stroke = 0.5)+
    scale_y_continuous(name = bquote("Relative"~italic("alkB")~"activity"~bgroup("(",over(RFU[Treatment],RFU[LB]),")")),
                       limits = c(0, 17))+
    scale_x_continuous(name = "Time [d]", limits = c(-0.1, 30), breaks = c(1, 6, 28))+
    scale_fill_manual(name = "Treatment", labels = trt_invitro[1:2], values = palette_invitro[1:2])+
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_invitro[1:2], stroke = 0)))+
    theme_bw()+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = NA),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25)) +
    add_pvalue(p_values_summary, 
               xmin = "xmin", xmax = "xmax",
               label.size = 5,
               tip.length = 0.01,
               bracket.size = 0.3)

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
    setNames(label_trt_rep$id) %>% 
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
    geom_point(size = 2, alpha = 0.5, stroke = 0)+
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

## Wrap plots
wrap_plots(plt_od, plot_microscopy_lb_bhb, plot_microscopy_summary, ncol = 1)+
    plot_annotation(tag_levels = "A") &
    theme(legend.justification = "left",
          axis.title.y = element_text(hjust = 0.5, vjust = 0.2))
ggsave(here('results', 'figure1.png'), dpi = 300, width = 6, height = 8)

plot_lb
ggsave(here('results', 'figureS1.png'), dpi = 300, width = 4, height = 4)

plot_microscopy_rep
ggsave(here('results', 'figureS2.png'), dpi = 300, width = 8, height = 3)
