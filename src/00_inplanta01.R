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

#   Summary
cfu %>% 
    group_by(time_d) %>% 
    summarise(
        mean = mean(CFU_gFW),
        sd = sd(CFU_gFW))

#   Plot
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

#   Labels
label_trt_plant <- 
    plant_cell %>% 
    group_by(time_d, treatment) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, sep = "-"))

#   Threshold (T0)
summary_t0 <-
    plant_cell %>% 
    filter(time_d == 0 & treatment == "PFF2") %>% 
    summarise(
        n = n(),
        mean_rfu = mean(rfu),
        sd_rfu = sd(rfu),
        iqr_rfu = IQR(rfu),
        q3 = quantile(rfu, probs = 0.75),
        I = q3 + 1.5*iqr_rfu,
        mean_sd = mean_rfu + 2*sd_rfu,
        skewness = skewness(rfu),
        kurtosis = kurtosis(rfu)
    )

threshold <-
    summary_t0 %>% 
    pull(I)

#   Summary
plant_cell %>% 
    group_by(treatment, time_d) %>% 
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

#   Probability distribution
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
label_trt_plant2 <- 
    plant_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    tally %>% 
    mutate(id = paste(time_d, treatment, replicate, sep = "-"))

prob_plant_replicate <- 
    plant_cell %>% 
    group_by(time_d, treatment, replicate) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(label_trt_plant2$id) %>% 
    bind_rows(.id = "id") %>% 
    separate(id, into = c('time_d', 'trt', 'rep'), sep = "-", remove = FALSE)
prob_plant_replicate$time_d <- factor(prob_plant_replicate$time_d, levels = c("0", "2", "7"))
prob_plant_replicate$trt <- factor(prob_plant_replicate$trt, levels = c("inoculum", "PFF2"))

plt_prob_rep <- 
    prob_plant_replicate %>% 
    filter(ecdf < 1 & trt == "PFF2") %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, fill = time_d))+
    facet_grid(cols = vars(time_d))+
    geom_point(data = prob_plant_replicate %>% 
                   filter(time_d == 0 & ecdf < 1 & trt == "PFF2") %>% 
                   select(ecdf, x),
               fill = "gray", size = 1, alpha = 0.7, stroke = 0.1, pch = 21)+
    geom_point(size = 1, alpha = 0.8, stroke = 0.1, pch = 21)+
    scale_fill_manual(name = "Time [dpi]", 
                      values = c("grey", "#41AE78", "#00441B"), 
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

#   Median plot
#   Stat
summary_t0 <-
    prob_plant_replicate %>% 
    filter(trt == "PFF2" & time_d == 0 & ecdf < 1) %>% 
    summarise(
        n = n(),
        mean_rfu = mean(x),
        sd_rfu = sd(x),
        iqr_rfu = IQR(x),
        q3 = quantile(x, probs = 0.75),
        I = q3 + 1.5*iqr_rfu,
        mean_sd = mean_rfu + 2*sd_rfu,
        skewness = skewness(x),
        kurtosis = kurtosis(x)
    )

threshold <-
    summary_t0 %>% 
    pull(I)

#   Median and Interquantile range (IQR)
prob_plant_replicate %>% 
    filter(trt == "PFF2") %>% 
    group_by(time_d, rep) %>% 
    summarise(
        n = n(),
        median = median(x),
        q1 = quantile(x, probs = 0.25),
        q3 = quantile(x, probs = 0.75),
        IQR = paste0(sprintf("%.1f", q1), "--", sprintf("%.1f", q3)))

prob_plant_replicate %>% 
    filter(trt == "PFF2") %>% 
    mutate(rel_fl = x/threshold) %>% 
    group_by(time_d) %>% 
    summarise(
        n = n(),
        median = median(x),
        IQR = IQR(x),
        # Proportion of highly fluorescent cells
        prop_high = mean(x > threshold) * 100,
        rel_fl = mean(rel_fl),
        # Skewness and Kurtosis
        skewness = skewness(x),
        kurtosis = kurtosis(x),
        # Coefficient of Variation
        cv = (sd(x) / mean(x)) * 100)

#   Stats
inplanta_cell_sum <- 
    prob_plant_replicate %>% 
    filter(trt == "PFF2") %>% 
    group_by(time_d, rep) %>% 
    summarise(m_rfu = median(x),
              rel_rfu = m_rfu/threshold,
              .groups = "drop") %>% 
    mutate(time_d2 = case_when(
        time_d == "0" ~ 0,
        time_d == "2" ~ 2,
        time_d == "7" ~ 7))

#   Kruskal-Walis on treatments by day  
inplanta_cell_sum %>% 
    kruskal_test(rel_rfu ~ time_d2)

inplanta_pvalue <-
    inplanta_cell_sum %>%
    dunn_test(rel_rfu ~ time_d2, p.adjust.method = "BH") %>% 
    add_xy_position(step.increase = 0.05, stack = FALSE) %>% 
    mutate(xmin = as.numeric(group1),
           xmax = as.numeric(group2))

#   Plot
plot_sum_total <- 
    inplanta_cell_sum %>% 
    ggplot(aes(x = time_d2, y = rel_rfu))+
    geom_line(aes(group = time_d), linewidth = 0.5)+
    geom_jitter(aes(fill = time_d), width = 0.05,
                pch = 21, stroke = 0.5)+
    add_pvalue(data = inplanta_pvalue, 
               xmin = "xmin", xmax = "xmax",
               label.size = 3,
               tip.length = 0.01,
               bracket.size = 0.3)+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote("Relative"~italic("alkB")~"activity"~bgroup("(",over(RFU,RFU[t0]),")")),
                       limits = c(0, 4))+
    scale_fill_manual(name = "Time [dpi]", 
                      values = c("gray", "#41AE78", "#00441B"), 
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
    plot_layout(heights = c(1.2,1))

wrap_plots(plt_cfu, plt_prob_rep, plot_sum_total, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.2,1,1))

ggsave(here("results", "inplanta_exp1.png"), dpi = 300, width = 6, height = 8)
