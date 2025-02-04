library(here)
library(tidyverse)
library(tidymodels)
library(rstatix)
library(ggpubr)
set.seed(14071990)

#####   CFU #####
##      Data
cfu <- read.csv(here('data', 'inplanta_cfu_e2.csv'))
str(cfu)

##      Stats
cfu %>% 
    wilcox_test(log_CFU_g ~ time, p.adjust.method = "BH")

##      Plots
plt_cfu_2 <- 
    cfu %>% 
    ggplot(aes(x = time, y = log_CFU_g))+
    geom_point(alpha = 0.5, position = position_jitter(width = 0.02))+
    stat_summary(fun = mean, geom = "line", linetype = "dashed")+
    stat_summary(fun.y = mean, 
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x),
                 shape = 21, geom = "pointrange", fill = "white", color = "black", linewidth = 0.8)+
    scale_y_continuous(name = bquote("Bacterial load ["~log[10]~CFU~gFW^-1~"]"), limits = c(4, 8))+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 7.2), breaks = c(0, 2, 7))+
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

## Load model
cell_lr_fit <- readRDS(here('results', 'lr_model.rds'))

## Load data
data <- read.csv(here('data', 'singlecell_data.csv'))
data_select <- 
    data %>% 
    select(-X)

background <- read.csv(here('data', 'inplanta_background.csv'))

bg_df <- background %>% 
    na.omit() %>% 
    select(Label, Mean) %>% 
    mutate(Label = str_replace(Label, ".czi:.*", "")) %>% 
    separate(Label, into = c("Time", "Rep", "Img"), sep = "-")

bg_values <- bg_df %>%
    group_by(Time, Rep, Img) %>% 
    summarise(median = median(Mean),
              iqr = quantile(Mean, probs = 0.25))

##  Prediction
cell_data <- 
    cell_lr_fit %>% 
    extract_workflow() %>% 
    predict(new_data = data_select) %>% 
    cbind(data, .) %>% 
    filter(.pred_class == "Yes")

#   Clean data
sc_fl_df <- 
    cell_data %>% 
    mutate(Label = str_replace(Label, ".czi:.*", "")) %>% 
    separate(Label, into = c("Channel", "Time", "Rep", "Img"), sep = "-") %>% 
    select(Time, Rep, Img, Area, Mean) %>% 
    inner_join(., bg_values, by = c("Time", "Rep", "Img")) %>% 
    mutate(rfu = Mean - iqr) %>% 
    filter(rfu > 0)

sc_fl_df %>% 
    group_by(Time, Rep, Img) %>% 
    tally

sc_fl_df %>% 
    ggplot(aes(sample = rfu, group = Rep))+
    facet_grid(cols = vars(Time))+
    geom_qq()

#   Normal probability plots
##  Label
sc_fl_label <- sc_fl_df %>% 
    group_by(Time, Rep) %>% 
    tally() %>% 
    mutate(label = paste(Time, Rep, sep = "_"))

##  Total
prob_plant_total <- sc_fl_df %>% 
    group_by(Time) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(nm = c("T0", "T1", "T2")) %>% 
    bind_rows(.id = "Time")
prob_plant_total$Time <- factor(prob_plant_total$Time, 
                          levels = c("T0", "T1", "T2"),
                          labels = c(0, 2, 7))

### Plot
prob_plant_total %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, colour = Time))+
    geom_point(size = 1, alpha = 0.75, stroke = 0)+
    #scale_color_manual(name = "Time [dpi]", values = palette_d_plant, labels = c("0 (Inoculum)", "0", "2", "7"))+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    #guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d_plant, stroke = 0)))+
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

##  Per replicate
prob_plant_rep <- sc_fl_df %>% 
    group_by(Time, Rep) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(sc_fl_label$label) %>% 
    bind_rows(.id = "label") %>% 
    separate(label, into = c('Time', 'Rep'), sep = "_", remove = FALSE)
prob_plant_rep$Time <- factor(prob_plant_rep$Time, 
                          levels = c("T0", "T1", "T2"),
                          labels = c(0, 2, 7))

### Plot
plt_prob_2 <- prob_plant_rep %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, fill = Time, group = Rep))+
    facet_grid(cols = vars(Time))+
    geom_point(size = 1, alpha = 0.5, stroke = 0.1, pch = 21)+
    scale_fill_manual(name = "Time [dpi]", 
                       values = c("#FED976", "#41AE76", "#00441B"), 
                       labels = c("0", "2", "7"))+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4,4))+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    guides(fill = "none")+
    theme(
        axis.text = element_text(size = 5, color = "black"),
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
prob_plant_rep_df <- prob_plant_rep %>% 
    filter(ecdf > 0.9) %>% 
    group_by(Time, Rep) %>%
    summarise(logrfu = median(log(x)), .groups = "drop") %>% 
    mutate(time_d = case_when(
        Time == "0" ~ 0,
        Time == "2" ~ 2,
        Time == "7" ~ 7))

prob_plant_rep_df %>% 
    wilcox_test(logrfu ~ Time, p.adjust.method = "BH")

label_stats <- prob_plant_rep_df %>% 
    group_by(time_d) %>% 
    tally %>% 
    mutate(y.position = 5.2,
           grp = c("a", "b", "b"))

#   Plot
plt_alkb_2 <- prob_plant_rep_df %>% 
    ggplot(aes(x = time_d, y = logrfu))+
    geom_boxplot(aes(fill = Time), color = "black", width = 0.4, linewidth = 0.2, outlier.alpha = 0, alpha = 0.5)+
    geom_point(position = position_jitter(width = 0.1), alpha = 0.8, size = 1)+
    geom_text(data = label_stats, aes(y = y.position, label = grp), size = 3)+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote(italic("alkB")~"activity ["~log[10]~RFU~"]"),
                       limits = c(3, 5.2))+
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
p <- plt_prob_2 + plt_alkb_2 + plot_layout(widths = c(3,1))

wrap_plots(plt_cfu_2, p, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.5,1))
