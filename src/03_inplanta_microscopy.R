library(here)
library(tidyverse)
library(tidymodels)
library(moments)
library(rstatix)
library(ggprism)
library(patchwork)
set.seed(14071990)

#####   CFU     #####
##      Data
cfu_files <- 
    list.files(path = here('data'), 
               pattern = "\\w{3}\\d{1}_inplanta_cfu", 
               full.names = TRUE)

cfu_data <- 
    lapply(seq_along(cfu_files), function (i) {
        df <- read.csv(cfu_files[i]) %>% 
            select(replicate, time, CFU_gFW, logCFU) %>% 
            group_by(replicate, time) %>% 
            summarise(
                CFU_gFW = mean(CFU_gFW),
                logCFU = mean(logCFU), .groups = "drop") %>% 
            arrange(time)
        df$ID <- c("rep1", "rep2")[i]  # Add ID column
        return(df)
    })

cfu_df <-
    bind_rows(cfu_data)

cfu_df %>% 
    group_by(ID, time) %>% 
    summarise(median = median(CFU_gFW),
              q75 = quantile(CFU_gFW, probs = .75),
              iqr = IQR(CFU_gFW),
              I = q75 + 1.5 * iqr,
              .groups = "drop")

##      Stats
cfu_df %>% 
    group_by(ID) %>% 
    wilcox_test(logCFU ~ time, p.adjust.method = "holm")

##      Plots
plt_cfu <- 
    cfu_df %>% 
    ggplot(aes(x = time, y = logCFU, group = ID))+
    geom_point(alpha = 0.5, position = position_jitter(width = 0.02))+
    stat_summary(fun = mean, geom = "line", linetype = "dashed")+
    stat_summary(aes(shape = ID),
                 fun.y = mean, 
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x),
                 geom = "pointrange", fill = "white", color = "black", linewidth = 0.8)+
    scale_y_continuous(name = bquote("Bacterial load ["~log[10]~CFU~gFW^-1~"]"), limits = c(4, 8))+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 7.2), breaks = c(0, 2, 7))+
    scale_shape_manual(name = "Experiment", labels = c("Exp1", "Exp2"), values = c(21, 23))+
    scale_colour_manual(values = c("black", "grey"))+
    guides(color = "none",
           shape  = guide_legend(position = "inside"))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position.inside = c(0.8, 0.2),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"))

#####   Microscopy     #####
##  Load model
cell_lr_fit <- readRDS(here('results', 'lr_model.rds'))

##  Load data
single_cell_files <- 
    list.files(path = here('data'), 
               pattern = "\\w{3}\\d{1}_single_cell", 
               full.names = TRUE)

single_cell_data <- 
    lapply(seq_along(single_cell_files), function (i) {
        df <- read.csv(single_cell_files[i]) %>% 
            mutate(Label = str_replace(Label, ".czi:.*", "")) %>%
            separate(Label, into = c("Channel", "Time", "Rep", "Img"), sep = "-") %>% 
            select(-X, -Channel)
        df$ID <- c("rep1", "rep2")[i]  # Add ID column
        return(df)
    })

single_cell_df <-
    bind_rows(single_cell_data)

##  Prediction
cell_data <- 
    cell_lr_fit %>% 
    extract_workflow() %>% 
    predict(new_data = single_cell_df) %>% 
    cbind(single_cell_df, .) %>% 
    filter(.pred_class == "Yes")

cell_data %>% 
    group_by(ID, Time) %>% 
    tally()

##  Background
background_files <- 
    list.files(path = here('data'), 
               pattern = "\\w{3}\\d{1}_background", 
               full.names = TRUE)

background_data <- 
    lapply(seq_along(background_files), function (i) {
        df <- read.csv(background_files[i]) %>% 
            mutate(Label = str_replace(Label, ".czi:.*", "")) %>%
            separate(Label, into = c("Time", "Rep", "Img"), sep = "-") %>% 
            select(-X)
        df$ID <- c("rep1", "rep2")[i]  # Add ID column
        return(df)
    })

background_df <- 
    background_data %>% 
    bind_rows() %>% 
    group_by(ID, Time, Rep, Img) %>% 
    summarise(median = median(Mean),
              q75 = quantile(Mean, probs = .75),
              iqr = IQR(Mean),
              I = q75 + 1.5 * iqr,
              .groups = "drop")

##  Exploratory
background_data %>% 
    bind_rows() %>% 
    ggplot(aes(x = Mean))+
    facet_grid(cols = vars(Time), rows = vars(ID))+
    geom_histogram()


##  Background-corrected data
sc_fl_df <- 
    cell_data %>% 
    select(ID, Time, Rep, Img, Area, Mean) %>% 
    inner_join(., background_df, 
               by = c("ID", "Time", "Rep", "Img")) %>% 
    mutate(rfu = Mean - median) %>% 
    filter(rfu > 0)
sc_fl_df <- sc_fl_df[-which(sc_fl_df$ID == "rep1" & sc_fl_df$Time == "T0" & sc_fl_df$rfu > 25), ]

sc_fl_df %>% 
    group_by(ID, Time, Rep) %>% 
    tally() %>% 
    group_by(ID, Time) %>% 
    tally()

sc_fl_df %>% 
    ggplot(aes(sample = rfu, group = Rep))+
    facet_grid(cols = vars(Time), rows = vars(ID))+
    geom_qq()

#   Normal probability plots
##  Labels
exp_lab <- c("rep1", "rep2")
time_lab <- c("T0", "T1", "T2")
label_rep_time <-
    expand.grid(exp_lab, time_lab) %>% 
    mutate(label = paste(Var1, Var2, sep = "_"))
    
sc_fl_label <- 
    sc_fl_df %>% 
    group_by(ID, Time, Rep) %>% 
    tally() %>% 
    mutate(label = paste(ID, Time, Rep, sep = "_"))

##  Total
prob_plant_total <- 
    sc_fl_df %>% 
    group_by(ID, Time) %>% 
    mutate(ecdf_prob = ecdf(rfu)(rfu),
           time = case_when(
               Time == "T0" ~ 0,
               Time == "T1" ~ 2,
               Time == "T2" ~ 7)) %>% 
    ungroup()

### Plot
prob_plant_total %>% 
    filter(ecdf_prob < 1) %>% 
    ggplot(aes(x = qnorm(ecdf_prob), y = rfu, colour = Time))+
    facet_grid(cols = vars(ID))+
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
prob_plant_total

##  Per replicate
prob_plant_rep <- 
    sc_fl_df %>% 
    group_by(ID, Time, Rep) %>% 
    mutate(ecdf_prob = ecdf(rfu)(rfu),
           time = case_when(
               Time == "T0" ~ 0,
               Time == "T1" ~ 2,
               Time == "T2" ~ 7)) %>% 
    ungroup()

### Plot
lab_Time = c(T0 = "0 dpi", T1 = "2 dpi", T2 = "7 dpi")
lab_ID = c(rep1 = "Exp 1", rep2 = "Exp 2")

plt_prob <- 
    prob_plant_rep %>% 
    filter(ecdf_prob < 1) %>% 
    ggplot(aes(x = qnorm(ecdf_prob), y = rfu, group = ID))+
    facet_grid(cols = vars(Time), rows = vars(ID), 
               labeller = labeller(.cols = lab_Time, .rows = lab_ID))+
    geom_vline(xintercept = qnorm(0.9), linetype = "dotdash", alpha = 0.5, linewidth = 0.3)+
    geom_point(data = prob_plant_rep %>% 
                   filter(time == 0 & ecdf_prob < 1) %>% 
                   select(ID, ecdf_prob, rfu),
               fill = "gray", size = 1, alpha = 0.5, stroke = 0.1, pch = 21)+
    geom_point(aes(fill = Time, group = Rep), size = 1, alpha = 0.5, stroke = 0.1, pch = 21)+
    scale_fill_manual(name = "Time [dpi]", 
                       values = c("gray", "#41AE78", "#00441B"), 
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
        strip.background = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

#   Median plot
#   Stat
summary_t0 <-
    prob_plant_rep %>% 
    group_by(ID) %>% 
    filter(time == 0 & ecdf_prob < 1) %>% 
    summarise(
        n = n(),
        mean_rfu = mean(rfu),
        sd_rfu = sd(rfu),
        iqr_rfu = IQR(rfu),
        q3 = quantile(rfu, probs = 0.75),
        I_thr = q3 + 1.5 * iqr_rfu,
        mean_sd = mean_rfu + 2*sd_rfu,
        skewness = skewness(rfu),
        kurtosis = kurtosis(rfu)
    )

threshold <-
    summary_t0 %>% 
    group_by(ID) %>% 
    select(ID, I_thr)

#   Median and Interquantile range (IQR)
prob_plant_rep %>% 
    group_by(ID, Time, Rep) %>% 
    summarise(
        n = n(),
        median = median(rfu),
        q1 = quantile(rfu, probs = 0.25),
        q3 = quantile(rfu, probs = 0.75),
        iqr = IQR(rfu),
        iqr_range = paste0(sprintf("%.1f", q1), "-", sprintf("%.1f", q3)))

prob_plant_rep %>% 
    group_by(ID, Rep, Time) %>% 
    tally

prob_plant_rep %>% 
    inner_join(., threshold, by = "ID") %>% 
    mutate(rel_fl = rfu/I_thr) %>% 
    group_by(ID, Time) %>% 
    summarise(
        n = n(),
        median = median(rfu),
        IQR = IQR(rfu),
        # Proportion of highly fluorescent cells
        prop_high = mean(rfu > I_thr) * 100,
        rel_fl = median(rel_fl),
        # Skewness and Kurtosis
        skewness = skewness(rfu),
        kurtosis = kurtosis(rfu))

#   Stats
inplanta_cell_sum <- 
    prob_plant_rep %>%
    group_by(ID, Time, time, Rep) %>% 
    summarise(m_rfu = median(rfu), .groups = "drop") %>% 
    inner_join(., threshold, by = "ID") %>% 
    mutate(rel_rfu = m_rfu / I_thr)

#   Kruskal-Walis on treatments by day  
inplanta_cell_sum %>% 
    group_by(ID) %>% 
    kruskal_test(rel_rfu ~ time)

inplanta_pvalue <-
    inplanta_cell_sum %>%
    group_by(ID) %>% 
    dunn_test(rel_rfu ~ time, p.adjust.method = "holm") %>% 
    add_xy_position(x = "ID", step.increase = 0.2, stack = FALSE) %>% 
    mutate(xmin = as.numeric(group1),
           xmax = as.numeric(group2)) %>% 
    filter(p.adj < 0.1)

#   Plot
plot_sum_total <- 
    inplanta_cell_sum %>% 
    ggplot(aes(x = time, y = rel_rfu))+
    facet_grid(rows = vars(ID), labeller = labeller(.rows = lab_ID))+
    geom_line(aes(group = Time), linewidth = 0.5)+
    geom_jitter(aes(fill = Time), width = 0.05,
               pch = 21, stroke = 0.5)+
    add_pvalue(data = inplanta_pvalue, 
               xmin = "xmin", xmax = "xmax",
               label.size = 4,
               tip.length = 0.01,
               bracket.size = 0.3)+
    scale_x_continuous(name = "Time [dpi]", 
                       limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote("Relative"~italic("alkB")~"activity"~bgroup("(",over(RFU,RFU[t0]),")")),
                       limits = c(0, 1.1), breaks = seq(0, 1, 0.2))+
    scale_fill_manual(name = "Time [dpi]", 
                      values = c("gray", "#41AE78", "#00441B"), 
                      labels = c("0", "2", "7"))+
    labs(title = "100% of cells")+
    guides(fill = "none")+
    theme(
        title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

##  Top 10% of cells
inplanta_cell_sum_top <- 
    prob_plant_rep %>% 
    filter(ecdf_prob > 0.9) %>% 
    group_by(ID, time, Time, Rep) %>% 
    summarise(m_rfu = median(rfu), .groups = "drop") %>% 
    inner_join(., threshold, by = "ID") %>% 
    mutate(rel_rfu = m_rfu / I_thr)

#   Kruskal-Walis on treatments by day  
inplanta_cell_sum_top %>% 
    group_by(ID) %>% 
    kruskal_test(rel_rfu ~ time)

inplanta_pvalue_top <-
    inplanta_cell_sum_top %>%
    group_by(ID) %>% 
    dunn_test(rel_rfu ~ time, p.adjust.method = "holm") %>% 
    add_xy_position(x = "ID", step.increase = 0.1, stack = FALSE) %>% 
    mutate(xmin = as.numeric(group1),
           xmax = as.numeric(group2)) %>% 
    filter(p.adj < 0.05)
    
plot_sum_top <- 
    inplanta_cell_sum_top %>%
    ggplot(aes(x = time, y = rel_rfu))+
    facet_grid(rows = vars(ID), labeller = labeller(.rows = lab_ID))+
    geom_hline(yintercept = 1, linetype = "dotdash", alpha = 0.5, linewidth = 0.3)+
    geom_line(aes(group = Time), linewidth = 0.5)+
    geom_jitter(aes(fill = Time), width = 0.05,
                pch = 21, stroke = 0.5)+
    add_pvalue(data = inplanta_pvalue_top, 
               xmin = "xmin", xmax = "xmax",
               label.size = 4,
               tip.length = 0.01,
               bracket.size = 0.3)+
    scale_x_continuous(name = "Time [dpi]", limits = c(-0.2, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote("Relative"~italic("alkB")~"activity"~bgroup("(",over(RFU,RFU[t0]),")")),
                       limits = c(0, 6.5))+
    scale_fill_manual(name = "Time [dpi]", 
                      values = c("gray", "#41AE78", "#00441B"), 
                      labels = c("0", "2", "7"))+
    guides(fill = "none")+
    theme(axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25))

## Wrap plots
p <- plot_sum_total + plot_sum_top + plot_layout(widths = c(1,1))

wrap_plots(plt_cfu, plt_prob_2, p, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.2,1,1))

p2 <- plt_prob_2 + plot_sum_top +
    plot_layout(widths = c(3,1))

wrap_plots(plt_cfu, p2, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1,1.2))

wrap_plots(plt_cfu, plt_prob, plot_sum_top, ncol = 3)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(widths = c(1, 3, 1))

ggsave(here("results", "inplanta_exp2.png"), dpi = 300, width = 6, height = 8)
