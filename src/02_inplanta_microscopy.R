library(here)
library(tidyverse)
library(tidymodels)
library(moments)
library(rstatix)
library(ggprism)
library(patchwork)
set.seed(14071990)

#####   Microscopy     #####
##  Load model
cell_rf_fit <- readRDS(here('results', 'inplanta_rf.rds'))

##  Load data
single_cell_files <- 
    list.files(path = here('data'), 
               pattern = "\\w{3}\\d{1}_inplanta_single_cell", 
               full.names = TRUE)

single_cell_data <- 
    lapply(seq_along(single_cell_files), function (i) {
        df <- read.csv(single_cell_files[i]) %>% 
            mutate(Label = str_replace(Label, ".czi:.*", "")) %>%
            mutate(Label = str_replace_all(Label, "PFF2_", "")) %>% 
            mutate(Label = str_replace_all(Label, "_", "-")) %>%
            separate(Label, into = c("Channel", "Time", "Rep", "Img"), sep = "-") %>% 
            select(-X, -Channel)
        df$ID <- c("rep1", "rep2")[i]  # Add ID column
        return(df)
    })

single_cell_df <-
    bind_rows(single_cell_data) %>% 
    na.exclude %>% 
    mutate(Time = case_when(
        Time == "Day0" ~ "T0",
        Time == "Day2" ~ "T1",
        Time == "Day7" ~ "T2",
        TRUE ~ Time))

##  Prediction
cell_data <- 
    cell_rf_fit %>% 
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
               pattern = "\\w{3}\\d{1}_inplanta_background", 
               full.names = TRUE)

background_data <- 
    lapply(seq_along(background_files), function (i) {
        df <- read.csv(background_files[i]) %>% 
            mutate(Label = str_replace(Label, ".czi:.*", "")) %>%
            mutate(Label = str_replace_all(Label, "PFF2_", "")) %>% 
            mutate(Label = str_replace_all(Label, "_", "-")) %>%
            separate(Label, into = c("Time", "Rep", "Img"), sep = "-") %>% 
            mutate(Time = case_when(
                Time == "Day0" ~ "T0",
                Time == "Day2" ~ "T1",
                Time == "Day7" ~ "T2",
                TRUE ~ Time)) %>% 
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
inplanta_fluorescence_df <- 
    cell_data %>% 
    select(ID, Time, Rep, Img, Area, Mean) %>% 
    inner_join(., background_df, 
               by = c("ID", "Time", "Rep", "Img")) %>% 
    mutate(rfu = Mean - median) %>% 
    filter(rfu > 0) %>% 
    mutate(time = case_when(
        Time == "T0" ~ 0,
        Time == "T1" ~ 2,
        Time == "T2" ~ 7))
#sc_fl_df <- sc_fl_df[-which(sc_fl_df$ID == "rep1" & sc_fl_df$Time == "T0" & sc_fl_df$rfu > 25), ]

##      Exploratory analysis of background corrected data
#       Summary
inplanta_fluorescence_df %>% 
    group_by(ID, Time) %>% 
    summarise(
        median_fl = median(rfu),
        skewness = skewness(rfu),
        kurtosis = kurtosis(rfu)) %>% 
    arrange(Time)
    
#       Number of replicates
inplanta_fluorescence_df %>% 
    group_by(ID, Time, Rep) %>% 
    tally() %>% 
    group_by(ID, Time) %>% 
    tally()

#       Q-Q plot
inplanta_fluorescence_df %>% 
    ggplot(aes(sample = rfu, group = Rep))+
    facet_grid(cols = vars(Time), rows = vars(ID))+
    geom_qq()

inplanta_fluorescence_df %>% 
    ggplot(aes(x = rfu, group = Rep))+
    facet_grid(cols = vars(Time), rows = vars(ID))+
    geom_histogram()

#####   Normal probability plots    #####
##      Labels
exp_lab <- 
    c("rep1" = "Rep 1",
      "rep2" = "Rep 2")
time_lab <- 
    c("T0" = "0 dpi",
      "T1" = "2 dpi",
      "T2" = "7 dpi")

##      Total Cumulative Distribution
prob_inplanta_total <- 
    inplanta_fluorescence_df %>% 
    group_by(ID, Time, time, Rep) %>% 
    #   Cumulative distribution
    mutate(ecdf_prob = ecdf(rfu)(rfu)) %>% 
    ungroup()

### Plot
plot_inplanta_distribution<-
    prob_inplanta_total %>% 
    filter(ecdf_prob < 1) %>% 
    ggplot(aes(x = qnorm(ecdf_prob), y = rfu, group = ID))+
    facet_grid(cols = vars(Time), rows = vars(ID), 
               labeller = labeller(.cols = lab_Time, .rows = exp_lab))+
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

###     Normalisation by control
##      T0 control
inplanta_control <-
    inplanta_fluorescence_df %>% 
    filter(Time == "T0") %>%    # Filter BHB no diesel
    group_by(ID) %>% 
    summarise(
        Q75_ctrl = quantile(rfu, probs = 0.75),  # 75th percentile
        Q99_ctrl = quantile(rfu, probs = 0.99),   # 99th percentile
        I = Q75_ctrl + 1.5 * IQR(rfu),
        .groups = "drop")

##      Normalised single cell data by control
inplanta_normalised_df <- 
    inplanta_fluorescence_df %>%
    left_join(., inplanta_control, 
              by = "ID") %>% 
    group_by(Time, time, Rep) %>% 
    mutate(relative_rfu = rfu / I.y,
           ecdf_prob = ecdf(relative_rfu)(relative_rfu))

##  Descriptive statistics of relative fluorescence
inplanta_normalised_df %>% 
    group_by(ID, Time) %>% 
    tally()

inplanta_normalised_df %>% 
    group_by(ID, Time, time, Rep) %>% 
    summarise(
        n = n(),
        relative_fl = mean(relative_rfu),
        prop_high = mean(rfu > I.y) * 100,
        skewness = skewness(relative_rfu),
        kurtosis = kurtosis(relative_rfu),
        .groups = "drop") %>% 
    group_by(ID, time) %>% 
    summarise(
        median_rel_fl = median(relative_fl),
        IQR_rel_fl = IQR(relative_fl),
        median_prop_high = median(prop_high),
        median_skeweness = median(skewness),
        median_kurtosis = median(kurtosis))

##      Plot normalised probability distribution
plot_inplanta_normalised_distribution <-
    inplanta_normalised_df %>% 
    filter(ecdf_prob < 1) %>% 
    ggplot(aes(x = qnorm(ecdf_prob), y = relative_rfu))+
    facet_grid(rows = vars(ID),
               cols = vars(Time), 
               labeller = labeller(ID = exp_lab, Time = time_lab))+
    geom_vline(aes(xintercept = qnorm(0.90)), linetype = "dotdash", alpha = 0.5, linewidth = 0.3)+
    geom_point(aes(colour = as.factor(time)),
               size = 1.5, alpha = 0.75, stroke = 0)+
    scale_colour_manual(name = "Time [dpi]", 
                      values = c("gray", "#41AE78", "#00441B"), 
                      labels = c("0", "2", "7"))+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4.2,4.2))+
    scale_y_continuous(name = "Relative fluorescence")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, stroke = 0)))+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        strip.background = element_blank())

##      Top 10% fluorescence values per group
inplanta_top_percent <- 
    inplanta_normalised_df %>%
    group_by(ID, time, Rep) %>%
    filter(ecdf_prob > 0.90) %>%  # Select top 10%
    summarise(
        n = n(),
        median_top = median(relative_rfu), 
        mean_top = mean(relative_rfu), 
        q99 = quantile(relative_rfu, 0.99), 
        q95 = quantile(relative_rfu, 0.95), 
        .groups = "drop") %>%
    mutate(Time = factor(time))

#       Test for normality
inplanta_top_percent %>% 
    group_by(ID, Time) %>% 
    shapiro_test(median_top)

#       Test for homoscedasticity
inplanta_top_percent %>% 
    group_by(ID) %>% 
    levene_test(median_top ~ Time)

##      Are they different to zero?
inplanta_top_percent %>% 
    group_by(ID, Time) %>% 
    t_test(median_top ~ 0, mu = 1, var.equal = TRUE)

##      Are they different between sampling time
inplanta_top_percent %>% 
    group_by(ID) %>% 
    kruskal_test(median_top ~ Time)

inplanta_p_value_top <- 
    inplanta_top_percent %>% 
    group_by(ID) %>% 
    dunn_test(median_top ~ Time, p.adjust.method = "holm") %>% 
    add_xy_position() %>% 
    mutate(xmin = as.numeric(group1),
           xmax = as.numeric(group2)) %>% 
    filter(p.adj < 0.05)

#       Plot 
plot_inplanta_top_percent<-
    inplanta_top_percent %>% 
    ggplot(aes(x = time, y = median_top))+
    facet_grid(rows = vars(ID),
               labeller = labeller(ID = exp_lab))+
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, alpha = 0.2)+
    geom_line(aes(group = time),
              linewidth = 0.5)+
    geom_point(position = position_jitter(width = 0.1),
               alpha = 0.8)+
    add_pvalue(data = inplanta_p_value_top, 
               xmin = "xmin", xmax = "xmax",
               label.size = 4, tip.length = 0.01,
               bracket.size = 0.3, step.increase = 0.001)+
    scale_x_continuous(name = "Time [d]", limits = c(-0.5, 8), breaks = c(0, 2, 7))+
    scale_y_continuous(name = bquote(atop("Fold change fluorescence", "(" * Q[90] * ", normalised by T0)")),
                       limits = c(0, 7))+
    guides(fill = "none")+
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())

##  Plots
plot_inplanta_distribution
ggsave(here("results", "inplanta_distribution.png"), dpi = 300, width = 6, height = 4)

plot_inplanta_normalised_distribution + plot_inplanta_top_percent +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(3, 1))
ggsave(here("results", "inplanta_normalised.png"), dpi = 300, width = 10, height = 4)