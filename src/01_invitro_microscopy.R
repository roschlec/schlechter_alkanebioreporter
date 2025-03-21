library(here)
library(tidyverse)
library(tidymodels)
library(moments)
library(rstatix)
library(ggprism)
library(patchwork)
set.seed(14071990)

#####   Single-cell data clean up     #####
##      Load random forest model
cell_rf_fit <- readRDS(here('results', 'invitro_rf.rds'))

##      Load data
single_cell_data <- 
    read.csv(here('data', 'invitro_single_cell.csv')) %>% 
    mutate(Label = str_replace(Label, ".czi:.*", ""),
           Label = str_replace_all(Label, "-", "_")) %>% 
    separate(Label, into = c("Channel", "Time", "Strain", "Treatment", "Rep", "Img"), sep = "_") %>% 
    select(-X, -Channel, -Strain)

##      Classify cells
cell_data <- 
    cell_rf_fit %>% 
    extract_workflow() %>% 
    predict(new_data = single_cell_data) %>% 
    cbind(single_cell_data, .) %>% 
    filter(.pred_class == "Yes") %>% 
    mutate(Treatment = str_replace(Treatment, "nodiesel", "no_diesel"))

cell_data %>% 
    group_by(Treatment, Time, Rep) %>% 
    tally() # Number of identified cells per treatment and time point

##      Background data
background_data <- 
    read.csv(here('data', 'invitro_background.csv')) %>% 
    mutate(Label = str_replace(Label, ".czi:.*", ""),
           Label = str_replace(Label, "nodiesel", "no_diesel")) %>% 
    separate(Label, into = c("Time", "Strain", "Treatment", "Rep", "Img"), sep = "-") %>% 
    select(-X, -Strain)

##      Summary statistics of background data
background_df <- 
    background_data %>% 
    group_by(Time, Treatment, Rep, Img) %>% 
    summarise(
        mode = mode(Mean),
        median = median(Mean),
        q75 = quantile(Mean, probs = .75),
        iqr = IQR(Mean),
        I = q75 + 1.5 * iqr,
        .groups = "drop")

##      Exploratory analysis of background data
background_data %>% 
    ggplot(aes(x = Mean))+
    facet_grid(cols = vars(Time), rows = vars(Treatment))+
    geom_histogram(bins = 50)

##  Background-corrected data
invitro_fluorescence_df <- 
    cell_data %>% 
    select(Time, Treatment, Rep, Img, Mean) %>% 
    inner_join(., background_df, 
               by = c("Time", "Treatment", "Rep", "Img")) %>% 
    #   Correct fluorescence signal by median background fluorescence per image
    mutate(rfu = Mean - median) %>%
    filter(rfu >= 0) %>% 
    mutate(time = case_when(
        Time == "20220902" ~ 1,
        Time == "20220907" ~ 6,
        Time == "20220928" ~ 28))

#       One replicate only has four cells, so it is excluded from the analysis
invitro_fluorescence_df <-
    invitro_fluorescence_df[-which(
        invitro_fluorescence_df$Treatment == "no_diesel" & 
            invitro_fluorescence_df$Time == "20220928" & 
            invitro_fluorescence_df$Rep == "R3"), ]

##      Exploratory analysis of background corrected data
#       Number of replicates
invitro_fluorescence_df %>% 
    group_by(time, Treatment, Rep) %>% 
    tally() %>% 
    group_by(time, Treatment) %>% 
    tally()

#       Q-Q plot
invitro_fluorescence_df %>% 
    ggplot(aes(sample = rfu))+
    facet_grid(cols = vars(time), rows = vars(Treatment))+
    geom_qq()
    
#####   Normal probability plots    #####
##      Labels
treatment_lab <- 
    c("diesel" = "BHB + 1% v/v Diesel", 
      "no_diesel" = "BHB",
      "lb" = "LB")

time_lab <- 
    c("20220902" = "1 dpi",
      "20220907" = "6 dpi",
      "20220928" = "28 dpi")

##      Total Cumulative Distribution
prob_invitro_total <- 
    invitro_fluorescence_df %>% 
    group_by(Treatment, Time, time, Rep) %>% 
    #   Cumulative distribution
    mutate(ecdf_prob = ecdf(rfu)(rfu)) %>% 
    ungroup()

##  Summary statistics
##  Proportion of non-fluorescent cells
prob_invitro_total %>% 
    group_by(Treatment, time, Rep) %>% 
    filter(rfu < 256*0.1) %>%   # filter by 10% of max fluorescence
    summarise(prob = max(ecdf_prob)) %>% 
    group_by(Treatment, time) %>% 
    summarise(median_prob = median(prob))

prob_invitro_total %>% 
    filter(Treatment != "lb") %>% 
    group_by(Treatment, time, Rep) %>% 
    summarise(
        skewness = skewness(rfu),
        kurtosis = kurtosis(rfu)) %>% 
    group_by(Treatment) %>% 
    summarise(
        skewness = mean(skewness),
        kurtosis = mean(kurtosis))

##      Cumulative Distribution of BHB-grown cells
prob_invitro_bhb <-
    prob_invitro_total %>% 
    filter(Treatment != "lb")   # Exclude LB data

##      Cumulative Distribution of LB-grown cells
prob_invitro_lb <-
    invitro_fluorescence_df %>% 
    filter(Treatment == "lb") %>%    # Filters LB data
    mutate(ecdf_prob = ecdf(rfu)(rfu)) %>% 
    select(Time, rfu, ecdf_prob, Treatment)

prob_invitro_lb %>% 
    ggplot(aes(x = ecdf_prob, y = rfu))+
    geom_point()
    
###     Probability distribution plot
plot_invitro_distribution <-
    prob_invitro_bhb %>% 
    filter(ecdf_prob < 1) %>% 
    #   Probability distribution
    ggplot(aes(x = qnorm(ecdf_prob), y = rfu))+
    #   Separate by sampling point and BHB growth condition
    facet_grid(cols = vars(Time),
               rows = vars(Treatment),
               labeller = labeller(
                   Time = time_lab,
                   Treatment = treatment_lab))+
    geom_vline(aes(xintercept = qnorm(0.99)), linetype = "dotdash", alpha = 0.5, linewidth = 0.3)+
    #   Include LB data
    geom_point(data = prob_invitro_lb %>% 
                   filter(ecdf_prob < 1) %>%
                   mutate(Control = Treatment) %>% 
                   select(ecdf_prob, rfu, Control),
               aes(shape = Control),
               size = 1.5, alpha = 0.5, stroke = 0, color = "darkgray")+
    #   Probability distribution
    geom_point(aes(colour = as.factor(time)),
               size = 1.5, alpha = 0.75, stroke = 0)+
    #   Scales
    scale_color_manual(name = "Time [dpi]", values = palette_d)+
    scale_shape_manual(name = "Control", values = 18, label = "LB")+
    scale_y_continuous(name = "Single-cell fluorescence [a.u.]")+
    scale_x_continuous(name = "Cumulative distribution [%]", 
                       breaks = qnorm(b), 
                       labels = scales::percent(b, suffix = ""), 
                       limits = c(-4.2,4.2))+
    #   Guides
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d, stroke = 0)),
           shape = guide_legend(override.aes = list(alpha = 1, size = 4, stroke = 0)))+
    #   Themes
    theme(
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25),
        strip.background = element_blank())

###     Normalisation by control
##      BHB control
control <-
    invitro_fluorescence_df %>% 
    filter(Treatment == "no_diesel") %>%    # Filter BHB no diesel
    group_by(Time, time) %>% 
    summarise(
        Q75_ctrl = quantile(rfu, probs = 0.75),  # 75th percentile
        Q99_ctrl = quantile(rfu, probs = 0.99),   # 99th percentile
        I = Q75_ctrl + 1.5 * IQR(rfu),
        .groups = "drop")

##      Normalised single cell data by control
normalised_df <- 
    invitro_fluorescence_df %>% 
    filter(Treatment == "diesel") %>% 
    left_join(., control, 
              by = c("Time", "time")) %>% 
    group_by(Time, time, Rep) %>% 
    mutate(relative_rfu = rfu / I.y,
           ecdf_prob = ecdf(relative_rfu)(relative_rfu))

##  Descriptive statistics of relative fluorescence
normalised_df %>% 
    group_by(Time) %>% 
    tally()

normalised_df %>% 
    group_by(Time, Rep) %>% 
    filter(ecdf_prob == 1) %>% 
    select(time, Treatment, relative_rfu) %>% 
    group_by(Time) %>% 
    summarise(mean_rel = mean(relative_rfu))

normalised_df %>% 
    group_by(Time, time, Rep) %>% 
    summarise(
        n = n(),
        relative_fl = median(relative_rfu),
        prop_high = mean(rfu > I.y) * 100,
        skewness = skewness(relative_rfu),
        kurtosis = kurtosis(relative_rfu),
        .groups = "drop") %>% 
    group_by(Time) %>% 
    summarise(
        median_rel_fl = median(relative_fl),
        IQR_rel_fl = IQR(relative_fl),
        median_prop_high = median(prop_high),
        median_skeweness = median(skewness),
        median_kurtosis = median(kurtosis))

##      Plot normalised probability distribution
plot_invitro_normalised_distribution <-
    normalised_df %>% 
    filter(ecdf_prob < 1) %>% 
    ggplot(aes(x = qnorm(ecdf_prob), y = relative_rfu))+
    facet_grid(cols = vars(Time), labeller = labeller(Time = time_lab))+
    geom_vline(aes(xintercept = qnorm(0.99)), linetype = "dotdash", alpha = 0.5, linewidth = 0.3)+
    geom_point(aes(colour = as.factor(time)),
               size = 1.5, alpha = 0.75, stroke = 0)+
    scale_color_manual(name = "Time [dpi]", values = palette_d)+
    scale_x_continuous(name = "Cumulative distribution [%]", breaks = qnorm(b), labels = scales::percent(b, suffix = ""), limits = c(-4.2,4.2))+
    scale_y_continuous(name = "Relative fluorescence")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, fill = palette_d, stroke = 0)))+
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

##      Top 1% fluorescence values per group
top_percent <- 
    normalised_df %>%
    group_by(Treatment, time, Rep) %>%
    filter(ecdf_prob > 0.99) %>%  # Select top 1%
    summarise(
        n = n(),
        median_top = median(relative_rfu), 
        mean_top = mean(relative_rfu), 
        q99 = quantile(relative_rfu, 0.99), 
        q95 = quantile(relative_rfu, 0.95), 
        .groups = "drop") %>% 
    mutate(Time = factor(time))

#       Test for normality
top_percent %>% 
    group_by(Time) %>% 
    shapiro_test(q99)

#       Test for homoscedasticity
top_percent %>% 
    levene_test(q99 ~ Time)
fligner.test(q99 ~ Time, top_percent)

##      Are they different to zero?
top_percent %>% 
    group_by(Time) %>% 
    t_test(q99 ~ 0, mu = 1, var.equal = TRUE)
    
##      Are they different between sampling time
top_percent %>% 
    anova_test(q99 ~ Time)

p_value_top <- 
    top_percent %>% 
    tukey_hsd(q99 ~ Time) %>% 
    add_xy_position() %>% 
    mutate(xmin = as.numeric(group1),
           xmax = as.numeric(group2),
           y.position = y.position - 0.8) %>% 
    filter(p.adj < 0.05)

#       Plot 
plot_top_percent <-
    top_percent %>% 
    ggplot(aes(x = time, y = q99))+
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, alpha = 0.2)+
    geom_line(aes(group = interaction(time, Treatment)),
              position = position_dodge(width = 3),
              linewidth = 0.5)+
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.8)+
    add_pvalue(data = p_value_top, 
               xmin = "xmin", xmax = "xmax",
               label.size = 4, tip.length = 0.01,
               bracket.size = 0.3, step.increase = 0.001)+
    scale_x_continuous(name = "Time [d]", limits = c(-0.1, 30), breaks = c(1, 6, 28))+
    scale_y_continuous(name = bquote(atop("Fold change fluorescence", "(" * Q[99] * ", normalised by BHB)")),
                       limits = c(0, 7.5))+
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
plot_invitro_distribution
ggsave(here('results', 'figure_S3.pdf'), 
       dpi = 600, width = 180, height = 100, units = "mm")

plot_invitro_normalised_distribution + plot_top_percent +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(3, 1)) &
    theme(legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.spacing = unit(0, "pt"),
          legend.key.width = unit(0.5, "lines"),     # Width of the color box
          legend.key.height = unit(0.5, "lines"))
ggsave(here('results', 'figure_2.pdf'), 
       dpi = 600, width = 180, height = 60, units = "mm")
