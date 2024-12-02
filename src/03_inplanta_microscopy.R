library(here)
library(tidyverse)
library(tidymodels)
set.seed(14071990)

## Load model
cell_lr_fit <- readRDS(here('results', 'lr_model.rds'))

## Load data
data <- read.csv(here('data', 'inplanta_single_cell.csv'))
data_select <- 
    data %>% 
    select(-X, Circ = Circ., Perim = Perim.)

background <- read.csv(here('data', 'inplanta_background.csv'))
bg_values <- background %>% 
    na.omit() %>% 
    select(Label, Mean) %>% 
    mutate(Label = str_replace(Label, ".czi:.*", "")) %>% 
    group_by(Label) %>% 
    summarise(value = median(Mean)) %>% 
    separate(Label, into = c("Time", "Rep", "Img"), sep = "-")

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
    mutate(rfu = Mean - value) %>% 
    filter(rfu > 0)

sc_fl_df %>% 
    ggplot(aes(sample = rfu, group = Rep))+
    facet_grid(cols = vars(Time))+
    geom_qq()

# Normal probability
sc_fl_label <- sc_fl_df %>% 
    group_by(Time, Rep) %>% 
    tally() %>% 
    mutate(label = paste(Time, Rep, sep = "_"))

prob_plant <- sc_fl_df %>% 
    group_by(Time, Rep) %>% 
    group_map(~ df_ecdf(.x$rfu)) %>% 
    setNames(sc_fl_label$label) %>% 
    bind_rows(.id = "label") %>% 
    separate(label, into = c('Time', 'Rep'), sep = "_", remove = FALSE)
prob_plant$time_d <- factor(prob_plant$time_d, levels = c("0", "2", "7"))
prob_plant$trt <- factor(prob_plant$trt, levels = c("inoculum", "PFF2"))

prob_plant %>% 
    filter(ecdf < 1) %>% 
    ggplot(aes(x = qnorm(ecdf), y = x, colour = Time, group = Rep))+
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