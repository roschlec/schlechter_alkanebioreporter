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

background <- read.csv(here('data', 'background.csv'))

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
    mutate(Label = str_replace(Label, "-\\d{2}\\.czi:.*", "")) %>% 
    separate(Label, into = c("Channel", "Time", "Rep"), sep = "-") %>% 
    select(Time, Rep, Area, Mean)


