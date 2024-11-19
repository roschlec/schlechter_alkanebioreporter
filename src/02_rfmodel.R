library(tidyverse)
library(tidymodels)
library(here)
set.seed(14071990)

### Load data
data <- read.csv(here('data', 'train_data', 'train_data.csv')) %>% 
    select(-Label, -ROI)
head(data)

### Data Split
data_split <- initial_split(data, prop = 0.75)

#   Create training and testing data
data_train <- training(data_split)
data_test <- testing(data_split)

### Recipe
cell_rec <-
    recipe(Cell ~ ., data = data_train)

### Define a model
rf_mod <-
    rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>% 
    set_engine("ranger", importance = "impurity") %>% 
    set_mode("classification")

### Workflow
cell_rf_wkflw <-
    workflow() %>% 
    add_model(rf_mod) %>% 
    add_recipe(cell_rec)
cell_rf_wkflw

### Cross-validation
folds <- vfold_cv(data_train, v = 10)
folds

### Tuning hyperparameters
rf_res <- 
    cell_rf_wkflw %>% 
    tune_grid(folds,
              grid = 25,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))

rf_res %>% 
    show_best(metric = "roc_auc")

rf_best <- 
    rf_res %>% 
    select_best(metric = "roc_auc")

## Final model
cell_rf_wkflw_final <- 
    cell_rf_wkflw %>% 
    finalize_workflow(rf_best)

### Fit model
cell_rf_fit <- 
    cell_rf_wkflw_final %>%
    last_fit(data_split) 

rf_metrics <- cell_rf_fit %>%
    collect_metrics()

rf_roc_curve <- cell_rf_fit %>%
    collect_predictions() %>% 
    roc_curve(Cell, .pred_No) %>% 
    autoplot()+
    theme_bw()

rf_importance <- cell_rf_fit %>% 
    extract_fit_parsnip() %>% 
    vip() +
    theme_bw()

##  Save files
#   MODEL
saveRDS(cell_rf_fit, here('results', 'rf_model.rds'))
#   METRICS
write_csv(rf_metrics, here('results', 'rf_metrics.csv'))
#   ROC CURVE PLOT
ggsave(plot = rf_roc_curve, filename = here('results', 'rf_roc.pdf'), 
       dpi = 300, width = 5, height = 5)
#   IMPORTANCE
ggsave(plot = rf_importance, filename = here('results', 'rf_importance.pdf'), 
       dpi = 300, width = 7, height = 5)



cell_rf_fit %>% 
    extract_workflow() %>% 
    augment(data_test) %>% 
    select(Cell, .pred_class, .pred_No, .pred_Yes)
