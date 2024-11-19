library(tidyverse)
library(tidymodels)
library(here)
library(vip)
set.seed(14071990)

### Load data
data <- read.csv(here('data', 'train_data', 'train_data.csv')) %>% 
    select(-Label, -ROI, -Mean)
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
lr_mod <-
    logistic_reg(penalty = tune(), mixture = 1) %>% 
    set_engine("glmnet", family = "binomial")

### Workflow
cell_lr_wkflw <-
    workflow() %>% 
    add_model(lr_mod) %>% 
    add_recipe(cell_rec)
cell_lr_wkflw

### Cross-validation
folds <- vfold_cv(data_train, v = 10)
folds

### Tuning hyperparameters
lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))

lr_res <- 
    cell_lr_wkflw %>% 
    tune_grid(folds,
              grid = lr_reg_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))

#   Evaluate best penalty
lr_plot <- 
    lr_res %>% 
    collect_metrics() %>% 
    ggplot(aes(x = penalty, y = mean)) + 
    geom_point() + 
    geom_line() + 
    ylab("Area under the ROC Curve") +
    scale_x_log10(labels = scales::label_number())
lr_plot

lr_res %>% 
    show_best(metric = "roc_auc", n = 1)

best_lr <- 
    lr_res %>% 
    select_best(metric = "roc_auc")

cell_lr_wkflw_final <- 
    cell_lr_wkflw %>% 
    finalize_workflow(best_lr)

### Fit model
cell_lr_fit <- 
    cell_lr_wkflw_final %>%
    last_fit(data_split) 

lr_metrics <- cell_lr_fit %>%
    collect_metrics()

lr_roc_curve <- cell_lr_fit %>%
    collect_predictions() %>% 
    roc_curve(Cell, .pred_No) %>% 
    autoplot()

lr_importance <- cell_lr_fit %>% 
    extract_fit_parsnip() %>% 
    vip()+
    theme_bw()

##  Save files
#   MODEL
saveRDS(cell_lr_fit, here('results', 'lr_model.rds'))
#   METRICS
write_csv(lr_metrics, here('results', 'lr_metrics.csv'))
#   ROC CURVE PLOT
ggsave(plot = lr_roc_curve, filename = here('results', 'lr_roc.pdf'), 
       dpi = 300, width = 5, height = 5)
#   IMPORTANCE
ggsave(plot = lr_importance, filename = here('results', 'lr_importance.pdf'), 
       dpi = 300, width = 7, height = 5)
