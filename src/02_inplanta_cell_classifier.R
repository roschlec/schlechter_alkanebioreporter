###     02 In planta - Single cell classifier

###     Libraries
library(tidyverse)
library(tidymodels)
library(here)
library(vip)
library(patchwork)
set.seed(14071990)

###     Load data
data <- 
    read.csv(here('data', 'train_data_in_planta.csv')) %>% 
    select(-X, -Label, -Mean, -ROI) %>% 
    filter(Cell.Status %in% c("Yes", "No"))

### Data Split
data_split <- 
    initial_split(data, prop = 0.75)

###   Create training and testing data
data_train <- training(data_split)
data_test <- testing(data_split)

### Cross-validation
folds <- vfold_cv(data_train, v = 10)

### Recipe
cell_rec <-
    recipe(Cell.Status ~ ., data = data_train) %>% 
    step_dummy(all_nominal_predictors()) %>% 
    step_zv(all_predictors()) %>% 
    step_normalize(all_predictors())

### Logistic Regression Model
### Define a model
lr_mod <-
    logistic_reg(penalty = tune(), mixture = 1) %>% 
    set_engine("glmnet", family = "binomial")

### Workflow
inplanta_lr_wkflw <-
    workflow() %>% 
    add_model(lr_mod) %>% 
    add_recipe(cell_rec)

### Tuning hyperparameters
lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))

lr_res <- 
    inplanta_lr_wkflw %>% 
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

best_lr <- 
    lr_res %>% 
    select_best(metric = "roc_auc")

inplanta_lr_wkflw_final <- 
    inplanta_lr_wkflw %>% 
    finalize_workflow(best_lr)

### Fit model
inplanta_lr_fit <- 
    inplanta_lr_wkflw_final %>%
    last_fit(data_split,
             metrics = metric_set(
                 recall, precision, f_meas,
                 accuracy, kap,
                 roc_auc, sens, spec)) 

lr_metrics <- 
    inplanta_lr_fit %>%
    collect_metrics()

lr_roc_curve <- 
    inplanta_lr_fit %>%
    collect_predictions() %>% 
    roc_curve(Cell.Status, .pred_No) %>% 
    autoplot()+
    theme_classic()+
    labs(title = "LR")+
    theme(aspect.ratio = 1)

lr_importance <- 
    inplanta_lr_fit %>% 
    extract_fit_parsnip() %>% 
    vip()+
    theme_classic()+
    labs(title = "LR")+
    theme(aspect.ratio = 1)

### Random Forest Model
rf_mod <-
    rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>% 
    set_engine("ranger", importance = "impurity") %>% 
    set_mode("classification")

### Workflow
inplanta_rf_wkflw <-
    workflow() %>% 
    add_model(rf_mod) %>% 
    add_recipe(cell_rec)

### Tuning hyperparameters
rf_res <- 
    inplanta_rf_wkflw %>% 
    tune_grid(folds,
              grid = 25,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))
rf_best <- 
    rf_res %>% 
    select_best(metric = "roc_auc")

## Final model
inplanta_rf_wkflw_final <- 
    inplanta_rf_wkflw %>% 
    finalize_workflow(rf_best)

### Fit model
inplanta_rf_fit <- 
    inplanta_rf_wkflw_final %>%
    last_fit(data_split,
             metrics = metric_set(
                 recall, precision, f_meas,
                 accuracy, kap,
                 roc_auc, sens, spec))

rf_metrics <- 
    inplanta_rf_fit %>%
    collect_metrics()

rf_roc_curve <- 
    inplanta_rf_fit %>%
    collect_predictions() %>% 
    roc_curve(Cell.Status, .pred_No) %>% 
    autoplot()+
    theme_classic()+
    labs(title = "RF")+
    theme(aspect.ratio = 1)

rf_importance <- 
    inplanta_rf_fit %>% 
    extract_fit_parsnip() %>% 
    vip() +
    theme_classic()+
    labs(title = "RF")+
    theme(aspect.ratio = 1)

###     Compare Models
metrics_both <-
    left_join(
        collect_metrics(inplanta_lr_fit),
        collect_metrics(inplanta_rf_fit),
        by = c(".metric", ".estimator"),
        suffix = c(".LR", ".RF")) %>% 
    select(.metric, .estimate.LR, .estimate.RF)


###     Plots
inplanta_plt_metrics <-
    metrics_both %>% 
    pivot_longer(cols = .estimate.LR:.estimate.RF) %>% 
    ggplot(aes(x = .metric, y = value, color = name))+
    geom_point(
        size = 3,
        position = position_dodge(width = 0.6))+
    geom_segment(
        aes(y = 0.5, yend = value),
        linewidth = 1,
        position = position_dodge(width = 0.6))+
    coord_flip()+
    theme_classic()+
    theme(aspect.ratio = 1)+
    scale_color_manual(name = "Model", 
                       values = c("#FF6B6B", "#009E9E"), 
                       label = c("Logistic Regression", "Random Forest"))+ 
    scale_y_continuous(
        name = "Mean Performance") +
    scale_x_discrete(
        name = "Metric",
        label = c("Accuracy", "F1-Score", "Kappa", "Precision",
                  "Recall", "ROC AUC", "Sensitivity", "Specificity"))

inplanta_plt_roc <-
    lr_roc_curve + rf_roc_curve

inplanta_plt_imp <-
    lr_importance + rf_importance

inplanta_plt_metrics + inplanta_plt_imp +
    plot_annotation(tag_levels = "A")+
    plot_layout(width = c(0.8, 2)) &
    theme(
        aspect.ratio = 1,
        plot.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.spacing = unit(0, "pt"),
        legend.key.width = unit(0.5, "lines"),     # Width of the color box
        legend.key.height = unit(0.5, "lines"))
ggsave(here('results', 'figure_S4.pdf'), 
       dpi = 600, width = 300, height = 100, units = "mm")

### Save models
saveRDS(inplanta_lr_fit, here('results', 'inplanta_lr.rds'))
saveRDS(inplanta_rf_fit, here('results', 'inplanta_rf.rds'))
