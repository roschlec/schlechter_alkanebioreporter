library(here)
library(tidyverse)
library(patchwork)

#####   CFU     #####
##      Data
cfu_files <- 
    list.files(path = here('data'), 
               pattern = "\\w{3}\\d{1}_inplanta_cfu", 
               full.names = TRUE)

cfu_data <- 
    lapply(seq_along(cfu_files), function (i) {
        df <- read.csv(cfu_files[i]) %>% 
            group_by(replicate, time_d) %>% 
            summarise(
                CFU_gFW = mean(CFU_gFW),
                logCFU = mean(logCFU), .groups = "drop") %>% 
            arrange(time_d)
        df$ID <- c("rep1", "rep2")[i]  # Add ID column
        return(df)
    })

cfu_df <-
    bind_rows(cfu_data) %>% 
    na.exclude()

cfu_df %>% 
    group_by(ID, time_d) %>% 
    summarise(median = median(CFU_gFW),
              q75 = quantile(CFU_gFW, probs = .75),
              iqr = IQR(CFU_gFW),
              I = q75 + 1.5 * iqr,
              .groups = "drop")

##      Stats
cfu_df %>% 
    group_by(ID) %>% 
    wilcox_test(logCFU ~ time_d, p.adjust.method = "holm")

##      Plots
plt_cfu <- 
    cfu_df %>% 
    ggplot(aes(x = time_d, y = logCFU, group = ID))+
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

plt_cfu
ggsave(here("results", "inplanta_cfu.png"), dpi = 300, width = 4, height = 4)
