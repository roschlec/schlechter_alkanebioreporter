##  Clean data
library(tidyverse)
library(here)
library(scales)
library(RColorBrewer)
library(ggpubr)

##  Colour palette
palette_invitro <- c("diesel" = "#1F6683", "no_diesel" = "#DD3C51", "sterile" = "#D1C7B5")
palette_d <- c("1" = "#A6BDDB", "6" = "#02818A", "28" = "#014636")
palette_d_plant <- c("0-inoculum" = "#FED976", "0-PFF2" = "#99D8C9", "2-PFF2" = "#41AE76", "7-PFF2" = "#00441B")

##  Labels
trt_invitro <- c("diesel" = "BHB + 1% v/v Diesel",
                 "no_diesel" = "BHB",
                 "sterile" = "BHB (sterile)")

## Create breaks for x-axis (probability plots)
b <- c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999)

##  Functions
#   Probability distribution function
df_ecdf <- function(x){
    temp_ecdf <- ecdf(x)
    sort_prob <- temp_ecdf(sort(x))
    df <- data.frame(ecdf = sort_prob, x = sort(x))
    return(df)
}

### Function
mode <- function(codes){
    which.max(tabulate(codes))
}
