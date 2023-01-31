setwd("/Users/nimraaslam/Documents/GitHub/spliced-reads/lcm")
getwd()
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)

# Data import -------------------------------------------------------------

spliced_reads_orig <- read.csv("spliced_read_from_lcm.csv")
View(spliced_reads_orig)
spliced_reads <- spliced_reads_orig

# Looking at data ---------------------------------------------------------

spliced_reads |> 
  distinct(sample_name) #21 unique samples

spliced_reads |> 
  distinct(junction_name) #770 splice junctions

group_column <- spliced_reads |> 
  select(junction_name, sample_name, n_spliced_reads) |> 
  complete(junction_name, sample_name,fill=list(n_spliced_reads = 0)) |> 
  mutate(disease = ifelse(grepl("ALS", sample_name),
                          "ALS",
                          "control"))

View(group_column)

group_column |>  filter(grepl("STMN2",junction_name)) |> 
  ggplot(aes(x = n_spliced_reads,
             fill = disease)) +
  stat_density() + 
  scale_x_continuous(trans = scales::pseudo_log_trans())   

STMN2 <- group_column |>  filter(grepl("STMN2",junction_name)) 

k <- wilcox.test(n_spliced_reads ~ disease, STMN2) |> 
  broom::tidy()

View(k)

# Create a table of the average count of each junction in each disease

group_column |> 
  group_by(junction_name,disease) |> 
  mutate(mean_n_spliced_reads = mean(n_spliced_reads)) |> 
  ungroup() |> 
  select(junction_name,disease,mean_n_spliced_reads) |> 
  unique() |> #have mean_n_spliced_reads for each junction_name separately for control and ALS
  pivot_wider(names_from = 'disease',
              values_from = 'mean_n_spliced_reads')

# Nest data by junction_name

group_column_nested <- group_column |> 
  group_by(junction_name) |> 
  nest()
  
View(group_column_nested)

# Run wilcoxon test on all the nested data

wilcox_tested_nested = group_column_nested |> 
  mutate(wc_res = map(data,~{broom::tidy(wilcox.test(.x$n_spliced_reads ~ .x$disease, exact=FALSE))})) |> 
  unnest(wc_res) |> 
  arrange(p.value)

View(wilcox_tested_nested)










