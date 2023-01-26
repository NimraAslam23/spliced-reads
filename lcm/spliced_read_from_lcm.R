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

    #? why does code below not work
spliced_reads |> 
  filter(sample_name, contains("junction")) |> 
  head()

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






   
   

complete(group_column, junction_name, sample_name)



# "Group" column ----------------------------------------------------------

spliced_reads |> 
  mutate(group =
    group_by(sample_name) |> 
           ALS = select(sample_name, contains("ALS") & controls = select(sample_name, contains("Controls")))
    .before = sample_name
  )

spliced_reads |> 
  mutate(group =
           case_when(sample_name, contains("ALS") ~ 'ALS',
                     sample_name, contains("Control") ~ 'control')
  )

                        
# Are any junctions significantly higher in ALS vs controls? --------------

    #? new column - number of spliced reads for each junction?

test <- t.test(formula = n_spliced_reads ~ group,
               data = spliced_reads)




