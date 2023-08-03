library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(rstatix)
library(knitr)

pdf(file = "figures_spliced_read_from_lcm.pdf")

## Introduction

# read in data
spliced_reads_orig <- read.csv("spliced_read_from_lcm.csv")
spliced_reads <- spliced_reads_orig


# unique samples and junctions
n_samp <- spliced_reads |> 
  distinct(sample_name) |> 
  nrow()

n_splice <- spliced_reads |> 
  distinct(junction_name) |> nrow() 

# disease group column
group_column <- spliced_reads |> 
  select(junction_name, sample_name, n_spliced_reads) |> 
  complete(junction_name, sample_name,fill=list(n_spliced_reads = 0)) |> 
  mutate(disease = ifelse(grepl("ALS", sample_name),
                          "ALS",
                          "control"))

# density plot p-value

als_group_column <- group_column |> 
  filter(disease=="ALS")

control_group_column <- group_column |> 
  filter(disease=="control")

als_control_group_column_wilcox_test <- wilcox.test(als_group_column$n_spliced_reads, control_group_column$n_spliced_reads) 
group_column |> mutate(p.value = als_control_group_column_wilcox_test$p.value)

# density plot: number of reads vs disease group

a_spliced_read_from_lcm <- group_column |> 
  ggplot(aes(x = n_spliced_reads, fill = disease)) +
  stat_density() + 
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  labs(
    title = "Number of Spliced Reads in ALS and Control LCM Neurons",
    x = "Number of Spliced Reads",
    y = "Density",
    fill = "Disease Group"
  ) +
  scale_fill_manual(labels=c('ALS', 'Control'),
                    values=c('ALS' = "tomato", 'control' = "dodgerblue")) +
  annotate("text", x = 5, y = 2.5,
           label = paste0("p-value = ", round(als_control_group_column_wilcox_test$p.value, 3)),
           hjust = 1, vjust = -1,
           color = "black", size = 4) +  #fontface = "bold") +
  theme(legend.position = "bottom")

print(a_spliced_read_from_lcm)

# mean spliced reads for each junction separately per disease group
mean_per_junction <- group_column |> 
  group_by(junction_name,disease) |> 
  mutate(mean_n_spliced_reads = mean(n_spliced_reads)) |> 
  ungroup() |> 
  select(junction_name,disease,mean_n_spliced_reads) |> 
  unique() |> 
  pivot_wider(names_from = 'disease',
              values_from = 'mean_n_spliced_reads')


## Wilcoxon test


# nesting
group_column_nested <- group_column |> 
  group_by(junction_name) |> 
  nest()

# wilcoxon test
wilcox_tested_nested = group_column_nested |>
  mutate(wc_res = map(data,~{broom::tidy(wilcox.test(
    .x$n_spliced_reads ~ .x$disease, exact=FALSE))})) |> 
  unnest(wc_res) |> 
  arrange(p.value)

significant_junctions <- wilcox_tested_nested |> 
  filter(p.value < 0.05)

### Significantly expressed junctions

# fraction of junctions higher in ALS or controls
mean_per_junction_sig <- mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control, "ALS", "control")) |> 
  janitor::tabyl(higher_value)

#kable(mean_per_junction_sig, caption = "A greater proportion of junctions are more highly expressed in controls.")


# density plot - significant p value
sig_reads_groups <- group_column |> 
  semi_join(significant_junctions, by=("junction_name")) 

# calculating p value
als_sig_reads_groups <- sig_reads_groups |> 
  filter(disease=="ALS")

control_sig_reads_groups <- sig_reads_groups |> 
  filter(disease=="control")

als_control_sig_reads_groups_wilcox_test <- wilcox.test(als_sig_reads_groups$n_spliced_reads, control_sig_reads_groups$n_spliced_reads) 

# density plot - significant spliced reads in ALS and controls

b_spliced_read_from_lcm <- sig_reads_groups |> 
  ggplot(aes(x = n_spliced_reads, 
             color = disease, 
             fill = disease)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  labs(
    title = stringr::str_wrap("Number of Spliced Reads in Significantly Differentially Expressed Junctions in ALS and Control LCM Neurons", width = 85),
    x = "Number of Spliced Reads",
    y = "Density",
    color = "Disease Group",
    fill = "Disease Group"
  ) +
  scale_fill_manual(labels=c('ALS', 'Control'),
                    values=c('ALS' = "tomato", 'control' = "dodgerblue")) +
  scale_color_manual(labels=c('ALS', 'Control'),
                     values=c('ALS' = "tomato", 'control' = "dodgerblue")) +
  annotate("text", x = 5, y = 1,
           label=paste0("p-value < 0.01"),
           hjust = 1, vjust= -1,
           color = "black", size = 4) +
  theme(legend.position = "bottom")

print(b_spliced_read_from_lcm)

## Cryptic expression of different junctions in ALS and control samples

c_spliced_read_from_lcm <- mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control,"ALS","control")) |> 
  separate(junction_name, sep = '\\|',into = c("gene","junc_cat","n_datasets_junction_found"),convert = TRUE) |> 
  ggplot(aes(x = junc_cat, fill = higher_value)) +
  geom_bar() +
  labs(
    title = "Expression of Junction Events in ALS and Control LCM Neurons",
    x = "Junction Category",
    y = "Number of Junction Events",
    fill = "Disease Group",
  ) +
  scale_fill_manual(labels=c('ALS', 'Control'),
                    values=c('ALS' = "tomato", 'control' = "dodgerblue")) +
  theme(legend.position = "bottom")

print(c_spliced_read_from_lcm)

## Datasets supporting cryptic and annotated STMN2 events

datasets_junction_found <- mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control,"ALS","control")) |> group_by(higher_value) |> 
  separate(junction_name, sep = '\\|',into = c("gene","junc_cat","n_datasets_junction_found"),convert = TRUE) |> 
  summarize(mean_n_datasets = mean(n_datasets_junction_found))

#kable(datasets_junction_found, caption = "Average number of datatsets the *STMN2* junction events are found in.")

## Investigating *STMN2* events

ptdp_orig <- read.csv("ptdp_mn_death.csv")

spliced_reads$sample_name <- gsub(".SJ.out", "", as.character(spliced_reads$sample_name))

ptdp_orig <- ptdp_orig |> 
  rename("sample_name" = "sample")

spliced_reads_pdtp <- spliced_reads |> 
  left_join(ptdp_orig, by = "sample_name")

spliced_reads_pdtp_select <- spliced_reads_pdtp |> 
  select(sample_name, n_spliced_reads, junction_name, pTDP.43, MN_death, pTDP_category)

### *STMN2* events are more highly expressed in low accumulation of pTDP-43

# STMN2 events vs pTDP43 and MN_death
d_spliced_read_from_lcm <- spliced_reads_pdtp_select |> 
  filter(grepl("STMN2", junction_name)) |> 
  drop_na() |> 
  ggplot(aes(x = pTDP_category, 
             y = n_spliced_reads
             #fill = pTDP_category
  )) +
  geom_boxplot() +
  theme(legend.position = "none") +
  labs(
    title = "Effect of pTDP-43 on STMN2 Splicing", 
    x = "Level of pTDP-43 Accumulation",
    y = "Number of STMN2 Spliced Reads") +
  stat_compare_means(hjust = 1, vjust = 14) 
#scale_fill_manual(labels=c('high', 'low'),
#values=c('tomato', 'dodgerblue'))

print(d_spliced_read_from_lcm)

e_spliced_read_from_lcm <- spliced_reads_pdtp_select |> 
  filter(grepl("STMN2", junction_name)) |> 
  drop_na() |> 
  ggplot(aes(x = pTDP_category,
             y = MN_death,
             fill = pTDP_category)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  labs(
    title = "Effect of pTDP-43 on Motor Neuron Viability", 
    x = "Level of pTDP-43 Accumulation",
    y = "Motor Neuron Death (%)") +
  stat_compare_means(hjust = 1.5, vjust = 7) +
  scale_fill_manual(labels=c('high', 'low'),
                    values=c('tomato', 'dodgerblue'))

print(e_spliced_read_from_lcm)

# reading in STMN2 data
stmn2_all_junctions_orig <- read_tsv("stmn2_all_junctions.bed",
                                     col_names = c("chrosome", "start", "end", "sample_name", 
                                                   "n_spliced_reads", "strand"))
stmn2_all_junctions <- unique(stmn2_all_junctions_orig)
stmn2_all_junctions$sample_name <- gsub(".SJ.out", "", as.character(stmn2_all_junctions$sample_name))

spliced_reads_pdtp_stmn2 <- stmn2_all_junctions |> 
  left_join(ptdp_orig, by=c("sample_name"))

spliced_reads_pdtp_select |> 
  filter(grepl("STMN2", junction_name))

stmn2_all_junctions |> 
  distinct(end)

spliced_reads_pdtp_stmn2 <- spliced_reads_pdtp_stmn2 |> 
  mutate(junction = ifelse(grepl("79616822", end),
                           "cryptic",
                           "annotated")) |> 
  mutate(disease = ifelse(grepl("ALS", sample_name),
                          "ALS",
                          "control")) |> 
  relocate(disease, .after = sample_name) |> 
  relocate(junction, .after = disease)

# STMN2 reads vs disease group and pTDP category
f_spliced_read_from_lcm <- spliced_reads_pdtp_stmn2 |> 
  ggplot(aes(x = junction, 
             y = n_spliced_reads, 
             fill = disease)) +
  geom_boxplot() +
  labs(
    title = "Expression of Annotated and Cryptic STMN2 Events in ALS and Control \nLCM Neurons",
    x = "Junction Event",
    y = "Number of Spliced Reads",
    fill = "Disease Group"
  ) +
  theme(legend.position = "bottom") + 
  scale_fill_manual(labels=c('ALS', 'Control'),
                    values=c('tomato', 'dodgerblue')) +
  annotate("text", x = 1.1, y = 75,
           label=paste0("Wilcoxon, p-value = 0.0014"),
           hjust = 1, vjust= -1,
           color = "black", size = 4) +
  theme(legend.position = "bottom",
        plot.title = element_text(size=13)) 

print(f_spliced_read_from_lcm)

g_spliced_read_from_lcm <- spliced_reads_pdtp_stmn2 |> 
  drop_na() |> 
  ggplot(aes(x = junction, 
             y = n_spliced_reads, 
             fill = pTDP_category)) +
  geom_boxplot() +
  labs(
    title = "Effect of pTDP-43 Accumulation on the Expression of Annotated and Cryptic \nSTMN2 Events",
    x = "Junction Event",
    y = "Number of Spliced Reads",
    fill = "Level of pTDP-43 Accumulation"
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(size=12)) + 
  stat_compare_means() +
  scale_fill_manual(labels=c('high', 'low'),
                    values=c('tomato', 'dodgerblue')) 

print(g_spliced_read_from_lcm)

# STMN2 reads vs MN death in annotated and cryptic events
h_spliced_read_from_lcm <- spliced_reads_pdtp_stmn2 |> 
  drop_na() |> 
  ggplot(aes(x = MN_death, 
             y = n_spliced_reads, 
             color = junction)) +
  geom_point()+
  labs(
    x = "Motor Neuron Death (%)",
    y = "Number of Spliced Reads",
    color = "Junction",
    title = "Relationship between STMN2 Splice Events and Motor Neuron Viability"
  ) +
  stat_cor(show.legend = FALSE) +
  scale_color_manual(labels=c('annotated', 'cryptic'),
                     values=c('dodgerblue', 'tomato')) +
  theme(legend.position = "bottom")

print(h_spliced_read_from_lcm)

# psi table
psi_stmn2 <- spliced_reads_pdtp_stmn2 |> 
  select(-start,-end) |> 
  unique() |>
  pivot_wider(names_from = 'junction',
              values_from = 'n_spliced_reads',
              values_fill = 0) 

# calculating psi values for STMN2 junctions
psi_stmn2 <- psi_stmn2 |> 
  group_by(sample_name) |> 
  mutate(total_counts = annotated + cryptic) |> 
  mutate(psi = cryptic / total_counts)

# STMN2 psi value vs pTDP and MN death
i_spliced_read_from_lcm <- psi_stmn2 |> 
  drop_na() |> 
  ggplot(aes(x = pTDP_category, 
             y = psi)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  labs(
    title = "Effect of pTDP-43 on STMN2 Cryptic Splicing",
    x = "Level of pTDP-43 Accumulation",
    y = "psi"
  ) + 
  stat_compare_means()

print(i_spliced_read_from_lcm)

j_spliced_read_from_lcm <- psi_stmn2 |> 
  drop_na() |> 
  ggplot(aes(x = pTDP_category, 
             y = psi,
             fill = pTDP_category)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  labs(
    title = "Effect of pTDP-43 on STMN2 Cryptic Splicing",
    x = "Level of pTDP-43 Accumulation",
    y = "psi"
  ) + 
  stat_compare_means() +
  scale_fill_manual(labels=c("high","low"),
                    values=c("tomato","dodgerblue"))

print(j_spliced_read_from_lcm)

k_spliced_read_from_lcm <- psi_stmn2 |> 
  drop_na() |> 
  ggplot(aes(x = MN_death, 
             y = psi, 
             color = pTDP_category)) +
  geom_point() +
  labs(
    title = "Relationship between STMN2 Cryptic Splicing and Motor Neuron Viability",
    x = "Motor neuron death (%)",
    y = "psi",
    color = "Level of pTDP-43 Accumulation"
  ) + 
  stat_cor(show.legend = FALSE) +
  scale_color_manual(labels=c('high','low'),
                     values=c('tomato','dodgerblue')) +
  theme(legend.position = "bottom")

print(k_spliced_read_from_lcm)

dev.off()
