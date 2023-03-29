library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(repurrrsive)
library(jsonlite)
library(clipr)
library(naniar)
library(rstatix)
library(knitr)
library(snapcount)
library(ggsignif)

# raw cryptic counts and psi column ---------------------------------------

STMN2_cryptic_rawcounts_psi <- STMN2_query |> 
  select(gdc_cases.submitter_id, anno_count, cryptic_count, jir) |> 
  mutate(total_count = anno_count + cryptic_count) |> 
  mutate(STMN2_psi = cryptic_count / total_count) |> 
  select(-total_count) |> 
  rename("STMN2_anno_count" = "anno_count",
         "STMN2_cryptic_count" = "cryptic_count",
         "STMN2_jir" = "jir")

ARHGAP32_cryptic_rawcounts_psi <- ARHGAP32_query |> 
  select(gdc_cases.submitter_id, anno_count, cryptic_count, jir) |> 
  mutate(total_count = anno_count + cryptic_count) |> 
  mutate(ARHGAP32_psi = cryptic_count / total_count) |> 
  select(-total_count) |> 
  rename("ARHGAP32_anno_count" = "anno_count",
         "ARHGAP32_cryptic_count" = "cryptic_count",
         "ARHGAP32_jir" = "jir")

SYNJ2_cryptic_rawcounts_psi <- SYNJ2_query |> 
  select(gdc_cases.submitter_id, anno_count, cryptic_count, jir) |> 
  mutate(total_count = anno_count + cryptic_count) |> 
  mutate(SYNJ2_psi = cryptic_count / total_count) |> 
  select(-total_count) |> 
  rename("SYNJ2_anno_count" = "anno_count",
         "SYNJ2_cryptic_count" = "cryptic_count",
         "SYNJ2_jir" = "jir")

# combined table with STMN2 ARHGAP32 SYNJ2 raw jir psi --------------------

STMN2_ARHGAP32_SYNJ2_combined <- STMN2_cryptic_rawcounts_psi |> 
  full_join(ARHGAP32_cryptic_rawcounts_psi, by = "gdc_cases.submitter_id") |> 
  full_join(SYNJ2_cryptic_rawcounts_psi, by = "gdc_cases.submitter_id") |> 
  relocate(ARHGAP32_psi, .before = SYNJ2_psi) |> 
  relocate(STMN2_psi, .before = ARHGAP32_psi) |> 
  relocate(ARHGAP32_jir, .before = SYNJ2_jir) |> 
  relocate(STMN2_jir, .before = ARHGAP32_jir) |> 
  unique()

# do patients with high cryptic STMN2 have high cryptic expression --------

STMN2_ARHGAP32_SYNJ2_combined |> 
  ggplot(aes(x = STMN2_psi, y = ARHGAP32_psi)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE)

STMN2_ARHGAP32_SYNJ2_combined |> 
  ggplot(aes(x = STMN2_psi, y = SYNJ2_psi)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE)

STMN2_ARHGAP32_SYNJ2_combined |> 
  ggplot(aes(x = SYNJ2_psi, y = ARHGAP32_psi)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE)
