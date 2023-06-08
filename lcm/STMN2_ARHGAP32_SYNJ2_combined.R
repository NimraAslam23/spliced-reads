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


# any cases with all three or two cryptic events?

intersect(STMN2_cryptic_cBio$case_submitter_id, ARHGAP32_cryptic_cBio$case_submitter_id) # none
intersect(STMN2_cryptic_cBio$case_submitter_id, SYNJ2_cryptic_cBio$case_submitter_id) # none
common_patients_ARHGAP32_SYNJ2 <- intersect(ARHGAP32_cryptic_cBio$case_submitter_id, SYNJ2_cryptic_cBio$case_submitter_id) 
      # 9 cases with both ARHGAP32 and SYNJ2 cryptic events

# which cancer type do you detect which cryptic event in? -----------------

    # df with only the 9 ARHGAP32+SYNJ2 cases 

ARHGAP32_cryptic_cBio |> 
       distinct() |> 
       filter(case_submitter_id %in% common_patients_ARHGAP32_SYNJ2) |> 
  left_join(SYNJ2_cryptic_cBio, by = "case_submitter_id") |> 
  distinct() |> 
  filter(sample_type.x == "Primary Tumor") |> 
  filter(sample_type.y == "Primary Tumor") |> 
  View()


common_patients_ARHGAP32_SYNJ2

# joining all three df together (STMN2 cryptic, ARHGAP32 cryptic, SYNJ2 cryptic) --------
    # to get the cases that have all three cryptic events

STMN2_cryptic_cBio <- STMN2_cryptic_cBio |> 
  select(c(STMN2_cryptic_coverage, STMN2_annotated_coverage, jir, case_submitter_id, cancer_type, gdc_primary_site,
           junction_coverage, junction_avg_coverage, rpm, Months.of.disease.specific.survival, Mutation.Count)) |> 
  rename("STMN2_jir" = "jir") |> 
  rename("STMN2_junction_coverage" = "junction_coverage") |> 
  rename("STMN2_junction_avg_coverage" = "junction_avg_coverage") |> 
  rename("STMN2_rpm" = "rpm") |> 
  rename("STMN2_months_survival" = "Months.of.disease.specific.survival") |> 
  rename("STMN2_mutation_count" = "Mutation.Count") |> 
  distinct()

ARHGAP32_cryptic_cBio <- ARHGAP32_cryptic_cBio |> 
  select(c(cryptic_count, anno_count, jir, case_submitter_id, junction_coverage, junction_avg_coverage, cancer_abbrev, 
           disease_survival_months, mutation_count)) |> 
  rename("ARHGAP32_cryptic_count" = "cryptic_count") |> 
  rename("ARHGAP32_anno_count" = "anno_count") |> 
  rename("ARHGAP32_jir" = "jir") |> 
  rename("ARHGAP32_junction_coverage" = "junction_coverage") |> 
  rename("ARHGAP32_junction_avg_coverage" = "junction_avg_coverage") |> 
  rename("ARHGAP32_months_survival" = "disease_survival_months") |> 
  rename("ARHGAP32_mutation_count" = "mutation_count") |> 
  distinct()

SYNJ2_cryptic_cBio <- SYNJ2_cryptic_cBio |> 
  select(c(cryptic_count, anno_count, jir, junction_coverage, junction_avg_coverage, case_submitter_id, 
           disease_survival_months, mutation_count)) |> 
  rename("SYNJ2_cryptic_count" = "cryptic_count") |> 
  rename("SYNJ2_anno_count" = "anno_count") |> 
  rename("SYNJ2_jir" = "jir") |> 
  rename("SYNJ2_junction_coverage" = "junction_coverage") |> 
  rename("SYNJ2_junction_avg_coverage" = "junction_avg_coverage") |> 
  rename("SYNJ2_months_survival" = "disease_survival_months") |> 
  rename("SYNJ2_mutation_count" = "mutation_count") |> 
  distinct()

STMN2_ARHGAP32_SYNJ2_cryptic <- STMN2_cryptic_cBio |> 
  inner_join(SYNJ2_cryptic_cBio, by="case_submitter_id") |> 
  inner_join(ARHGAP32_cryptic_cBio, by="case_submitter_id")

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
