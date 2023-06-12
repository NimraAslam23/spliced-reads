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

common_patients_STMN2_ARHGAP32 <- intersect(STMN2_cryptic_cBio$case_submitter_id, ARHGAP32_cryptic_cBio$case_submitter_id) 
      # 5 cases
intersect(STMN2_cryptic_cBio$case_submitter_id, SYNJ2_cryptic_cBio$case_submitter_id) # none
common_patients_ARHGAP32_SYNJ2 <- intersect(ARHGAP32_cryptic_cBio$case_submitter_id, SYNJ2_cryptic_cBio$case_submitter_id) 
      # 31 cases with both ARHGAP32 and SYNJ2 cryptic events

# which cancer type do you detect which cryptic event in? -----------------

# df with ARHGAP32 SYNJ2 cases ------------------------------------------

common_ARHGAP32_SYNJ2 <- ARHGAP32_cryptic_cBio |> 
  filter(case_submitter_id %in% common_patients_ARHGAP32_SYNJ2) |> 
  select(case_submitter_id, start, end, strand, anno_count, cryptic_count, jir) |> 
  left_join(SYNJ2_cryptic_cBio, by = "case_submitter_id") |> 
  distinct() |> 
  select(-cancer) |> 
  filter(sample_type == "Primary Tumor") |> 
  janitor::clean_names() |> 
  rename_all(~str_replace_all(.x, "_x", "_ARHGAP32")) |> 
  rename_all(~str_replace_all(.x, "_y", "_SYNJ2")) 

  # within each patient, what is the fraction of cryptic to annotated reads?
common_ARHGAP32_SYNJ2 |> 
  select(case_submitter_id,gdc_primary_site,disease_survival_months,mutation_count,jir_ARHGAP32,jir_SYNJ2) |> 
  pivot_longer(cols = c(jir_ARHGAP32,jir_SYNJ2),
               names_to = "gene",
               values_to = "cryptic_jir") |>  
  mutate(gene = str_replace(gene, "jir_", "")) |> 
  mutate(cryptic_jir = 1 - cryptic_jir) |> 
  ggplot(aes(x = gene,y = cryptic_jir)) + 
  geom_boxplot()

  # within each patient, is ARHGAP32 or SYNJ2 cryptic expressed at a higher rate?
common_ARHGAP32_SYNJ2 |> 
  select(case_submitter_id,gdc_primary_site,disease_survival_months,mutation_count,jir_ARHGAP32,jir_SYNJ2) |> 
  pivot_longer(cols = c(jir_ARHGAP32,jir_SYNJ2),
               names_to = "gene",
               values_to = "cryptic_jir") |>  
  mutate(gene = str_replace(gene, "jir_", "")) |> 
  mutate(cryptic_jir = 1 - cryptic_jir) |> 
  ggplot(aes(x = gene,y = cryptic_jir,group = case_submitter_id)) + 
  geom_point() +
  geom_line()

common_ARHGAP32_SYNJ2 |> 
  select(case_submitter_id,gdc_primary_site,disease_survival_months,mutation_count,jir_ARHGAP32,jir_SYNJ2) |> 
  mutate(jir_ARHGAP32 = 1 - jir_ARHGAP32 , jir_SYNJ2 = 1 - jir_SYNJ2) |> 
  mutate(arh_ratio = log2(jir_ARHGAP32 / jir_SYNJ2)) |> 
  ggplot(aes(y = arh_ratio,x = 1)) + 
  geom_violin() + 
  geom_hline(yintercept = 0) + 
  ylab("Log2Fold Ratio of ARHGAP32 to SYNJ2 cryptic expression")

# df with STMN2 ARHGAP32 cases --------------------------------------------

common_STMN2_ARHGAP32 <- STMN2_cryptic_cBio |> 
  filter(sample_type == "Primary Tumor",
         case_submitter_id %in% common_patients_STMN2_ARHGAP32) |> 
  select(case_submitter_id, stmn2_cryptic_coverage, stmn2_annotated_coverage, jir) |> 
  left_join(ARHGAP32_cryptic_cBio, by = "case_submitter_id") |> 
  distinct() |> 
  filter(sample_type == "Primary Tumor") |> 
  janitor::clean_names() |> 
  select(-c(chromosome, start, end, strand, rail_id, junction_coverage, junction_avg_coverage)) |> 
  rename_all(~str_replace_all(.x, "_x", "_STMN2")) |> 
  rename_all(~str_replace_all(.x, "_y", "_ARHGAP32")) |> 
  rename("cryptic_count_STMN2" = "stmn2_cryptic_coverage",
         "anno_count_STMN2" = "stmn2_annotated_coverage",
         "anno_count_ARHGAP32" = "anno_count",
         "cryptic_count_ARHGAP32" = "cryptic_count")

# df with all common cases ------------------------------------------------

all_common_cases <- common_STMN2_ARHGAP32 |> 
  full_join(common_ARHGAP32_SYNJ2) 

all_common_cases |> 
  pivot_longer(cols = c(jir_STMN2, jir_ARHGAP32, jir_SYNJ2),
               names_to = "gene",
               values_to = "cryptic_jir") |> 
  mutate(gene = str_replace(gene, "jir_", "")) |> 
  mutate(cryptic_jir = 1 - cryptic_jir) |> 
  ggplot(aes(x = gene, 
             y = cryptic_jir,
             group = case_submitter_id)) +
  geom_point() +
  geom_line()


# case set of all 14 common patients --------------------------------------

all_common_cases <- all_common_cases |> 
  mutate(case_id = paste0(study_id, ":", case_submitter_id))

write_clip(unique(all_common_cases$case_id)) # copies all common patients to clipboard
write_clip(unique(all_common_cases$case_submitter_id))




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
