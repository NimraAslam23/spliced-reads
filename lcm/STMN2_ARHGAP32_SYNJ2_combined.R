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

# files to run for this script:
    # STMN2_cryptic_cBio
    # ARHGAP32_cryptic_cBio
    # SYNJ2_cryptic_cBio
    # all_mutation_data
    # mutation_clinical_data


# functions in this script ------------------------------------------------

read_save_mutation_data <- function(filename) {
  df <- janitor::clean_names(read.csv(filename, sep = "\t"))
  return(df)
}
# read_save_mutation_data("TCGA-VQ-A91N_mutations.tsv") |> View()

combine_mutation_data <- function(folder_path, pattern = "^TCGA.*_mutations.tsv$") {
  mutation_files <- list.files(folder_path,
                               pattern = pattern,
                               full.names = TRUE)
  #print(mutation_files)
  df <- purrr::map(mutation_files, read_save_mutation_data)
  patient_id <- purrr::simplify(purrr::map(mutation_files, basename))
  patient_id = gsub("_mutations.tsv", "", patient_id)
  non_empty_df <- df[sapply(df, function(x) nrow(x) > 0)]
  non_empty_patient_id <- patient_id[sapply(df, function(x) nrow(x) > 0)]
  df <- purrr::map2(non_empty_df, non_empty_patient_id, ~cbind(.x, case_submitter_id = .y))
  df <- data.table::rbindlist(df)
  return(df)
}

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
  left_join(SYNJ2_cryptic_cBio, by = "case_submitter_id", suffix = c("_ARHGAP32","_SYNJ2")) |> 
  distinct() |> 
  drop_na() |> 
  select(-cancer) |> 
  filter(sample_type == "Primary Tumor") |> 
  janitor::clean_names() 

  # within each patient, what is the fraction of cryptic to annotated reads?
common_ARHGAP32_SYNJ2 |> 
  select(case_submitter_id,gdc_primary_site,disease_survival_months,mutation_count,jir_arhgap32,jir_synj2) |> 
  pivot_longer(cols = c(jir_arhgap32,jir_synj2),
               names_to = "gene",
               values_to = "cryptic_jir") |>  
  mutate(gene = str_replace(gene, "jir_", "")) |> 
  mutate(cryptic_jir = 1 - cryptic_jir) |> 
  ggplot(aes(x = gene,y = cryptic_jir)) + 
  geom_boxplot()

  # within each patient, is ARHGAP32 or SYNJ2 cryptic expressed at a higher rate?
common_ARHGAP32_SYNJ2 |> 
  select(case_submitter_id,gdc_primary_site,disease_survival_months,mutation_count,jir_arhgap32,jir_synj2) |> 
  pivot_longer(cols = c(jir_arhgap32,jir_synj2),
               names_to = "gene",
               values_to = "cryptic_jir") |>  
  mutate(gene = str_replace(gene, "jir_", "")) |> 
  mutate(cryptic_jir = 1 - cryptic_jir) |> 
  ggplot(aes(x = gene,y = cryptic_jir,group = case_submitter_id)) + 
  geom_point() +
  geom_line()

common_ARHGAP32_SYNJ2 |> 
  select(case_submitter_id,gdc_primary_site,disease_survival_months,mutation_count,jir_arhgap32,jir_synj2) |> 
  mutate(jir_arhgap32 = 1 - jir_arhgap32 , jir_synj2 = 1 - jir_synj2) |> 
  mutate(arh_ratio = log2(jir_arhgap32 / jir_synj2)) |> 
  ggplot(aes(y = arh_ratio,x = 1)) + 
  geom_violin() + 
  geom_hline(yintercept = 0) + 
  ylab("Log2Fold Ratio of ARHGAP32 to SYNJ2 cryptic expression")

# df with STMN2 ARHGAP32 cases --------------------------------------------

common_STMN2_ARHGAP32 <- STMN2_cryptic_cBio |> 
  filter(sample_type == "Primary Tumor",
         case_submitter_id %in% common_patients_STMN2_ARHGAP32) |> 
  select(case_submitter_id, stmn2_cryptic_coverage, stmn2_annotated_coverage, jir) |> 
  left_join(ARHGAP32_cryptic_cBio, by = "case_submitter_id",suffix = c("_STMN2","_ARHGAP32")) |>   
  distinct() |> 
  filter(sample_type == "Primary Tumor") |> 
  janitor::clean_names() |> 
  select(-c(chromosome, start, end, strand, rail_id, junction_coverage, junction_avg_coverage)) 

# df with all common cases ------------------------------------------------

all_common_cases <- common_STMN2_ARHGAP32 |> 
  full_join(common_ARHGAP32_SYNJ2) 

all_common_cases |> 
  pivot_longer(cols = c(jir_stmn2, jir_arhgap32, jir_synj2),
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
  mutate(case_id = paste0(study_id, ":", case_submitter_id)) |> 
  distinct()

write_clip(unique(all_common_cases$case_id)) # copies all common patients to clipboard
write_clip(unique(all_common_cases$case_submitter_id))


# filter mutation_clinical_data for 14 common patients only ---------------

x <- mutation_clinical_data |> 
  filter(case_submitter_id %in% all_common_cases)

intersect(all_common_cases$case_submitter_id, mutation_clinical_data$case_submitter_id)


intersect(STMN2_cryptic_cBio$case_submitter_id, SYNJ2_cryptic_cBio$case_submitter_id)







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
