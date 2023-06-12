# gene_name = "SYNJ2" 
# snapcount_coords_cryptic = "chr6:158017291-158019983"
# snapcount_coords_annotated = "chr6:158017291-158028755"
# strand_code = "+"

library(tidyverse)
library(repurrrsive)
library(jsonlite)
library(clipr)
library(naniar)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
library(snapcount)
library(ggsignif)

# query - cryptic and anno ------------------------------------------------

SYNJ2_cryptic_query <- count_query(gene_name = "SYNJ2", snapcount_coords = "chr6:158017291-158019983", strand_code = "+")
SYNJ2_anno_query <- count_query(gene_name = "SYNJ2", snapcount_coords = "chr6:158017291-158028755", strand_code = "+")

# join cryptic and anno df, add jir column --------------------------------

SYNJ2_cryptic_query <- SYNJ2_cryptic_query |> rename("cryptic_count" = "count")
SYNJ2_anno_query <- SYNJ2_anno_query |> rename("anno_count" = "count")

SYNJ2_query <- SYNJ2_cryptic_query |> 
  select(gdc_cases.submitter_id, cryptic_count) |> 
  left_join(SYNJ2_anno_query, by = "gdc_cases.submitter_id")

SYNJ2_query <- SYNJ2_query |> 
  relocate(cryptic_count, .after = anno_count)

SYNJ2_query <- SYNJ2_query |> 
  mutate(jir = anno_count/(anno_count + cryptic_count)) |> 
  relocate(jir, .after = cryptic_count)

write_clip(SYNJ2_query$gdc_cases.submitter_id)

# clinical data from TCGA (case set) --------------------------------------

SYNJ2_clinical_orig <- read.csv("SYNJ2_clinical.tsv", sep = "\t", header = TRUE,
                                   na.strings = "", fill = TRUE)

SYNJ2_clinical <- SYNJ2_clinical_orig

SYNJ2_clinical <- SYNJ2_clinical|> 
  select(case_submitter_id, project_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage) |> 
  separate(project_id, into = c("project", "cancer_type")) |> 
  select(-project)

# bar plot - cancer types with SYNJ2 events ----------------------------

SYNJ2_clinical |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_type)))) +
  geom_bar() +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = ""
  ) +
  theme(plot.title = element_text(size=10))

# join query table with clinical table ------------------------------------

SYNJ2_clinical_jir <- SYNJ2_query |> rename("case_submitter_id" = "gdc_cases.submitter_id") |> 
  left_join(SYNJ2_clinical, by=c("case_submitter_id")) |> 
  select(-gdc_cases.diagnoses.tumor_stage, -age_at_index, -gdc_cases.demographic.gender, 
         -ajcc_pathologic_stage, -ethnicity, -race, -cancer_type, -gender) |> 
  rename("cancer_type" = "gdc_cases.project.name") |> 
  rename("gdc_primary_site" = "gdc_cases.project.primary_site") |> 
  rename("sample_type" = "gdc_cases.samples.sample_type") |> 
  rename("cgc_primary_site" = "cgc_case_primary_site")

SYNJ2_clinical_jir_cryptic <- SYNJ2_clinical_jir |> 
  filter(cryptic_count > 2) 

# primary sites of cancers with cryptic SYNJ2 events -------------------

SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(gdc_primary_site)))) +
  geom_bar() +
  coord_flip() +
  labs(
    x = "Primary Site of Cancer",
    y = "Number of Cases",
    title = ""
  ) +
  theme(plot.title = element_text(size=10))

# cryptic SYNJ2 expression in different cancer sites -------------------

SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = cryptic_count,
             y = fct_reorder(gdc_primary_site, cryptic_count, median))) +
  geom_boxplot() +
  labs(
    x = "SYNJ2 Cryptic Coverage",
    y = "Primary Site of Cancer"
  ) +
  theme(plot.title = element_text(size=10))

# which cancers have the most cryptic SYNJ2 events? --------------------

SYNJ2_events_different_cancers <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(cancer_type) |> arrange(-percent)

# where are the cancers with cryptic SYNJ2 events located? -------------

SYNJ2_events_primary_sites <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(gdc_primary_site) |> arrange(-percent)

# cBioPortal clinical data, join with SYNJ2 table ----------------------

SYNJ2_cryptic_cBio <- SYNJ2_clinical_jir_cryptic |> 
  left_join(cBio_clinical, by = "case_submitter_id") |> 
  select(-cgc_primary_site, -Cancer.Type, -Overall.Survival..Months.) |>
  rename("disease_survival_months" = "Months.of.disease.specific.survival") |> 
  relocate(case_submitter_id, .after = Sample.ID) |> 
  relocate(cancer_type, .after = case_submitter_id) |> 
  relocate(gdc_primary_site, .after = cancer_abbrev) |> 
  relocate(sample_type, .after = gdc_primary_site) |> 
  janitor::clean_names()

# boxplot - n mutations for each cancer type ------------------------------

SYNJ2_cryptic_cBio <- SYNJ2_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric)

SYNJ2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 500) |> 
  ggplot(aes(x = fct_reorder(cancer_abbrev, mutation_count, median),
             y = mutation_count)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Mutation Count"
  ) +
  geom_signif(comparisons = list(c("ESCA", "HNSC")),
              map_signif_level = TRUE,
              y_position = c(450))

# cryptic coverage vs mutation count --------------------------------------

SYNJ2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 300) |> 
  ggplot(aes(x = as.factor(cryptic_count),
             y = mutation_count)) +
  geom_boxplot() +
  labs(
    x = "Number of SYNJ2 cryptic events",
    y = "Mutation Count"
  ) + 
  geom_signif(comparisons = list(c("8", "9"), c("7","9")),
              map_signif_level = TRUE,
              y_position = c(230,210))

SYNJ2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 500) |> 
  ggplot(aes(x = cryptic_count, 
             y = mutation_count)) +
  geom_hex() +
  labs(
    x = "Number of SYNJ2 cryptic events",
    y = "Mutation Count"
  )

# Mutational burden

total_mutations_each_cancer <- cBio_clinical |> 
  drop_na(Mutation.Count) |> 
  filter(grepl("^\\d+$", Mutation.Count)) |>
  group_by(cancer_abbrev) |> 
  summarise(total_mutations = sum(as.numeric(Mutation.Count)))

total_mutations_each_cancer_with_SYNJ2cryptic <- SYNJ2_cryptic_cBio |> 
  drop_na(mutation_count) |> 
  filter(mutation_count < 2500) |> 
  group_by(cancer_abbrev) |> 
  summarise(total_mutations_cryptic = sum(mutation_count))

mutations_each_cancer_general_vs_SYNJ2cryptic <- total_mutations_each_cancer |> 
  left_join(total_mutations_each_cancer_with_SYNJ2cryptic, by=c("cancer_abbrev")) |> 
  mutate(percent_with_cryptic = (total_mutations_cryptic/total_mutations)*100)
#16% of all mutations in LGG cancer are in cases with STMN2 cryptic events

# Fraction of each cancer that has cryptic STMN2 events 

total_each_cancer <- cBio_clinical |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev = n())

total_each_cancer_with_SYNJ2cryptic <- SYNJ2_cryptic_cBio |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev_with_cryptic = n())

total_each_cancer_general_vs_SYNJ2cryptic <- total_each_cancer |> 
  left_join(total_each_cancer_with_SYNJ2cryptic, by=c("cancer_abbrev")) |> 
  mutate(percent_with_cryptic = (total_cancer_abbrev_with_cryptic/total_cancer_abbrev)*100)

total_each_cancer_general_vs_SYNJ2cryptic |> 
  drop_na() |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent_with_cryptic), y = percent_with_cryptic)) +
  geom_bar(aes(fill = cancer_abbrev), stat = "identity") +
  labs(
    x = "Cancer Type",
    y = "Percentage of cases with cryptic SYNJ2 events"
  ) +
  theme(legend.position = "none")

# calculating mutation burden

cBio_clinical |>
  left_join(SYNJ2_clinical_jir,by = c("case_submitter_id")) |> 
  rename("synj2_cryptic_coverage" = "cryptic_count") |> 
  rename("synj2_annotated_coverage" = "anno_count") |> 
  janitor::clean_names() |> 
  select(case_submitter_id, study_id, cancer_abbrev, mutation_count, synj2_cryptic_coverage, synj2_annotated_coverage) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  unique() |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(synj2_cryptic_detected = synj2_cryptic_coverage >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_synj2 = sum(synj2_cryptic_detected,na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(synj2_cryptic_detected)) |> 
  filter(n_detected_synj2 > 2) |> 
  ggplot(aes(x = synj2_cryptic_detected,
             y = mutation_count)) + 
  geom_boxplot(aes(fill = synj2_cryptic_detected)) + 
  labs(x = "Cancer Type",
       y = "Mutation Count") +
  facet_wrap(~study_start) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.format")


# Survival comparisons ----------------------------------------------------

survival_SYNJ2_cryptic <- cBio_clinical |>
  left_join(SYNJ2_clinical_jir,by = c("case_submitter_id")) |> 
  rename("synj2_cryptic_coverage" = "cryptic_count") |> 
  rename("synj2_annotated_coverage" = "anno_count") |> 
  janitor::clean_names() |> 
  select(case_submitter_id, study_id, sample_id, cancer_abbrev, mutation_count, synj2_cryptic_coverage, synj2_annotated_coverage, months_of_disease_specific_survival, overall_survival_months, disease_specific_survival_status) |> 
  mutate(months_of_disease_specific_survival = as.numeric(months_of_disease_specific_survival)) |> 
  filter(!is.na(months_of_disease_specific_survival) & months_of_disease_specific_survival != "NA") |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(synj2_cryptic_detected = synj2_cryptic_coverage >= 2) |> 
  group_by(cancer_abbrev) |> 
  mutate(n_detected_synj2 = sum(synj2_cryptic_detected,na.rm = TRUE)) |>
  ungroup() |> 
  filter(!is.na(synj2_cryptic_detected)) |> 
  filter(n_detected_synj2 > 2) |> 
  mutate(synj2_cryptic_detected = as.logical(synj2_cryptic_detected)) |> #changes FALSE and TRUE to 0 and 1 respectively
  distinct()

# survival status column
survival_SYNJ2_cryptic <- survival_SYNJ2_cryptic |> 
  filter(!is.na(disease_specific_survival_status) & disease_specific_survival_status != "NA") |> 
  mutate(disease_specific_survival_status = str_extract(disease_specific_survival_status, "^[^:]+"),
         disease_specific_survival_status = as.numeric(disease_specific_survival_status))
# 1 = dead with tumor
# 0 = alive or dead tumor free

# survival KM curve for cases with SYNJ2 cryptic
survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ 1, 
        data = subset(survival_SYNJ2_cryptic, synj2_cryptic_detected==TRUE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: SYNJ2 Cryptic Events"
  ) +
  add_confidence_interval() 

# survival KM curve for cases with no cryptic
survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ 1, 
        data = subset(survival_SYNJ2_cryptic, synj2_cryptic_detected==FALSE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: No SYNJ2 Cryptic Events"
  ) +
  add_confidence_interval() 

# survival KM curve - both cryptic and non
survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ synj2_cryptic_detected, 
        data = survival_SYNJ2_cryptic) |> 
  ggsurvfit() +
  labs(
    x = "Months with disease",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: All Cancer Types (SYNJ2)"
  ) +
  add_confidence_interval() 


# density plot - months survival in cryptic vs non
survival_SYNJ2_cryptic |> 
  mutate(synj2_cryptic_detected = as.logical(synj2_cryptic_detected)) |> 
  ggplot(aes(x = months_of_disease_specific_survival, colour = synj2_cryptic_detected)) +
  geom_density() +
  labs(
    x = "Survival months with disease"  
  )
# same survival - cryptic vs non-cryptic (all cancers combined)

# number of SYNJ2 detected and non in each cancer type
SYNJ2_n_detected_non_cancer_abbrev <- survival_SYNJ2_cryptic |> 
  group_by(cancer_abbrev, synj2_cryptic_detected) |> 
  summarise(count = n()) |> 
  ungroup() |> 
  pivot_wider(
    names_from = synj2_cryptic_detected,
    values_from = count
  ) |> 
  rename("synj2_cryptic_true" = "TRUE") |> 
  rename("synj2_cryptic_false" = "FALSE")