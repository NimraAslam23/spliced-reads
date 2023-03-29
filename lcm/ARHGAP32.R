# gene_name = "ARHGAP32" 
# snapcount_coords_cryptic = "chr11:128992047-128998318"
# snapcount_coords_annotated = "chr11:128988126-128998318"
# strand_code = "-"

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

ARHGAP32_cryptic_query <- count_query(gene_name = "ARHGAP32", snapcount_coords = "chr11:128992047-128998318", strand_code = "-")
ARHGAP32_anno_query <- count_query(gene_name = "ARHGAP32", snapcount_coords = "chr11:128988126-128998318", strand_code = "-")

# join cryptic and anno df, add jir column --------------------------------

ARHGAP32_cryptic_query <- ARHGAP32_cryptic_query |> rename("cryptic_count" = "count")
ARHGAP32_anno_query <- ARHGAP32_anno_query |> rename("anno_count" = "count")

ARHGAP32_query <- ARHGAP32_cryptic_query |> 
  select(gdc_cases.submitter_id, cryptic_count) |> 
  left_join(ARHGAP32_anno_query, by = "gdc_cases.submitter_id")

ARHGAP32_query <- ARHGAP32_query |> 
  relocate(cryptic_count, .after = anno_count)

ARHGAP32_query <- ARHGAP32_query |> 
  mutate(jir = anno_count/(anno_count + cryptic_count)) |> 
  relocate(jir, .after = cryptic_count)

write_clip(ARHGAP32_query$gdc_cases.submitter_id)

# clinical data from TCGA (case set) --------------------------------------

ARHGAP32_clinical_orig <- read.csv("ARHGAP32_clinical.tsv", sep = "\t", header = TRUE,
                                   na.strings = "", fill = TRUE)

ARHGAP32_clinical <- ARHGAP32_clinical_orig

ARHGAP32_clinical <- ARHGAP32_clinical|> 
  select(case_submitter_id, project_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage) |> 
  separate(project_id, into = c("project", "cancer_type")) |> 
  select(-project)

# bar plot - cancer types with ARHGAP32 events ----------------------------

ARHGAP32_clinical |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_type)))) +
  geom_bar() +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = "ARHGAP32 is expressed in most breast cancers"
  ) +
  theme(plot.title = element_text(size=10))

# join query table with clinical table ------------------------------------

ARHGAP32_clinical_jir <- ARHGAP32_query |> rename("case_submitter_id" = "gdc_cases.submitter_id") |> 
  left_join(ARHGAP32_clinical, by=c("case_submitter_id")) |> 
  select(-gdc_cases.diagnoses.tumor_stage, -age_at_index, -gdc_cases.demographic.gender, 
         -ajcc_pathologic_stage, -ethnicity, -race, -cancer_type, -gender) |> 
  rename("cancer_type" = "gdc_cases.project.name") |> 
  rename("gdc_primary_site" = "gdc_cases.project.primary_site") |> 
  rename("sample_type" = "gdc_cases.samples.sample_type") |> 
  rename("cgc_primary_site" = "cgc_case_primary_site")

ARHGAP32_clinical_jir_cryptic <- ARHGAP32_clinical_jir |> 
  filter(cryptic_count > 5) 

# primary sites of cancers with cryptic ARHGAP32 events -------------------

ARHGAP32_clinical_jir_cryptic |> 
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

# cryptic ARHGAP32 expression in different cancer sites -------------------

ARHGAP32_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = cryptic_count,
             y = fct_reorder(gdc_primary_site, cryptic_count, median))) +
  geom_boxplot() +
  labs(
    x = "ARHGAP32 Cryptic Coverage",
    y = "Primary Site of Cancer"
  ) +
  theme(plot.title = element_text(size=10))

# which cancers have the most cryptic ARHGAP32 events? --------------------

ARHGAP32_events_different_cancers <- ARHGAP32_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(cancer_type) |> arrange(-percent)

# where are the cancers with cryptic ARHGAP32 events located? -------------

ARHGAP32_events_primary_sites <- ARHGAP32_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(gdc_primary_site) |> arrange(-percent)

# cBioPortal clinical data, join with ARHGAP32 table ----------------------

ARHGAP32_cryptic_cBio <- ARHGAP32_clinical_jir_cryptic |> 
  left_join(cBio_clinical, by = "case_submitter_id") |> 
  select(-cgc_primary_site, -Cancer.Type, -Overall.Survival..Months.) |>
  rename("disease_survival_months" = "Months.of.disease.specific.survival") |> 
  relocate(case_submitter_id, .after = Sample.ID) |> 
  relocate(cancer_type, .after = case_submitter_id) |> 
  relocate(gdc_primary_site, .after = cancer_abbrev) |> 
  relocate(sample_type, .after = gdc_primary_site) |> 
  janitor::clean_names()

# boxplot - n mutations for each cancer type ------------------------------

ARHGAP32_cryptic_cBio <- ARHGAP32_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric)

ARHGAP32_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 2000) |> 
  ggplot(aes(x = fct_reorder(cancer_abbrev, mutation_count, median),
             y = mutation_count)) +
  geom_boxplot() +
  labs(
    x = "Cancer Type",
    y = "Mutation Count",
    title = ""
  ) 

# cryptic coverage vs mutation count --------------------------------------

ARHGAP32_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 600) |> 
  ggplot(aes(x = as.factor(cryptic_count),
             y = mutation_count)) +
  geom_boxplot() +
  labs(
    x = "Number of ARHGAP32 cryptic events",
    y = "Mutation Count"
  )

ARHGAP32_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 500) |> 
  ggplot(aes(x = cryptic_count, 
             y = mutation_count)) +
  geom_hex() +
  labs(
    x = "Number of ARHGAP32 cryptic events",
    y = "Mutation Count"
  )

