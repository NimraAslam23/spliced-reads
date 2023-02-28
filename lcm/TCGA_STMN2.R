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

jir_new$gdc_cases.submitter_id
write_clip(jir_new$gdc_cases.submitter_id)

# Import gene data from TCGA ---------------------------------------------------

STMN2_events_genes_orig <- read.csv("STMN2_events_genes.tsv", sep = "\t", header = TRUE, na.strings="", fill = TRUE)
head(STMN2_events_genes_orig)

STMN2_events_genes <- STMN2_events_genes_orig

STMN2_events_genes <- STMN2_events_genes |> 
  rename("affected_cases_in_cohort" = "X..SSM.Affected.Cases.in.Cohort") |> 
    # number of cases where gene is mutated / number of cases tested for simple somatic mutations (3,259)
  rename("affected_cases_in_GDC" = "X..SSM.Affected.Cases.Across.the.GDC") |> 
    # number of cases where gene contains simple somatic mutations / number of cases tested for simple somatic mutations portal wide (13,714)
  rename("CNV_gain" = "X..CNV.Gain") |> 
    # number of cases where CNV gain observed in gene / number of cases tested for copy number alteration in gene (3,344)
  rename("CNV_loss" = "X..CNV.Loss") |> 
    # number of cases where CNV loss observed in gene / number of cases tested for copy number alteration in gene (3,344)
  rename("mutations" = "X..Mutations")
    # unique simple somatic mutations in the gene in cohort 

# Import mutation data from TCGA ------------------------------------------

STMN2_events_mutations_orig <- read.csv("STMN2_events_mutations.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)
head(STMN2_events_mutations_orig)

STMN2_events_mutations <- STMN2_events_mutations_orig

# Import clinical data from TCGA ------------------------------------------

STMN2_clinical_orig <- read.csv("STMN2_clinical.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)
head(STMN2_clinical_orig)

STMN2_clinical <- STMN2_clinical_orig

STMN2_clinical <- STMN2_clinical|> 
  select(case_submitter_id, project_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage) |> 
  separate(project_id, into = c("project", "cancer_type")) |> 
  select(-project)

# bar plot - cancer types with cryptic STMN2 events -----------------------

STMN2_clinical |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_type)))) +
  geom_bar(aes(fill = cancer_type)) +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = "STMN2 events are found in mostly breast and brain cancers" 
  ) +
  theme(legend.position = "none", plot.title = element_text(size=10)) 

# joining jir table with clinical table -----------------------------------

jir_new <- jir_new |> 
  rename("case_submitter_id" = "gdc_cases.submitter_id")

STMN2_clinical_jir <- jir_new |> 
  left_join(STMN2_clinical, by=c("case_submitter_id")) 
STMN2_clinical_jir <- STMN2_clinical_jir |> 
  select(-gdc_cases.diagnoses.tumor_stage, -age_at_index, -gdc_cases.demographic.gender, 
         -ajcc_pathologic_stage, -ethnicity, -race, -cancer_type, -gender)
 
STMN2_clinical_jir <- STMN2_clinical_jir |> 
  rename("cancer_type" = "gdc_cases.project.name") |> 
  rename("gdc_primary_site" = "gdc_cases.project.primary_site") |> 
  rename("sample_type" = "gdc_cases.samples.sample_type") |> 
  rename("cgc_primary_site" = "cgc_case_primary_site")

# filter for cryptic coverage > 2 -----------------------------------------

STMN2_clinical_jir_cryptic <- STMN2_clinical_jir |> 
  filter(STMN2_cryptic_coverage > 2)

# bar plot - primary sites of cancers with cryptic STMN2 events -----------

STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(gdc_primary_site)))) +
  geom_bar(aes(fill = gdc_primary_site)) +
  coord_flip() +
  labs(
    x = "Primary Site of Cancer",
    y = "Number of Cases",
    title = 
      "Primary sites of cancers with cryptic STMN2 events are mostly the adrenal gland and brain" 
  ) +
  theme(legend.position = "none", plot.title = element_text(size=9)) 

# boxplot - cryptic STMN2 coverage in different cancer sites --------------

STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  filter(cgc_primary_site != "") |> 
  ggplot(aes(x = junction_avg_coverage, y = fct_reorder(gdc_primary_site, 
                                                        junction_avg_coverage, median))) +
  geom_boxplot(aes(fill = gdc_primary_site)) +
  labs(
    x = "Junction Average Coverage",
    y = "Primary Site of Cancer",
  ) +
  theme(legend.position = "none", plot.title = element_text(size=10)) +
  geom_signif(comparisons = list(c("Stomach", "Breast"), c("Stomach", "Brain")),
              map_signif_level = TRUE,
              y_position = c(75, 80))

# # / fraction of cases with STMN2 events in cancer types -----------------

STMN2_events_different_cancers <- STMN2_clinical_jir_cryptic |> 
  drop_na() |>
  janitor::tabyl(cancer_type) |> arrange(-percent)

# # / fraction of cases with STMN2 events in cancer sites -----------------

STMN2_events_primary_sites1 <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(gdc_primary_site) |> arrange(-percent)

# All clinical data from TCGA ---------------------------------------------

TCGA_all_clinical_orig <- read.csv("TCGA_all_clinical.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)
head(TCGA_all_clinical_orig)

TCGA_clinical <- TCGA_all_clinical_orig

TCGA_clinical <- TCGA_clinical|> 
  select(case_submitter_id, project_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage) |> 
  separate(project_id, into = c("project", "cancer_type")) |> 
  select(-project)



