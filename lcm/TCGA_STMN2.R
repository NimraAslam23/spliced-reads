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
library(TCGAbiolinks)

jir_new <- read.csv("jir_new.csv")

jir_new$case_submitter_id
write_clip(jir_new$case_submitter_id)

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
    title = "STMN2 is expressed in mostly breast and brain cancers" 
  ) +
  theme(legend.position = "none", plot.title = element_text(size=10)) 

# joining jir table with clinical table -----------------------------------

STMN2_clinical_jir <- jir_new |> 
  left_join(STMN2_clinical, by=c("case_submitter_id")) 
STMN2_clinical_jir <- STMN2_clinical_jir |> 
  select(-gdc_cases.diagnoses.tumor_stage, -age_at_index, -gdc_cases.demographic.gender, 
         -ajcc_pathologic_stage, -ethnicity, -race, -gender)
 
STMN2_clinical_jir <- STMN2_clinical_jir |> 
  rename("cancer" = "gdc_cases.project.name") |> 
  rename("gdc_primary_site" = "gdc_cases.project.primary_site") |> 
  rename("sample_type" = "gdc_cases.samples.sample_type") |> 
  rename("cgc_primary_site" = "cgc_case_primary_site") |> 
  rename("cancer_abbrev" = "cancer_type")

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

# boxplot - overall STMN2 junction coverage in different cancer sites --------------

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

# boxplot - cryptic STMN2 expression in different cancer sites ------------

STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  filter(cgc_primary_site != "") |> 
  ggplot(aes(x = STMN2_cryptic_coverage, 
             y = fct_reorder(gdc_primary_site, STMN2_cryptic_coverage, median))) +
  geom_boxplot(aes(fill = gdc_primary_site)) +
  labs(
    x = "Junction Average Coverage",
    y = "Primary Site of Cancer",
  ) +
  theme(legend.position = "none", plot.title = element_text(size=10))

# num / fraction of cases with STMN2 events in cancer types -----------------

STMN2_events_different_cancers <- STMN2_clinical_jir_cryptic |> 
  drop_na() |>
  janitor::tabyl(cancer) |> arrange(-percent)

# num / fraction of cases with STMN2 events in cancer sites -----------------

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

# All clinical data from cBioPortal ---------------------------------------

cBio_all_clinical_orig <- read.csv("cBio_clinical_data.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)

cBio_clinical <- cBio_all_clinical_orig

cBio_clinical <- cBio_clinical |> 
  select(Study.ID, Patient.ID, Sample.ID, Aneuploidy.Score, Cancer.Type, TCGA.PanCanAtlas.Cancer.Type.Acronym, Cancer.Type.Detailed,
         Months.of.disease.specific.survival, Mutation.Count, Overall.Survival..Months.)

cBio_clinical <- cBio_clinical |> 
  rename("case_submitter_id" = "Patient.ID") |> 
  rename("cancer" = "Cancer.Type.Detailed") |> 
  rename("cancer_abbrev" = "TCGA.PanCanAtlas.Cancer.Type.Acronym")

# join cBio data with STMN2 cryptic table ---------------------------------

cBio_clinical_2 <- cBio_clinical |> 
  select(-Cancer.Type, -cancer, -Sample.ID, -Overall.Survival..Months., -cancer_abbrev)

STMN2_cryptic_cBio <- STMN2_clinical_jir_cryptic |> 
  left_join(cBio_clinical_2, by = "case_submitter_id") 

STMN2_cryptic_cBio <- STMN2_cryptic_cBio |> 
  select(-cgc_primary_site) |> 
  rename("disease_survival_months" = "Months.of.disease.specific.survival") |> 
  relocate(case_submitter_id, .after = sample_id) |> 
  relocate(cancer, .after = case_submitter_id) |> 
  relocate(gdc_primary_site, .after = cancer_abbrev) |> 
  relocate(sample_type, .after = gdc_primary_site) |> 
  janitor::clean_names()

# mutation counts ---------------------------------------------------------

STMN2_cryptic_cBio_mutations <- STMN2_cryptic_cBio |>
  drop_na() |> 
  filter(mutation_count !="NA") |> 
  group_by(cancer) |> 
  mutate(mean_mutation_count = mean(as.numeric(mutation_count))) |> 
  ungroup() |> 
  select(cancer, cancer_abbrev, mean_mutation_count) |> 
  unique() 

STMN2_cryptic_cBio_mutations |> 
  ggplot(aes(x = fct_reorder(cancer_abbrev, mean_mutation_count, median), y = mean_mutation_count)) +
  geom_bar(stat = 'identity', aes(fill = cancer_abbrev)) +
  labs(
    x = "Cancer Type",
    y = "Mean Mutation Count"
  ) +
  theme(legend.position = "none")

# mutation data of one patient --------------------------------------------

patient_mutations_orig <- read.csv("TCGA-A5-A0G2_mutations.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE) 

patient_mutations <- patient_mutations_orig |> 
  select(Gene, Protein.Change, Annotation, Functional.Impact, Chromosome, Start.Pos, End.Pos, Ref, Var, 
         HGVSg, Mutation.Type, Variant.Type, Allele.Freq)

    # number of mutations in each gene 

patient_mutations |> 
  count(Gene, sort = TRUE) |> 
  filter(n > 17) |> 
  ggplot(aes(x = fct_reorder(Gene, n, mean), y = n)) +
  geom_bar(stat = 'identity', aes(fill = Gene)) +
  labs(
    x = "Gene",
    y = "Mutation Count"
  ) +
  theme(legend.position = "none")

    # types of mutations in TTN

patient_mutations |> 
  filter(Gene == "TTN") |> 
  ggplot(aes(x = fct_rev(fct_infreq(Mutation.Type)))) +
  geom_bar(aes(fill = Mutation.Type)) +
  labs(
    x = "Mutation Type",
    y = "Count in TTN Gene"
  ) +
  theme(legend.position = "none")


# STMN2 cryptic coverage vs mutation count --------------------------------

STMN2_cryptic_cBio <- STMN2_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric) 

STMN2_cryptic_cBio |> 
  drop_na() |> 
  ggplot(aes(x = stmn2_cryptic_coverage, y = mutation_count)) +
  labs(
    x = "Number of STMN2 cryptic events",
    y = "Mutation Count"
  ) +
  geom_point()

STMN2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 400, stmn2_cryptic_coverage > 10) |> 
  ggplot(aes(x = as.factor(stmn2_cryptic_coverage), y = mutation_count)) +
  labs(
    x = "Number of STMN2 cryptic events",
    y = "Mutation Count"
  ) +
  geom_boxplot()

# COSMIC Cancer Gene Consensus --------------------------------------------

cosmic_cancer_genes <- read.csv("cosmic_cancer_gene_consensus.csv") |> 
  rename("Gene" = "Gene.Symbol")

# filter for cancer driver genes ------------------------------------------

cosmic_patient_mutations <- filter(patient_mutations, Gene %in% cosmic_cancer_genes$Gene)
                                                        
  # number of mutations in each gene 

cosmic_patient_mutations |> 
  count(Gene, sort = TRUE) |> 
  filter(n > 11) |> 
  ggplot(aes(x = fct_reorder(Gene, n, mean), y = n)) +
  geom_bar(stat = 'identity', aes(fill = Gene)) +
  labs(
    x = "Cancer Gene",
    y = "Mutation Count"
  ) +
  theme(legend.position = "none")

  # types of mutations in MUC16

cosmic_patient_mutations |> 
  filter(Gene == "MUC16") |> 
  ggplot(aes(x = fct_rev(fct_infreq(Mutation.Type)))) +
  geom_bar(aes(fill = Mutation.Type)) +
  labs(
    x = "Mutation Type",
    y = "Count in MUC16 Gene"
  ) +
  theme(legend.position = "none")

# type of cancer genes ----------------------------------------------------

cosmic_patient_mutations <- cosmic_patient_mutations |> 
  left_join(cosmic_cancer_genes, by = "Gene") 

cosmic_patient_mutations <- cosmic_patient_mutations |> 
  select(-Entrez.GeneId, -Tier, -Hallmark, -Chr.Band, -Somatic, -Germline, -Tumour.Types.Somatic., 
         -Tumour.Types.Germline., -Cancer.Syndrome, -Molecular.Genetics, -Translocation.Partner, 
         -Other.Germline.Mut, -Other.Syndrome, -Synonyms) 

  # fusion
      # Gene fusions, or translocations, resulting from chromosomal rearrangements are the most common mutation class. 
      #They lead to chimeric transcripts or to deregulation of genes through juxtapositioning of novel promoter or enhancer regions.

  # TSG = tumour suppressor gene

  # oncogene

cosmic_patient_mutations |> 
  drop_na() |> 
  filter(Role.in.Cancer !="") |> 
  ggplot(aes(x = Role.in.Cancer)) +
  geom_bar()

cosmic_patient_mutations <- cosmic_patient_mutations |> 
  mutate(TSG = ifelse(grepl("TSG", Role.in.Cancer), "yes", "no")) |> 
  mutate(oncogene = ifelse(grepl("oncogene", Role.in.Cancer), "yes", "no")) |> 
  mutate(fusion = ifelse(grepl("fusion", Role.in.Cancer), "yes", "no"))
       
n_TSG <- cosmic_patient_mutations |> 
  janitor::tabyl(TSG)

n_oncogene <- cosmic_patient_mutations |> 
  janitor::tabyl(oncogene)

n_fusion <- cosmic_patient_mutations |> 
  janitor::tabyl(fusion)
  
# TCGA biolinks -----------------------------------------------------------

library(TCGAbiolinks)

    # all cancer abbreviations

cancer_abbrev <- unique(STMN2_cryptic_cBio$cancer_abbrev, na.rm = TRUE)
cancer_abbrev


