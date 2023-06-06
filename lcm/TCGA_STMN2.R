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
library(ggsurvfit)
library(survival)
library(survminer)

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

# boxplot - overall STMN2 junction coverage in different cancer sites --------------

STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
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

# boxplot - cryptic STMN2 junction coverage in different cancer sites ------------

STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = STMN2_cryptic_coverage, 
             y = fct_reorder(gdc_primary_site, STMN2_cryptic_coverage, median))) +
  geom_boxplot(aes(fill = gdc_primary_site)) +
  labs(
    x = "Cryptic Coverage",
    y = "Primary Site of Cancer"
  ) +
  theme(legend.position = "none", plot.title = element_text(size=10))

  # adding reads per million (rpm) column for coverage

STMN2_clinical_jir_cryptic <- STMN2_clinical_jir_cryptic |> 
  mutate(rpm = (STMN2_cryptic_coverage/junction_coverage)*1000000)

STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = rpm,
             y = fct_reorder(gdc_primary_site, rpm, median))) +
  geom_boxplot(aes(fill = gdc_primary_site)) +
  labs(
    x = "reads per million",
    y = "Primary Site of Cancer"
  ) +
  theme(legend.position = "none")
  
# num / fraction of cases with STMN2 events in cancer types -----------------

STMN2_events_different_cancers <- STMN2_clinical_jir_cryptic |> 
  drop_na() |>
  janitor::tabyl(cancer) |> arrange(-percent)

STMN2_events_different_cancers |> 
  mutate(percentage = percent*100) |> 
  ggplot(aes(x = reorder(cancer_type, percentage), y = percentage)) +
  geom_bar(aes(fill = cancer_type), stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Percentage of Cases"
  ) +
  theme(legend.position = "none")

# num / fraction of cases with STMN2 events in cancer sites -----------------

STMN2_events_primary_sites1 <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(gdc_primary_site) |> arrange(-percent)

STMN2_events_primary_sites1 |> 
  mutate(percentage = percent*100) |> 
  ggplot(aes(x = reorder(gdc_primary_site, percentage), y = percentage)) +
  geom_bar(aes(fill = gdc_primary_site), stat = "identity") +
  coord_flip() +
  labs(
    x = "Primary Sites of Cancer",
    y = "Percentage of Cases",
    title = "Cancers with cryptic STMN2 events are primarily in the adrenal gland and brain"
  ) +
  theme(legend.position = "none", plot.title = element_text(size=8))
# 46% of cancers with cryptic STMN2 events are in the adrenal gland

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
  left_join(cBio_clinical_2, by = "case_submitter_id") |> 
  distinct()

STMN2_cryptic_cBio <- STMN2_cryptic_cBio |> 
  select(-cgc_primary_site) |> 
  rename("disease_survival_months" = "Months.of.disease.specific.survival") |> 
  relocate(case_submitter_id, .after = sample_id) |> 
  relocate(cancer, .after = case_submitter_id) |> 
  relocate(gdc_primary_site, .after = cancer_abbrev) |> 
  relocate(sample_type, .after = gdc_primary_site) |> 
  janitor::clean_names()

# mutation counts per cancer type ---------------------------------------------------------

STMN2_cryptic_cBio_mutations <- STMN2_cryptic_cBio |>
  drop_na() |> 
  filter(mutation_count !="NA") |> 
  group_by(cancer) |> 
  mutate(mean_mutation_count = mean(as.numeric(mutation_count))) |> 
  ungroup() |> 
  select(cancer, cancer_abbrev, mean_mutation_count) |> 
  unique() 

STMN2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 2500) |> 
  ggplot(aes(x = fct_reorder(cancer_abbrev, mutation_count, median), y = mutation_count)) +
  labs(
    x = "Cancer Type",
    y = "Mutation Count",
    title = "UCEC and LUSC cancers have the greatest number of mutations"
  #among the cancers with cryptic STMN2 events
  ) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("UCEC", "LUSC")), test = "wilcox.test",
              map_signif_level = TRUE)


# fraction of mutations in each cancer type that is found in cases with cryptic --------

total_mutations_each_cancer <- cBio_clinical |> 
  drop_na(Mutation.Count) |> 
  filter(grepl("^\\d+$", Mutation.Count)) |>
  group_by(cancer_abbrev) |> 
  summarise(total_mutations = sum(as.numeric(Mutation.Count)))

total_mutations_each_cancer_with_cryptic <- STMN2_cryptic_cBio |> 
  drop_na(mutation_count) |> 
  filter(mutation_count < 2500) |> 
  group_by(cancer_abbrev) |> 
  summarise(total_mutations_cryptic = sum(mutation_count))

mutations_each_cancer_general_vs_cryptic <- total_mutations_each_cancer |> 
  left_join(total_mutations_each_cancer_with_cryptic, by=c("cancer_abbrev")) |> 
  mutate(percent_with_cryptic = (total_mutations_cryptic/total_mutations)*100)
  #16% of all mutations in LGG cancer are in cases with STMN2 cryptic events

#calculating mutation burden (weighting total mutations against total cases) for general vs cryptic

mutations_each_cancer_general_vs_cryptic <- mutations_each_cancer_general_vs_cryptic |>
  mutate(total_cancer_abbrev = total_each_cancer_general_vs_cryptic$total_cancer_abbrev) |> 
  mutate(total_cancer_abbrev_with_cryptic = total_each_cancer_general_vs_cryptic$total_cancer_abbrev_with_cryptic) |> 
  mutate(avg_mutation_per_case = total_mutations/total_cancer_abbrev) |> 
  mutate(avg_mutation_per_case_cryptic = total_mutations_cryptic/total_cancer_abbrev_with_cryptic) |> 
  select(-total_cancer_abbrev, -total_cancer_abbrev_with_cryptic)

mutations_each_cancer_general_vs_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = reorder(cancer_abbrev, avg_mutation_per_case_cryptic))) +
  geom_bar(aes(y = avg_mutation_per_case, fill = "avg_mutation_per_case"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = avg_mutation_per_case_cryptic, fill = "avg_mutation_per_case_cryptic"), stat = "identity", position = "dodge") +
  labs(
    x = "Cancer Type",
    y = "Average number of mutations per case",
    fill = ""
  ) +
  scale_fill_manual(values = c("avg_mutation_per_case" = "violetred4", "avg_mutation_per_case_cryptic" = "royalblue1"),
                    labels = c("All Cases", "STMN2 Cryptic Cases only"))

cBio_clinical |>
  left_join(STMN2_clinical_jir,by = c("case_submitter_id")) |> 
  select('case_submitter_id',Study.ID,cancer_abbrev,Mutation.Count,STMN2_cryptic_coverage,STMN2_annotated_coverage) |> 
  separate(Study.ID,into = ('study_start'),remove = FALSE) |> 
  unique() |> 
  mutate(Mutation.Count = as.numeric(Mutation.Count)) |> 
  filter(!is.na(Mutation.Count) & Mutation.Count != "NA") |> 
  mutate(stmn2_cryptic_detected = STMN2_cryptic_coverage >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_stmn2 = sum(stmn2_cryptic_detected,na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(stmn2_cryptic_detected)) |> 
  filter(n_detected_stmn2 > 2) |> 
  ggplot(aes(x = study_start,
             y = Mutation.Count)) + 
  geom_boxplot(aes(fill = stmn2_cryptic_detected)) + 
  labs(x = "Cancer Type",
       y = "Mutation Count") +
  coord_flip() + 
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.format")


# survival comparisons ----------------------------------------------------

survival_STMN2_cryptic <- cBio_clinical |> 
  left_join(STMN2_clinical_jir, by = c("case_submitter_id")) |> 
  select(Study.ID, case_submitter_id, Sample.ID, cancer_abbrev, Months.of.disease.specific.survival, Mutation.Count, 
         Overall.Survival..Months., STMN2_cryptic_coverage, STMN2_annotated_coverage) |>
  mutate(Months.of.disease.specific.survival = as.numeric(Months.of.disease.specific.survival)) |> 
  filter(!is.na(Months.of.disease.specific.survival) & Months.of.disease.specific.survival != "NA") |> 
  mutate(Mutation.Count = as.numeric(Mutation.Count)) |> 
  filter(!is.na(Mutation.Count) & Mutation.Count != "NA") |> 
  mutate(stmn2_cryptic_detected = STMN2_cryptic_coverage >= 2) |> 
  group_by(cancer_abbrev) |> 
  mutate(n_detected_stmn2 = sum(stmn2_cryptic_detected,na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(stmn2_cryptic_detected)) |> 
  filter(n_detected_stmn2 > 2) |> 
  mutate(stmn2_cryptic_detected = as.logical(stmn2_cryptic_detected)) |> #changes FALSE and TRUE to 0 and 1 respectively
  distinct()
  
survfit(Surv(Months.of.disease.specific.survival, stmn2_cryptic_detected) ~ 1, 
        data = subset(survival_STMN2_cryptic, stmn2_cryptic_detected==TRUE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability"
  ) +
  add_confidence_interval() 

survfit(Surv(Months.of.disease.specific.survival, stmn2_cryptic_detected) ~ 1, 
        data = subset(survival_STMN2_cryptic, stmn2_cryptic_detected==FALSE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability"
  ) +
  add_confidence_interval() 

survfit(Surv(Months.of.disease.specific.survival, stmn2_cryptic_detected) ~ stmn2_cryptic_detected, 
        data = survival_STMN2_cryptic) |> 
  ggsurvfit() +
  labs(
    x = "Months with disease",
    y = "Overall Survival Probability",
    title = "All Cancer Types"
  ) +
  add_confidence_interval() 

# density plot - months survival in cryptic vs non-cryptic ----------------

survival_STMN2_cryptic |> 
  mutate(stmn2_cryptic_detected = as.logical(stmn2_cryptic_detected)) |> 
  ggplot(aes(x = Months.of.disease.specific.survival, colour = stmn2_cryptic_detected)) +
  geom_density() +
  labs(
    x = "Survival months with disease"  
    )
    # same survival - cryptic vs non-cryptic (all cancers combined)


# number of STMN2 detected and non in each cancer type --------------------

n_detected_non_cancer_abbrev <- survival_STMN2_cryptic |> 
  group_by(cancer_abbrev, stmn2_cryptic_detected) |> 
  summarise(count = n()) |> 
  ungroup() |> 
  pivot_wider(
    names_from = stmn2_cryptic_detected,
    values_from = count
  ) |> 
  rename("stmn2_cryptic_true" = "TRUE") |> 
  rename("stmn2_cryptic_false" = "FALSE")
      # PCPG has 130 false, 51 true
      # GBM has 154 false, 11 true
  

# PCPG survival analysis --------------------------------------------------

survfit(Surv(Months.of.disease.specific.survival, stmn2_cryptic_detected) ~ stmn2_cryptic_detected, 
        data = subset(survival_STMN2_cryptic, cancer_abbrev == "PCPG")) |> 
  ggsurvfit() +
  labs(
    x = "Months with disease",
    y = "Overall Survival Probability",
    title = "PCPG Survival"
  ) +
  add_confidence_interval() 


# GBM survival analysis ---------------------------------------------------

survfit(Surv(Months.of.disease.specific.survival, stmn2_cryptic_detected) ~ stmn2_cryptic_detected, 
        data = subset(survival_STMN2_cryptic, cancer_abbrev == "GBM")) |> 
  ggsurvfit() +
  labs(
    x = "Months with disease",
    y = "Overall Survival Probability",
    title = "GBM Survival"
  ) +
  add_confidence_interval() 

# fraction of each cancer that has cryptic STMN2 events -------------------

total_each_cancer <- cBio_clinical |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev = n())

total_each_cancer_with_cryptic <- STMN2_cryptic_cBio |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev_with_cryptic = n())

total_each_cancer_general_vs_cryptic <- total_each_cancer |> 
  left_join(total_each_cancer_with_cryptic, by=c("cancer_abbrev")) |> 
  mutate(percent_with_cryptic = (total_cancer_abbrev_with_cryptic/total_cancer_abbrev)*100)

total_each_cancer_general_vs_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent_with_cryptic), y = percent_with_cryptic)) +
  geom_bar(aes(fill = cancer_abbrev), stat = "identity") +
  labs(
    x = "Cancer Type",
    y = "Percentage of cases with cryptic STMN2 events"
    ) +
  theme(legend.position = "none")

# mutation data of one patient --------------------------------------------

patient_mutations_orig <- read.csv("TCGA-06-5416_mutations.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE) 

patient_mutations <- patient_mutations_orig |> 
  select(Gene, Protein.Change, Annotation, Functional.Impact, Chromosome, Start.Pos, End.Pos, Ref, Var, 
         HGVSg, Mutation.Type, Variant.Type, Allele.Freq)

    # number of mutations in each gene 

patient_mutations |> 
  count(Gene, sort = TRUE) |> 
  filter(n > 10) |> 
  ggplot(aes(x = fct_reorder(Gene, n, mean), y = n)) +
  geom_bar(stat = "identity", aes(fill = Gene)) +
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
                                                        
  # number of mutations in each cancer driver gene 

cosmic_patient_mutations |> 
  count(Gene, sort = TRUE) |> 
  filter(n > 5) |> 
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


