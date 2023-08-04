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

pdf(file = "figures_SYNJ2_analysis.pdf")

# functions
count_query <- function(gene_name, snapcount_coords, strand_code) {
  cryptic_query <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords) 
  
  if(strand_code == "+"){
    cryptic_query <- set_row_filters(cryptic_query, strand == "+") 
  }else if(strand_code == "-"){
    cryptic_query <- set_row_filters(cryptic_query, strand == "-") 
  }
  
  cryptic_query_exact <- set_coordinate_modifier(cryptic_query, Coordinates$Exact) 
  print("junction querying beginning")
  juncs_on <- query_jx(cryptic_query) 
  juncs_on_flat <- query_jx(cryptic_query,return_rse = FALSE) 
  samples_with <- juncs_on@colData |> 
    as.data.frame()
  print("junction querying done")
  samples_with <- samples_with |> 
    select(c("rail_id",
             "gdc_cases.demographic.gender",
             "gdc_cases.submitter_id",
             "gdc_cases.project.name",
             "gdc_cases.project.primary_site",
             "gdc_cases.diagnoses.tumor_stage",
             "gdc_cases.samples.sample_type",
             "cgc_case_primary_site",
             "junction_coverage",
             "junction_avg_coverage")) 
  
  juncs_on_flat <- juncs_on_flat |> 
    select(chromosome,start,end,strand,samples) |> 
    separate_rows(samples, sep = ',') |> 
    filter(samples != "") |> 
    separate(samples, into = c("rail_id","count")) |> 
    mutate(rail_id = as.integer(rail_id)) |> 
    mutate(count = as.numeric(count)) 
  
  gene_df <- juncs_on_flat |> 
    left_join(samples_with) |> 
    filter(end == str_split(snapcount_coords,'-',simplify = TRUE)[[2]]) |> 
    filter(start == str_split(str_split(snapcount_coords,'-',simplify = TRUE)[[1]],":",simplify = TRUE)[[2]])
  return(gene_df)
}

# query cryptic + anno and join into one df -------------------------------

combine_two_junctions <- function(gene_name, snapcount_coords_cryptic, 
                                  snapcount_coords_annotated, strand_code) {
  
  cryptic_query <- count_query(gene_name, snapcount_coords_cryptic, strand_code) |> 
    rename("cryptic_count" = "count") |> 
    select(gdc_cases.submitter_id, cryptic_count)
  anno_query <- count_query(gene_name, snapcount_coords_annotated, strand_code) |> 
    rename("anno_count" = "count")
  
  query <- anno_query |> 
    left_join(cryptic_query, by = "gdc_cases.submitter_id") |> 
    mutate(jir = cryptic_count/(anno_count + cryptic_count)) |> 
    relocate(cryptic_count, .after = anno_count) |> 
    relocate(jir, .after = cryptic_count) |> 
    mutate(gene_name = gene_name) |> 
    mutate(coords_cryptic = snapcount_coords_cryptic)
  
  return(query)
  
}



# SYNJ2 cryptic and annotated reads
gene_name = "SYNJ2" 
snapcount_coords_cryptic = "chr6:158017291-158019983"
snapcount_coords_annotated = "chr6:158017291-158028755"
strand_code = "+"

# query - cryptic and anno ------------------------------------------------

# querying for cryptic and annotated STMN2 reads, join cryptic and anno df, add jir column
SYNJ2_query <- combine_two_junctions("SYNJ2", "chr6:158017291-158019983", "chr6:158017291-158028755", "+")
write.table(SYNJ2_query, file="SYNJ2_query.txt", sep=",")
SYNJ2_query <- read.table("SYNJ2_query.txt", sep=",")

#SYNJ2_cryptic_query <- count_query(gene_name = "SYNJ2", snapcount_coords = "chr6:158017291-158019983", strand_code = "+")
#SYNJ2_anno_query <- count_query(gene_name = "SYNJ2", snapcount_coords = "chr6:158017291-158028755", strand_code = "+")

write_clip(SYNJ2_query$gdc_cases.submitter_id)

# clinical data from TCGA (case set) --------------------------------------

SYNJ2_clinical_orig <- read.csv("SYNJ2_clinical.tsv", sep = "\t", header = TRUE,
                                   na.strings = "", fill = TRUE)

SYNJ2_clinical <- SYNJ2_clinical_orig

SYNJ2_clinical <- SYNJ2_clinical|> 
  select(case_submitter_id, project_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage) |> 
  separate(project_id, into = c("project", "cancer_type")) |> 
  select(-project)


# join query table with clinical table ------------------------------------

SYNJ2_clinical_jir <- SYNJ2_query |> rename("case_submitter_id" = "gdc_cases.submitter_id") |> 
  left_join(SYNJ2_clinical, by=c("case_submitter_id")) |> 
  select(-gdc_cases.diagnoses.tumor_stage, -age_at_index, -gdc_cases.demographic.gender, 
         -ajcc_pathologic_stage, -ethnicity, -race, -gender) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_type") |> 
  rename("cancer_type" = "gdc_cases_project_name")

SYNJ2_clinical_jir_cryptic <- SYNJ2_clinical_jir |> 
  filter(cryptic_count > 2) 

# bar plot - cancer types with SYNJ2 events ----------------------------

a_SYNJ2_analysis <- SYNJ2_clinical_jir |> 
  distinct() |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_abbrev)),
             fill = cancer_abbrev)) +
  geom_bar() +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = "SYNJ2 Expression in Different Cancer Types"
  ) +
  theme(legend.position = "none")

print(a_SYNJ2_analysis)

# Cancers with cryptic ARHGAP32 events

# bar plot - cancer type with SYNJ2 events, primary sites with cryptic

b_SYNJ2_analysis <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_abbrev)),
             fill = cancer_abbrev)) +
  geom_bar() +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = "Cryptic SYNJ2 Expression in Different Cancer Types"
  ) +
  theme(legend.position = "none")

print(b_SYNJ2_analysis)

c_SYNJ2_analysis <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(gdc_cases_project_primary_site)))) +
  geom_bar(aes(fill = gdc_cases_project_primary_site)) +
  coord_flip() +
  labs(
    x = "Primary Site of Cancer",
    y = "Number of Cases",
    title = "Cryptic SYNJ2 Expression in Different Cancer Sites"
  ) +
  theme(legend.position = "none")

print(c_SYNJ2_analysis)

## Junction coverage

#overall SYNJ2 junction coverage in different cancer sites
d_SYNJ2_analysis <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = junction_avg_coverage, 
             y = fct_reorder(gdc_cases_project_primary_site, junction_avg_coverage, median))) +
  geom_boxplot(aes(fill = gdc_cases_project_primary_site)) +
  labs(
    x = "Junction Average Coverage",
    y = "Primary Site of Cancer",
    title = "Average SYNJ2 Junction Coverage in Different Cancer Sites"
  ) + 
  theme(legend.position = "none") 

print(d_SYNJ2_analysis)

#adding reads per million (rpm) column for cryptic coverage
SYNJ2_clinical_jir_cryptic <- SYNJ2_clinical_jir_cryptic |> 
  mutate(rpm = (cryptic_count/junction_coverage)*1000000)

e_SYNJ2_analysis <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = log10(rpm),
             y = fct_reorder(gdc_cases_project_primary_site, rpm, median))) +
  geom_boxplot(aes(fill = gdc_cases_project_primary_site)) +
  labs(
    x = "Log10 Reads per Million",
    y = "Primary Site of Cancer",
    title = "SYNJ2 Cryptic Coverage in Different Cancer Sites"
  ) +
  theme(legend.position = "none")

print(e_SYNJ2_analysis)

# cryptic SYNJ2 expression in different cancer sites -------------------

f_SYNJ2_analysis <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = cryptic_count,
             y = fct_reorder(gdc_cases_project_primary_site, cryptic_count, median))) +
  geom_boxplot(aes(fill = gdc_cases_project_primary_site)) +
  labs(
    x = "SYNJ2 Cryptic Coverage",
    y = "Primary Site of Cancer"
  ) +
  theme(plot.title = element_text(size=10),
        legend.position = "none")

print(f_SYNJ2_analysis)

# which cancers have the most cryptic SYNJ2 events? --------------------

SYNJ2_events_different_cancers <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(cancer_abbrev) |> arrange(-percent)

#kable(SYNJ2_events_different_cancers, caption = "STAD cancer has high cryptic SYNJ2 expression")

g_SYNJ2_analysis <- SYNJ2_events_different_cancers |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent), 
             y = percent,
             fill = cancer_abbrev)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Percentage of Cases",
    title = "Fraction of all cryptic SYNJ2 cases that are of each cancer type"
  ) +
  theme(legend.position = "none") + 
  scale_y_continuous(labels = scales::percent_format())

print(g_SYNJ2_analysis)

# where are the cancers with cryptic SYNJ2 events located? -------------

SYNJ2_events_primary_sites <- SYNJ2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(gdc_cases_project_primary_site) |> arrange(-percent)

#kable(SYNJ2_events_primary_sites, caption = "Cancers with cryptic SYNJ2 events are found primarily in the stomach.")

h_SYNJ2_analysis <- SYNJ2_events_primary_sites |> 
  ggplot(aes(x = reorder(gdc_cases_project_primary_site, percent), 
             y = percent,
             fill = gdc_cases_project_primary_site)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Primary Sites of Cancer",
    y = "Percentage of Cases",
    title = "Fraction of all cryptic SYNJ2 cases that are in each cancer site"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format())

print(h_SYNJ2_analysis)

# cBioPortal clinical data, join with SYNJ2 table ----------------------

cBio_all_clinical_orig <- read.csv("cBio_clinical_data.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)

cBio_clinical <- cBio_all_clinical_orig

cBio_clinical <- cBio_clinical |> 
  select(Study.ID, Patient.ID, Sample.ID, Aneuploidy.Score, Cancer.Type, TCGA.PanCanAtlas.Cancer.Type.Acronym, Cancer.Type.Detailed,
         Months.of.disease.specific.survival, Mutation.Count, Overall.Survival..Months., Disease.specific.Survival.status)

cBio_clinical <- cBio_clinical |> 
  rename("case_submitter_id" = "Patient.ID") |> 
  rename("cancer" = "Cancer.Type.Detailed") |> 
  rename("cancer_abbrev" = "TCGA.PanCanAtlas.Cancer.Type.Acronym")

SYNJ2_cryptic_cBio <- SYNJ2_clinical_jir_cryptic |> 
  left_join(cBio_clinical, by = "case_submitter_id") |>
  janitor::clean_names() |> 
  select(-c(cgc_case_primary_site, cancer_abbrev_y, cancer_type)) |> 
  rename("cancer_abbrev" = "cancer_abbrev_x")

write.table(SYNJ2_cryptic_cBio, file = "SYNJ2_cryptic_cBio.txt", sep=",")

# boxplot - n mutations for each cancer type ------------------------------

i_SYNJ2_analysis <- SYNJ2_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric) |> 
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

print(i_SYNJ2_analysis)

# cryptic coverage vs mutation count --------------------------------------

j_SYNJ2_analysis <- SYNJ2_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na() |> 
  filter(mutation_count < 300, cryptic_count > 1) |> 
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

print(j_SYNJ2_analysis)

k_SYNJ2_analysis <- SYNJ2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 500) |> 
  ggplot(aes(x = cryptic_count, 
             y = mutation_count)) +
  geom_hex() +
  labs(
    x = "Number of SYNJ2 cryptic events",
    y = "Mutation Count"
  )

print(k_SYNJ2_analysis)


# Fraction of each cancer that has cryptic SYNJ2 events -------------------

total_each_cancer <- cBio_clinical |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev = n())

total_each_cancer_with_SYNJ2cryptic <- SYNJ2_cryptic_cBio |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev_with_cryptic = n())

total_each_cancer_general_vs_SYNJ2cryptic <- total_each_cancer |> 
  left_join(total_each_cancer_with_SYNJ2cryptic, by=c("cancer_abbrev")) |> 
  mutate(percent_with_cryptic = (total_cancer_abbrev_with_cryptic/total_cancer_abbrev)*100)

l_SYNJ2_analysis <- total_each_cancer_general_vs_SYNJ2cryptic |> 
  drop_na() |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent_with_cryptic), y = percent_with_cryptic)) +
  geom_bar(aes(fill = cancer_abbrev), stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Percentage of Cases",
    title = "Percentage of Cases with Cryptic SYNJ2 Events"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format())

print(l_SYNJ2_analysis)

# Mutational burden

# fraction of each cancer that has cryptic SYNJ2 events

total_mutations_each_cancer <- cBio_clinical |> 
  janitor::clean_names() |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na(mutation_count) |> 
  filter(grepl("^\\d+$", mutation_count)) |>
  group_by(cancer_abbrev) |> 
  summarise(total_mutations = sum(mutation_count))

total_mutations_each_cancer_with_SYNJ2cryptic <- SYNJ2_cryptic_cBio |> 
  janitor::clean_names() |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na(mutation_count) |> 
  filter(mutation_count < 2500) |> 
  group_by(cancer_abbrev) |> 
  summarise(total_mutations_cryptic = sum(mutation_count))

mutations_each_cancer_general_vs_SYNJ2cryptic <- total_mutations_each_cancer |> 
  left_join(total_mutations_each_cancer_with_SYNJ2cryptic, by=c("cancer_abbrev")) |>   
  mutate(percent_with_cryptic = (total_mutations_cryptic/total_mutations)*100) |> 
  #16% of all mutations in LGG cancer are in cases with STMN2 cryptic events
  drop_na() |> 
  arrange(-percent_with_cryptic)

#kable(mutations_each_cancer_general_vs_SYNJ2cryptic, caption = "Mutational burden of cryptic SYNJ2 cases.")

# calculating mutation burden

m_SYNJ2_analysis <- cBio_clinical |>
  left_join(SYNJ2_clinical_jir,by = c("case_submitter_id")) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  select(case_submitter_id, study_id, cancer_abbrev, mutation_count, cryptic_count, anno_count) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  unique() |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(synj2_cryptic_detected = cryptic_count >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_synj2 = sum(synj2_cryptic_detected,na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(synj2_cryptic_detected)) |> 
  filter(n_detected_synj2 > 2) |> 
  ggplot(aes(x = synj2_cryptic_detected,
             y = log10(mutation_count))) + 
  geom_boxplot(aes(fill = synj2_cryptic_detected)) + 
  labs(x = "Cancer Type",
       y = "Log10 Mutation Count") +
  facet_wrap(~study_start) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.format") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Cryptic SYNJ2 detected")) +
  theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 25),
    strip.text = element_text(size = 25, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )

print(m_SYNJ2_analysis)

# Survival comparisons ----------------------------------------------------

survival_SYNJ2_cryptic <- cBio_clinical |>
  left_join(SYNJ2_clinical_jir,by = c("case_submitter_id")) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  select(case_submitter_id, study_id, sample_id, cancer_abbrev, mutation_count, cryptic_count, anno_count, months_of_disease_specific_survival, overall_survival_months, disease_specific_survival_status) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  mutate(study_start = as.factor(study_start)) |> 
  mutate(months_of_disease_specific_survival = as.numeric(months_of_disease_specific_survival)) |> 
  filter(!is.na(months_of_disease_specific_survival) & months_of_disease_specific_survival != "NA") |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(synj2_cryptic_detected = cryptic_count >= 2) |> 
  group_by(cancer_abbrev) |> 
  mutate(n_detected_synj2 = sum(synj2_cryptic_detected,na.rm = TRUE)) |>
  ungroup() |> 
  filter(!is.na(synj2_cryptic_detected)) |> 
  filter(n_detected_synj2 > 2) |> 
  mutate(synj2_cryptic_detected = as.logical(synj2_cryptic_detected)) |> #changes FALSE and TRUE to 0 and 1 respectively
  distinct() |> 
  filter(!is.na(disease_specific_survival_status) & disease_specific_survival_status != "NA") |> 
  mutate(disease_specific_survival_status = str_extract(disease_specific_survival_status, "^[^:]+"),
         disease_specific_survival_status = as.numeric(disease_specific_survival_status))
# 1 = dead with tumor
# 0 = alive or dead tumor free

# survival KM curve for cases with SYNJ2 cryptic
n_SYNJ2_analysis <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ 1, 
        data = subset(survival_SYNJ2_cryptic, synj2_cryptic_detected==TRUE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: SYNJ2 Cryptic Events"
  ) +
  add_confidence_interval() 

print(n_SYNJ2_analysis)

# survival KM curve for cases with no cryptic
o_SYNJ2_analysis <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ 1, 
        data = subset(survival_SYNJ2_cryptic, synj2_cryptic_detected==FALSE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: No SYNJ2 Cryptic Events"
  ) +
  add_confidence_interval() 

print(o_SYNJ2_analysis)

# survival KM curve - both cryptic and non
SYNJ2_KM_fit <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ synj2_cryptic_detected, 
                        data = survival_SYNJ2_cryptic) 

p_SYNJ2_analysis <- ggsurvplot(SYNJ2_KM_fit, data = survival_SYNJ2_cryptic,
           main = "Kaplan-Meier Curve: All Cancer Types (SYNJ2)",
           legend.title = "SYNJ2 Cryptic Detected",
           legend.labs = c("FALSE", "TRUE"),
           legend.position = "bottom",
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           xlab = "Time (months)",
           font.x = c(10, "bold"),
           ylab = "Survival Probability",
           font.y = c(10, "bold"),
           fontsize = 4
) 

print(p_SYNJ2_analysis)


# density plot - months survival in cryptic vs non
q_SYNJ2_analysis <- survival_SYNJ2_cryptic |> 
  mutate(synj2_cryptic_detected = as.logical(synj2_cryptic_detected)) |> 
  ggplot(aes(x = months_of_disease_specific_survival, colour = synj2_cryptic_detected)) +
  geom_density() +
  labs(
    x = "Survival months with disease"  
  )
# same survival - cryptic vs non-cryptic (all cancers combined)

print(q_SYNJ2_analysis)

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

# Comparing Aneuploidy cancer-by-cancer

aneuploidy_SYNJ2 <- cBio_clinical |>
  left_join(SYNJ2_clinical_jir,by = c("case_submitter_id")) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  select(case_submitter_id, study_id, cancer_abbrev, aneuploidy_score, cryptic_count, anno_count) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  unique() |> 
  mutate(aneuploidy_score = as.numeric(aneuploidy_score)) |> 
  filter(!is.na(aneuploidy_score) & aneuploidy_score != "NA") |> 
  mutate(synj2_cryptic_detected = cryptic_count >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_synj2 = sum(synj2_cryptic_detected,na.rm = TRUE)) |>
  ungroup() |> 
  filter(!is.na(synj2_cryptic_detected)) |> 
  filter(n_detected_synj2 > 2) |> 
  group_by(study_start) |> 
  mutate(wilcox_result = if (length(unique(synj2_cryptic_detected)) < 2) {NA} 
         else {
           list(broom::tidy(wilcox.test(aneuploidy_score ~ synj2_cryptic_detected, exact = FALSE)))
         }) |> 
  ungroup() |> 
  unnest(wilcox_result)

r_SYNJ2_analysis <- aneuploidy_SYNJ2 |> 
  ggplot(aes(x = synj2_cryptic_detected,
             y = aneuploidy_score)) + 
  geom_boxplot(aes(fill = synj2_cryptic_detected)) + 
  labs(x = "Cancer Type",
       y = "Aneuploidy Score") +
  facet_wrap(~study_start) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.format") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Cryptic SYNJ2 detected")) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

print(r_SYNJ2_analysis)

dev.off()

# survival - looping through each cancer type -----------------------------

cancer_abbrev <- unique(survival_SYNJ2_cryptic$cancer_abbrev)

SYNJ2_survival_plots <- list()

pdf(file = "survival_SYNJ2.pdf", width = 10, height = 10)

for(abbrev in cancer_abbrev) {
  
  message(glue::glue("Processing cancer type: {abbrev}"))
  
  plot_data <- survival_SYNJ2_cryptic |> 
    filter(cancer_abbrev == abbrev)
  
  cancer <- survival_SYNJ2_cryptic |> 
    filter(cancer_abbrev == abbrev) |> 
    pull(cancer_abbrev) |> 
    unique()
  
  plot_name <- glue::glue("{cancer}")
  
  if (length(unique(plot_data$synj2_cryptic_detected)) == 2) {
    KM_fit <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ synj2_cryptic_detected, data = plot_data)
    
    SYNJ2_survival_plots[[abbrev]] <- ggsurvplot(KM_fit, data = plot_data,
                                                    legend.title = "Cryptic Detected",
                                                    legend.labs = c("FALSE", "TRUE"),
                                                    legend.position = "bottom",
                                                    pval = TRUE,
                                                    conf.int = TRUE,
                                                    risk.table = TRUE,
                                                    tables.height = 0.2,
                                                    tables.theme = theme_cleantable(),
                                                    xlab = "Time (months)",
                                                    font.x = c(12, "bold"),
                                                    ylab = "Survival Probability",
                                                    font.y = c(12, "bold"),
                                                    fontsize = 6
    ) +
      ggtitle(plot_name)
    
    print(SYNJ2_survival_plots[[abbrev]])
  } else {
    message(glue::glue("Skipping cancer type: {abbrev} due to insufficient data."))
  }
}

dev.off()

