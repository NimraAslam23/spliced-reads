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

pdf(file = "figures_STMN2_analysis.pdf")

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



# STMN2 cryptic and annotated reads
gene_name <- "STMN2" 
snapcount_coords_cryptic <- "chr8:79611215-79616821"
snapcount_coords_annotated <- "chr8:79611215-79636801"
strand_code <- "+"

# querying for cryptic and annotated STMN2 reads, join cryptic and anno df, add jir column
# STMN2_query <- combine_two_junctions("STMN2", "chr8:79611215-79616821", "chr8:79611215-79636801", "+")
# write.table(STMN2_query, file="STMN2_query.txt", sep=",")
STMN2_query <- read.table("STMN2_query.txt", sep=",")

# write_clip(STMN2_query$gdc_cases.submitter_id)

# Import clinical data from TCGA

STMN2_clinical_orig <- read.csv("STMN2_clinical.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)

STMN2_clinical <- STMN2_clinical_orig

STMN2_clinical <- STMN2_clinical|> 
  select(case_submitter_id, project_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage) |> 
  separate(project_id, into = c("project", "cancer_type")) |> 
  select(-project)

head(STMN2_clinical)

# join query table with clinical table
STMN2_clinical_jir <- STMN2_query |> 
  rename("case_submitter_id" = "gdc_cases.submitter_id") |> 
  left_join(STMN2_clinical, by=c("case_submitter_id")) |> 
  janitor::clean_names()

STMN2_clinical_jir <- STMN2_clinical_jir |> 
  select(-c(gdc_cases_diagnoses_tumor_stage, age_at_index, gdc_cases_demographic_gender, ajcc_pathologic_stage, ethnicity, race, gender)) |>   rename("cancer_abbrev" = "cancer_type")

# filter for cryptic coverage 
STMN2_clinical_jir_cryptic <- STMN2_clinical_jir |> filter(cryptic_count > 2)

head(STMN2_clinical_jir_cryptic)


# Cancers with STMN2 expression (in general)

# bar plot - cancer types with general STMN2 expression
a_STMN2_analysis <- STMN2_clinical_jir |> 
  distinct() |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_abbrev)))) +
  geom_bar(aes(fill = cancer_abbrev)) +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = "STMN2 Expression in Different Cancer Types" 
  ) +
  theme(legend.position = "none", plot.title = element_text(size=14)) 
# BRCA = breast cancer; LGG = low-grade gliomas (brain tumours)
print(a_STMN2_analysis)

# Cancers with cryptic STMN2 events

# bar plot - cancer types with STMN2 events, primary sites with cryptic
b_STMN2_analysis <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_abbrev)))) +
  geom_bar(aes(fill = cancer_abbrev)) +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases",
    title = "Cryptic STMN2 Expression in Different Cancer Types" 
  ) +
  theme(
    #plot.title = element_text(size=15),
    # axis.title.x = element_text(size = 15, face = "bold"),
    # axis.title.y = element_text(size = 15, face = "bold"),
    # axis.text.x = element_text(size = 12),
    # axis.text.y = element_text(size = 12),
    legend.position = "none") 

print(b_STMN2_analysis)

c_STMN2_analysis <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(gdc_cases_project_primary_site)))) +
  geom_bar(aes(fill = gdc_cases_project_primary_site)) +
  coord_flip() +
  labs(
    x = "Primary Site of Cancer",
    y = "Number of Cases",
    title = "Cryptic STMN2 Expression in Different Cancer Sites" 
  ) +
  theme(
    #plot.title = element_text(size=15),
    # axis.title.x = element_text(size = 15, face = "bold"),
    # axis.title.y = element_text(size = 15, face = "bold"),
    # axis.text.x = element_text(size = 12),
    # axis.text.y = element_text(size = 12),
    legend.position = "none") 

print(c_STMN2_analysis)

## Junction coverage

# boxplot - STMN2 junction coverage
#overall STMN2 junction coverage in different cancer sites
d_STMN2_analysis <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = junction_avg_coverage, 
             y = fct_reorder(gdc_cases_project_primary_site, junction_avg_coverage, median))) +
  geom_boxplot(aes(fill = gdc_cases_project_primary_site)) +
  labs(
    x = "Junction Average Coverage",
    y = "Primary Site of Cancer",
    title = "Average STMN2 Junction Coverage in Different Cancer Sites"
  ) +
  theme(
    #plot.title = element_text(size=15),
    # axis.title.x = element_text(size = 15, face = "bold"),
    # axis.title.y = element_text(size = 15, face = "bold"),
    # axis.text.x = element_text(size = 12),
    # axis.text.y = element_text(size = 12),
    legend.position = "none") +
  geom_signif(comparisons = list(c("Stomach", "Brain"), c("Brain", "Lung")),
              map_signif_level = TRUE,
              y_position = c(75, 60))

print(d_STMN2_analysis)

#adding reads per million (rpm) column for cryptic coverage
STMN2_clinical_jir_cryptic <- STMN2_clinical_jir_cryptic |> 
  mutate(rpm = (cryptic_count/junction_coverage)*1000000)

e_STMN2_analysis <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = rpm,
             y = fct_reorder(gdc_cases_project_primary_site, rpm, median))) +
  geom_boxplot(aes(fill = gdc_cases_project_primary_site)) +
  labs(
    x = "reads per million",
    y = "Primary Site of Cancer",
    title = "STMN2 Cryptic Coverage in Different Cancer Sites"
  ) +
  theme(
    #plot.title = element_text(size=15),
    # axis.title.x = element_text(size = 15, face = "bold"),
    # axis.title.y = element_text(size = 15, face = "bold"),
    # axis.text.x = element_text(size = 12),
    # axis.text.y = element_text(size = 12),
    legend.position = "none")

print(e_STMN2_analysis)

# boxplot - cryptic STMN2 expression in different cancer sites
f_STMN2_analysis <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  ggplot(aes(x = cryptic_count, 
             y = fct_reorder(gdc_cases_project_primary_site, cryptic_count, median))) +
  geom_boxplot(aes(fill = gdc_cases_project_primary_site)) +
  labs(
    x = "STMN2 Cryptic Coverage",
    y = "Primary Site of Cancer",
  ) +
  theme(plot.title = element_text(size=10),
        legend.position = "none")

print(f_STMN2_analysis)
#On average, cancers in the stomach and breast have the greatest number of reads supporting cryptic STMN2 events. 
#Cancers of the brain have significantly fewer reads surpporting cryptic events. 

## Which cancers have the most cryptic STMN2 events?

STMN2_events_different_cancers <- STMN2_clinical_jir_cryptic |>
  drop_na(cancer_abbrev) |>
  janitor::tabyl(cancer_abbrev) |> arrange(desc(percent))

#kable(STMN2_events_different_cancers, caption = "PCPG and GBM cancers have high cryptic STMN2 expression.")

g_STMN2_analysis <- STMN2_events_different_cancers |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent), 
             y = percent,
             fill = cancer_abbrev)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = 
                  ifelse(cancer_abbrev %in% c("PCPG", "GBM", "STAD"), 
                         sprintf("%.1f%%", percent * 100), "")), 
            size = 3.2, vjust=0.5, hjust = 0) +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Percentage of Cases",
    title = "Fraction of all cryptic STMN2 cases that are of each cancer type"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format())

print(g_STMN2_analysis)

## Where are the cancers with cryptic STMN2 events located?

STMN2_events_primary_sites <- STMN2_clinical_jir_cryptic |> 
  drop_na() |> 
  janitor::tabyl(gdc_cases_project_primary_site) |> arrange(-percent)

#kable(STMN2_events_primary_sites, caption = "Cancers with cryptic STMN2 events are found primarily in the adrenal gland and brain")

h_STMN2_analysis <- STMN2_events_primary_sites |> 
  ggplot(aes(x = reorder(gdc_cases_project_primary_site, percent), 
             y = percent,
             fill = gdc_cases_project_primary_site)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = 
                  ifelse(gdc_cases_project_primary_site %in% c("Adrenal Gland", "Brain"), 
                         sprintf("%.1f%%", percent * 100), "")), 
            size = 3.1, vjust=0.5, hjust = 0) +
  coord_flip() +
  labs(
    x = "Primary Sites of Cancer",
    y = "Percentage of Cases",
    title = "Fraction of all cryptic STMN2 cases that are in each cancer site"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format())

print(h_STMN2_analysis)
# 46% of cancers with cryptic STMN2 events are in the adrenal gland

# Mutational burden

# data import - cBioPortal Clinical Data

cBio_all_clinical_orig <- read.csv("cBio_clinical_data.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)

cBio_clinical <- cBio_all_clinical_orig

cBio_clinical <- cBio_clinical |> 
  select(Study.ID, Patient.ID, Sample.ID, Aneuploidy.Score, Cancer.Type, TCGA.PanCanAtlas.Cancer.Type.Acronym, Cancer.Type.Detailed,
         Months.of.disease.specific.survival, Mutation.Count, Overall.Survival..Months., Disease.specific.Survival.status)

cBio_clinical <- cBio_clinical |> 
  rename("case_submitter_id" = "Patient.ID") |> 
  rename("cancer" = "Cancer.Type.Detailed") |> 
  rename("cancer_abbrev" = "TCGA.PanCanAtlas.Cancer.Type.Acronym")

# join cBio clinical data with STMN2 cryptic table
STMN2_cryptic_cBio <- STMN2_clinical_jir_cryptic |> 
  left_join(cBio_clinical, by = "case_submitter_id") |> 
  janitor::clean_names() |> 
  select(-c(cgc_case_primary_site, cancer_abbrev_y, cancer_type)) |> 
  rename("cancer_abbrev" = "cancer_abbrev_x")

write.table(STMN2_cryptic_cBio, file = "STMN2_cryptic_cBio.txt", sep=",")

# boxplot - n mutations for each cancer type
i_STMN2_analysis <- STMN2_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na() |>
  filter(mutation_count !="NA") |> 
  ggplot(aes(x = fct_reorder(cancer_abbrev, mutation_count, median),
             y = mutation_count)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Mutation Count",
    title = "UCEC and LUSC cancers have the greatest number of mutations"
  ) +
  geom_signif(comparisons = list(c("LUSC", "UCEC")),
              map_signif_level = TRUE,
              y_position = c(1000)
  ) 

print(i_STMN2_analysis)

#Figure 6: Among cancer patients with STMN2 cryptic expression, uterine corpus endometrial carcinoma (UCEC) has the greatest number of mutations.

# STMN2 cryptic coverage vs mutation count
j_STMN2_analysis <- STMN2_cryptic_cBio |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na() |> 
  filter(mutation_count < 400, cryptic_count > 1) |> 
  ggplot(aes(x = as.factor(cryptic_count), 
             y = mutation_count)) +
  geom_boxplot() +
  labs(
    x = "Number of STMN2 cryptic events",
    y = "Mutation Count"
  )

print(j_STMN2_analysis)

k_STMN2_analysis <- STMN2_cryptic_cBio |> 
  drop_na() |> 
  filter(mutation_count < 2500) |> 
  ggplot(aes(x = cryptic_count, y = mutation_count)) +
  geom_hex() +
  labs(
    x = "Number of STMN2 cryptic events",
    y = "Mutation Count"
  )

print(k_STMN2_analysis)

#Figure 7: There is no correlation between number of STMN2 cryptic events and number of mutations.

## Fraction of each cancer that has cryptic STMN2 events 

# barplot - fraction of cases of each cancer that have cryptic STMN2
total_each_cancer <- cBio_clinical |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev = n())

total_each_cancer_with_STMN2cryptic <- STMN2_cryptic_cBio |> 
  group_by(cancer_abbrev) |> 
  summarise(total_cancer_abbrev_with_cryptic = n())

total_each_cancer_general_vs_STMN2cryptic <- total_each_cancer |> 
  left_join(total_each_cancer_with_STMN2cryptic, by=c("cancer_abbrev")) |> 
  mutate(percent_with_cryptic = (total_cancer_abbrev_with_cryptic/total_cancer_abbrev))

l_STMN2_analysis <- total_each_cancer_general_vs_STMN2cryptic |> 
  drop_na() |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent_with_cryptic), y = percent_with_cryptic)) +
  geom_bar(aes(fill = cancer_abbrev), stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Percentage of Cases",
    title = "Percentage of Cases with Cryptic STMN2 Events"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format())

print(l_STMN2_analysis)

# fraction of mutations in each cancer type that is found in cases with cryptic
total_mutations_each_cancer <- cBio_clinical |> 
  janitor::clean_names() |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na(mutation_count) |> 
  filter(grepl("^\\d+$", mutation_count)) |>
  group_by(cancer_abbrev) |> 
  summarise(total_mutations = sum(mutation_count))

total_mutations_each_cancer_with_STMN2cryptic <- STMN2_cryptic_cBio |> 
  janitor::clean_names() |> 
  mutate_at("mutation_count", as.numeric) |> 
  drop_na(mutation_count) |> 
  filter(mutation_count < 2500) |> 
  group_by(cancer_abbrev) |> 
  summarise(total_mutations_cryptic = sum(mutation_count))

mutations_each_cancer_general_vs_STMN2cryptic <- total_mutations_each_cancer |> 
  left_join(total_mutations_each_cancer_with_STMN2cryptic, by=c("cancer_abbrev")) |>   
  mutate(percent_with_cryptic = (total_mutations_cryptic/total_mutations)*100) |> 
  #16% of all mutations in LGG cancer are in cases with STMN2 cryptic events
  drop_na() |> 
  arrange(-percent_with_cryptic)

#kable(mutations_each_cancer_general_vs_STMN2cryptic, caption = "Mutational burden of cryptic STMN2 cases.")

# calculating mutation burden for general vs cryptic
m_STMN2_analysis <- cBio_clinical |>
  left_join(STMN2_clinical_jir,by = c("case_submitter_id")) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  select(case_submitter_id, study_id, cancer_abbrev, mutation_count, cryptic_count, anno_count) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  unique() |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(stmn2_cryptic_detected = cryptic_count >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_stmn2 = sum(stmn2_cryptic_detected,na.rm = TRUE)) |>
  ungroup() |> 
  filter(!is.na(stmn2_cryptic_detected)) |> 
  filter(n_detected_stmn2 > 2) |> 
  ggplot(aes(x = stmn2_cryptic_detected,
             y = log10(mutation_count))) + 
  geom_boxplot(aes(fill = stmn2_cryptic_detected)) + 
  labs(x = "Cancer Type",
       y = "Log10 Mutation Count") +
  facet_wrap(~study_start) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), 
                     label = "p.format",
                     label.size = 6) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Cryptic STMN2 detected")) +
  theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 25),
    strip.text = element_text(size = 25, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
  )

print(m_STMN2_analysis)

# Survival Comparisons

survival_STMN2_cryptic <- cBio_clinical |> 
  left_join(STMN2_clinical_jir, by = c("case_submitter_id")) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  select(case_submitter_id, study_id, sample_id, cancer_abbrev, mutation_count, cryptic_count, anno_count, months_of_disease_specific_survival, overall_survival_months, disease_specific_survival_status) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  mutate(study_start = as.character(study_start)) |> 
  mutate(months_of_disease_specific_survival = as.numeric(months_of_disease_specific_survival)) |> 
  filter(!is.na(months_of_disease_specific_survival) & months_of_disease_specific_survival != "NA") |>
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(stmn2_cryptic_detected = cryptic_count >= 2) |> 
  group_by(study_start) |> 
  mutate(n_detected_stmn2 = sum(stmn2_cryptic_detected,na.rm = TRUE)) |>
  ungroup() |> 
  filter(!is.na(stmn2_cryptic_detected)) |> 
  filter(n_detected_stmn2 > 2) |> 
  mutate(stmn2_cryptic_detected = as.logical(stmn2_cryptic_detected)) |> #changes FALSE and TRUE to 0 and 1 respectively
  distinct() |> 
  filter(!is.na(disease_specific_survival_status) & disease_specific_survival_status != "NA") |> 
  mutate(disease_specific_survival_status = str_extract(disease_specific_survival_status, "^[^:]+"),
         disease_specific_survival_status = as.numeric(disease_specific_survival_status))
# 1 = dead with tumor
# 0 = alive or dead tumor free

# survival KM curve for cases with STMN2 cryptic
n_STMN2_analysis <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ 1, 
        data = subset(survival_STMN2_cryptic, stmn2_cryptic_detected==TRUE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: STMN2 Cryptic Events"
  ) +
  add_confidence_interval() 

print(n_STMN2_analysis)

# survival KM curve for cases with no STMN2 cryptic
o_STMN2_analysis <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ 1, 
        data = subset(survival_STMN2_cryptic, stmn2_cryptic_detected==FALSE)) |> 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall Survival Probability",
    title = "Kaplan-Meier Curve: No STMN2 Cryptic Events"
  ) +
  add_confidence_interval() 

print(o_STMN2_analysis)
# survival KM curves STMN2 - both cryptic and non

STMN2_KM_fit <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ stmn2_cryptic_detected, data = survival_STMN2_cryptic) 

p_STMN2_analysis <- ggsurvplot(STMN2_KM_fit, data = survival_STMN2_cryptic, 
           main = "Kaplan-Meier Curve: All Cancer Types (STMN2)",
           legend.title = "STMN2 Cryptic Detected",
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

print(p_STMN2_analysis)

# density plot - months survival in cryptic vs non-cryptic STMN2
q_STMN2_analysis <- survival_STMN2_cryptic |> 
  mutate(stmn2_cryptic_detected = as.logical(stmn2_cryptic_detected)) |> 
  ggplot(aes(x = months_of_disease_specific_survival, colour = stmn2_cryptic_detected)) +
  geom_density() +
  labs(
    x = "Survival months with disease"  
  )

print(q_STMN2_analysis)

# number of STMN2 detected and non in each cancer type
STMN2_n_detected_non_cancer_abbrev <- survival_STMN2_cryptic |> 
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

# Comparing aneuploidy cancer-by-cancer

# boxplot aneuploidy score cancer-by-cancer in cryptic vs none
aneuploidy_STMN2 <- cBio_clinical |>
  left_join(STMN2_clinical_jir,by = c("case_submitter_id")) |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  select(case_submitter_id, study_id, cancer_abbrev, aneuploidy_score, cryptic_count, anno_count) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |> 
  unique() |> 
  mutate(aneuploidy_score = as.numeric(aneuploidy_score)) |> 
  filter(!is.na(aneuploidy_score) & aneuploidy_score != "NA") |> 
  mutate(stmn2_cryptic_detected = cryptic_count >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_stmn2 = sum(stmn2_cryptic_detected,na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(stmn2_cryptic_detected)) |> 
  filter(n_detected_stmn2 > 2) |> 
  group_by(study_start) |> 
  mutate(wilcox_result = if (length(unique(stmn2_cryptic_detected)) < 2) {NA} 
         else {
           list(broom::tidy(wilcox.test(aneuploidy_score ~ stmn2_cryptic_detected, exact = FALSE)))
         }) |> 
  ungroup() |> 
  unnest(wilcox_result)

r_STMN2_analysis <- aneuploidy_STMN2 |> 
  ggplot(aes(x = stmn2_cryptic_detected,
             y = aneuploidy_score)) + 
  geom_boxplot(aes(fill = stmn2_cryptic_detected)) + 
  labs(x = "Cancer Type",
       y = "Aneuploidy Score") +
  facet_wrap(~study_start) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.format") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Cryptic STMN2 detected")) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

print(r_STMN2_analysis)

dev.off()

# survival - looping through each cancer type -----------------------------

cancer_abbrev <- unique(survival_STMN2_cryptic$cancer_abbrev)

STMN2_survival_plots <- list()

pdf(file = "survival_STMN2.pdf", width = 10, height = 10)

for(abbrev in cancer_abbrev) {
  
  message(glue::glue("Processing cancer type: {abbrev}"))
  
  plot_data <- survival_STMN2_cryptic |> 
    filter(cancer_abbrev == abbrev)
  
  cancer <- survival_STMN2_cryptic |> 
    filter(cancer_abbrev == abbrev) |> 
    pull(cancer_abbrev) |> 
    unique()
  
  plot_name <- glue::glue("{cancer}")
  
  if (length(unique(plot_data$stmn2_cryptic_detected)) == 2) {
    KM_fit <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ stmn2_cryptic_detected, data = plot_data)
    
    STMN2_survival_plots[[abbrev]] <- ggsurvplot(KM_fit, data = plot_data,
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
    
    print(STMN2_survival_plots[[abbrev]])
  } else {
    message(glue::glue("Skipping cancer type: {abbrev} due to insufficient data."))
  }
}

dev.off()


# number of STMN2 detected and non in each cancer type --------------------

STMN2_n_detected_non_cancer_abbrev <- survival_STMN2_cryptic |> 
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
