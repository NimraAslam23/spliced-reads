library(tidyverse)
library(repurrrsive)
library(jsonlite)
library(clipr)
library(naniar)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(tidyr)
library(knitr)
library(snapcount)
library(ggsignif)
library(TCGAbiolinks)
library(ggsurvfit)
library(survival)
library(survminer)
library(data.table)
library(dplyr)


# add STMN2, ARHGAP32, SYNJ2 to two_junc table ----------------------------

two_junc <- two_junc |> 
  add_row(excl = "chr8:79611215-79636801", strand = "+", incl = "chr8:79611215-79616821", annot.strand = "+", annot.gene_id = "STMN2") |> 
  add_row(excl = "chr11:128988126-128998318", strand = "-", incl = "chr11:128992047-128998318", annot.strand = "-", annot.gene_id = "ARHGAP32") |> 
  add_row(excl = "chr6:158017291-158028755", strand = "+", incl = "chr6:158017291-158019983", annot.strand = "+", annot.gene_id = "SYNJ2")

# query TCGA for additional cryptic events and combine into metatable --------

tmp = list(data.table(matrix()))

for (val in rownames(two_junc)){
  id = as.numeric(val)
  #print(id)
  #print(two_junc[id,]$incl) # take the row for that value of 'id' and give the inclusion results
  tmp[[id]] = combine_two_junctions(snapcount_coords_annotated = two_junc[id,]$excl,
                                    snapcount_coords_cryptic = two_junc[id,]$incl,
                                    strand_code = two_junc[id,]$annot.strand,
                                    gene_name = two_junc[id,]$annot.gene_id)
  
  tmp[[id]] = tmp[[id]] |> mutate(coords_cryptic = two_junc[id,]$incl) # inside loop
} 

is_null = purrr::map(tmp, function(df){is.null(dim(df))}) #anonymous function function(df)
tcga_cryptics_metatable = tmp[which(is_null == FALSE)] |> rbindlist() 
  

write.table(tcga_cryptics_metatable, file="tcga_cryptics_metatable.txt", sep=",")

#write_clip(tcga_cryptics_metatable$gdc_cases.submitter_id)

# join clinical data ------------------------------------------------------

tcga_cryptics_metatable <- read.table("tcga_cryptics_metatable.txt", sep=",")

tcga_cryptics_metatable <- tcga_cryptics_metatable |> 
  janitor::clean_names() |> 
  select(-c(rail_id, gdc_cases_diagnoses_tumor_stage, gdc_cases_samples_sample_type, cgc_case_primary_site)) |> 
  rename("case_submitter_id" = "gdc_cases_submitter_id") |> 
  left_join(cBio_clinical, by = "case_submitter_id") |> 
  janitor::clean_names()


# cancer types with cryptic [gene] expression -----------------------------

tcga_cryptics_metatable |> 
  filter(cryptic_count > 2) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |>
  drop_na(study_start) |> 
  ggplot(aes(x = fct_rev(fct_infreq(study_start)))) +
  geom_bar(aes(fill = study_start)) +
  coord_flip() +
  facet_wrap(~gene_name) +
  labs(
    x = "Cancer Type",
    y = "Number of Cases"
  )


# cryptic coverage (rpm)  -------------------------------------------------

tcga_cryptics_metatable |> 
  add_rpm_column() |> 
  filter(cryptic_count > 2) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |>
  drop_na(study_start) |> 
  ggplot(aes(x = rpm,
             y = fct_reorder(study_start, rpm, median))) +
  geom_boxplot(aes(fill = study_start)) +
  facet_wrap(~gene_name, nrow = 7) +
  labs(
    x = "Cryptic Coverage (reads per million)",
    y = "Primary Site of Cancer"
  )

tcga_cryptics_metatable |> 
  add_rpm_column() |> 
  filter(cryptic_count > 2) |> 
  filter(gene_name %in% c("AARS1", "ARHGAP32", "CAMK2B", "CDK7", "EPB41L4A", "KALRN", "KNDC1", 
                           "SETD5", "STMN2", "STRA6", "SYNJ2", "TMEM260", "TRAPPC12")) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |>
  drop_na(study_start) |> 
  ggplot(aes(x = rpm,
             y = study_start)) +
  geom_boxplot(aes(fill = study_start)) +
  facet_wrap(~gene_name, ncol = 3) +
  labs(
    x = "Cryptic Coverage (reads per million)",
    y = "Primary Site of Cancer"
  )

# which cancers have the most cryptic events? -----------------------------

table <- list(data.table(matrix()))

tcga_cryptics_metatable |> 
for (gene in tcga_cryptics_metatable$gene_name) {
  table[[gene]] = janitor::tabyl(cancer_abbrev) |> 
    mutate(percentage = percent*100)
}

is_null <- purr::map(table, function(df){is.null(dim(df))})
cryptic_events_different_cancers <- table[which(is_null == FALSE)] |> 
  cbindlist() 


tcga_cryptics_metatable |> 
  filter(cryptic_count > 2) |> 
  separate(study_id,into = ('study_start'),remove = FALSE) |>
  drop_na(study_start) |> 
  mutate()


pdf("my_plots.pdf")
for(i in thing){
  a_pplot = ggplot()....
  print(a_pplot)
}
dev.off()


# TMEM260 gene ----------------------------------------------------------------

TMEM260_clinical_jir <- tcga_cryptics_metatable |> 
  filter(gene_name == "TMEM260") |> 
  separate(study_id,into = ("study_start"),remove = FALSE)

# cancer types with general TMEM260 expression
TMEM260_clinical_jir |> 
  distinct() |> 
  drop_na() |> 
  ggplot(aes(x = fct_rev(fct_infreq(cancer_abbrev)))) +
  geom_bar(aes(fill = cancer_abbrev)) +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Number of Cases"
  ) +
  theme(legend.position = "none")
  
# TMEM260 junction coverage

    # overall TMEM260 junction coverage in different cancer sites
TMEM260_clinical_jir |> 
  filter(cryptic_count > 2) |> 
  drop_na() |> 
  ggplot(aes(x = junction_avg_coverage,
             y = fct_reorder(cancer_abbrev, junction_avg_coverage, median))) +
  geom_boxplot(aes(fill = cancer_abbrev)) +
  labs(
    x = "Junction Average Coverage",
    y = "Cancer Type"
  ) + 
  theme(legend.position = "none")

    # adding reads per million (rpm) column for cryptic coverage
TMEM260_clinical_jir |> add_rpm_column() |> 
  drop_na() |> 
  ggplot(aes(x = log10(rpm),
             y = fct_reorder(cancer_abbrev, rpm, median))) +
  geom_boxplot(aes(fill = cancer_abbrev)) +
  labs(
    x = "Cryptic Coverage (log10 reads per million)",
    y = "Cancer Type"
  ) +
  theme(legend.position = "none")

# cryptic TMEM260 expression

    # which cancers?
TMEM260_events_different_cancers <- TMEM260_clinical_jir |> 
  filter(cryptic_count > 2) |> 
  drop_na(cancer_abbrev) |> 
  janitor::tabyl(cancer_abbrev) |> arrange(desc(percent))

TMEM260_events_different_cancers |> 
  ggplot(aes(x = reorder(cancer_abbrev, percent),
             y = percent)) +
  geom_bar(aes(fill = cancer_abbrev), stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Percentage of Cases"
  ) +
  theme(legend.position = "none")

    # where are these cancers?
TMEM260_events_primary_sites <- TMEM260_clinical_jir |> 
  filter(cryptic_count > 2) |> 
  drop_na(gdc_cases_project_primary_site) |> 
  janitor::tabyl(gdc_cases_project_primary_site) |> arrange(desc(percent))

TMEM260_events_primary_sites |> 
  ggplot(aes(x = reorder(gdc_cases_project_primary_site, percent),
             y = percent)) +
  geom_bar(aes(fill = gdc_cases_project_primary_site), stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  labs(
    x = "Primary Site of Cancer",
    y = "Percentage of Cases"
  ) +
  theme(legend.position = "none")

# Mutational burden 

cBio_clinical |> 
  left_join(TMEM260_clinical_jir, by = "case_submitter_id") |> 
  janitor::clean_names() |> 
  rename("cancer_abbrev" = "cancer_abbrev_x") |> 
  rename("cancer" = "cancer_x") |> 
  select(study_id, case_submitter_id, aneuploidy_score, cancer, months_of_disease_specific_survival, 
         mutation_count, overall_survival_months, disease_specific_survival_status, cryptic_count, 
         jir, anno_count, gdc_cases_project_primary_site) |> 
  separate(study_id, into = ("study_start"), remove = FALSE) |> 
  unique() |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count !="NA") |> 
  mutate(tmem260_cryptic_detected = cryptic_count >= 2) |> 
  group_by(study_start) |> 
  mutate(n_total_samples = n_distinct(case_submitter_id)) |> 
  mutate(n_detected_tmem260 = sum(tmem260_cryptic_detected, na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(tmem260_cryptic_detected)) |> 
  filter(n_detected_tmem260 > 2) |> 
  ggplot(aes(x = tmem260_cryptic_detected,
             y = log10(mutation_count))) +
  geom_boxplot(aes(fill = tmem260_cryptic_detected)) +
  labs(x = "Cancer Type",
       y = "Log10 Mutation Count") +
  facet_wrap(~study_start) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Cryptic TMEM260 detected"))
