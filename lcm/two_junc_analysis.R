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

# need two_junc df 
# need cBio_clinical df

# functions ---------------------------------------------------------------

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

# query cryptic + anno and join into one df

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

# add rpm column

add_rpm_column <- function(df) {
  
  df |> 
    mutate(rpm = (cryptic_count/junction_coverage)*1000000)
  
  return(df)
}


# add STMN2, ARHGAP32, SYNJ2 to two_junc table ----------------------------

two_junc <- read.csv("two_junc.csv", sep = ",")

two_junc <- two_junc |> 
  add_row(excl = "chr8:79611215-79636801", strand = "+", incl = "chr8:79611215-79616821", annot.strand = "+", annot.gene_id = "STMN2") |> 
  add_row(excl = "chr11:128988126-128998318", strand = "-", incl = "chr11:128992047-128998318", annot.strand = "-", annot.gene_id = "ARHGAP32") |> 
  add_row(excl = "chr6:158017291-158028755", strand = "+", incl = "chr6:158017291-158019983", annot.strand = "+", annot.gene_id = "SYNJ2")

write.table(two_junc, file="two_junc.txt", sep=",", quote=FALSE)

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

pdf(file = "figures_two_junc_analysis.pdf")

# join clinical data ------------------------------------------------------

cBio_all_clinical_orig <- read.csv("cBio_clinical_data.tsv", sep = "\t", header = TRUE, na.strings = "", fill = TRUE)

cBio_clinical <- cBio_all_clinical_orig

cBio_clinical <- cBio_clinical |> 
  select(Study.ID, Patient.ID, Sample.ID, Aneuploidy.Score, Cancer.Type, TCGA.PanCanAtlas.Cancer.Type.Acronym, Cancer.Type.Detailed,
         Months.of.disease.specific.survival, Mutation.Count, Overall.Survival..Months., Disease.specific.Survival.status)

cBio_clinical <- cBio_clinical |> 
  rename("case_submitter_id" = "Patient.ID") |> 
  rename("cancer" = "Cancer.Type.Detailed") |> 
  rename("cancer_abbrev" = "TCGA.PanCanAtlas.Cancer.Type.Acronym")

tcga_cryptics_metatable <- read.table("tcga_cryptics_metatable.txt", sep=",")

# cryptic and non-cryptic sample sizes
cryptic_true_false_counts <- tcga_cryptics_metatable |> 
  janitor::clean_names() |> 
  select(-c(rail_id, gdc_cases_diagnoses_tumor_stage, gdc_cases_samples_sample_type, cgc_case_primary_site)) |> 
  rename("case_submitter_id" = "gdc_cases_submitter_id") |> 
  left_join(cBio_clinical, by = "case_submitter_id") |> 
  janitor::clean_names() |> 
  distinct() |> 
  mutate(cryptic_detected = cryptic_count >= 2) |> 
  group_by(cancer_abbrev) |> 
  mutate(n_detected_cryptic = sum(cryptic_detected, na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(cryptic_detected)) |> 
  filter(n_detected_cryptic > 2) |> 
  mutate(cryptic_detected = as.logical(cryptic_detected)) |> 
  distinct() |> 
  group_by(cancer_abbrev, cryptic_detected) |> 
  summarise(count = n()) |> 
  ungroup() |> 
  pivot_wider(
    names_from = cryptic_detected,
    values_from = count
  ) |> 
  rename("cryptic_true" = "TRUE") |> 
  rename("cryptic_false" = "FALSE")


tcga_cryptics_metatable <- tcga_cryptics_metatable |> 
  janitor::clean_names() |> 
  select(-c(rail_id, gdc_cases_diagnoses_tumor_stage, gdc_cases_samples_sample_type, cgc_case_primary_site)) |> 
  rename("case_submitter_id" = "gdc_cases_submitter_id") |> 
  left_join(cBio_clinical, by = "case_submitter_id") |> 
  janitor::clean_names() |> 
  distinct() |> 
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  mutate(cryptic_detected = cryptic_count >= 2) |> 
  group_by(cancer_abbrev) |> 
  mutate(n_detected_cryptic = sum(cryptic_detected, na.rm = TRUE)) |> 
  ungroup() |> 
  filter(!is.na(cryptic_detected)) |> 
  filter(n_detected_cryptic > 2) |> 
  mutate(cryptic_detected = as.logical(cryptic_detected)) 

# cancer sample sizes
cancer_sample_sizes <- tcga_cryptics_metatable |> 
  group_by(cancer_abbrev) |> 
  summarise(count = n(), na.rm=TRUE) |> 
  ungroup() 

tcga_cryptics_metatable <- tcga_cryptics_metatable |> 
  left_join(cancer_sample_sizes, by="cancer_abbrev") |> 
  rename("cancer_sample_size" = "count") |> 
  left_join(cryptic_true_false_counts, by="cancer_abbrev")

# cancer types with cryptic [gene] expression -----------------------------

a_two_junc_analysis <- tcga_cryptics_metatable |> 
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

print(a_two_junc_analysis)


# cryptic coverage (rpm)  -------------------------------------------------

b_two_junc_analysis <- tcga_cryptics_metatable |> 
  mutate(rpm = (cryptic_count/junction_coverage)*1000000) |> 
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

print(b_two_junc_analysis)

c_two_junc_analysis <- tcga_cryptics_metatable |> 
  mutate(rpm = (cryptic_count/junction_coverage)*1000000) |> 
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

print(c_two_junc_analysis)

dev.off()

# survival analysis -------------------------------------------------------

cryptic_events <- unique(tcga_cryptics_metatable$coords_cryptic)
survival_plots <- list()

tcga_cryptics_metatable <- tcga_cryptics_metatable |> 
  mutate(months_of_disease_specific_survival = as.numeric(months_of_disease_specific_survival)) |> 
  filter(!is.na(months_of_disease_specific_survival) & months_of_disease_specific_survival != "NA") |>
  mutate(mutation_count = as.numeric(mutation_count)) |> 
  filter(!is.na(mutation_count) & mutation_count != "NA") |> 
  filter(!is.na(cryptic_detected)) |> 
  filter(n_detected_cryptic > 2) |> 
  mutate(cryptic_detected = as.logical(cryptic_detected)) |> #changes FALSE and TRUE to 0 and 1 respectively
  distinct() |> 
  filter(!is.na(disease_specific_survival_status) & disease_specific_survival_status != "NA") |> 
  mutate(disease_specific_survival_status = str_extract(disease_specific_survival_status, "^[^:]+"),
         disease_specific_survival_status = as.numeric(disease_specific_survival_status))
# 1 = dead with tumor

pdf("survival_metatable.pdf", width = 10, height = 10)

for(event in cryptic_events) {
  
  gene_name <- tcga_cryptics_metatable |> 
    filter(coords_cryptic == event) |> 
    pull(gene_name) |> 
    unique()
  plot_name <- glue::glue("{event} - {gene_name}")
  
  message(glue::glue("Processing: {event} - {gene_name}"))
  
  plot_data <- tcga_cryptics_metatable |> # Use filtered data for the current event
    filter(coords_cryptic == event) 
  
  if (length(unique(plot_data$cryptic_detected)) == 2) {
  KM_fit <- survfit(Surv(months_of_disease_specific_survival, disease_specific_survival_status) ~ cryptic_detected, data = plot_data)

  survival_plots[[event]] <- ggsurvplot(KM_fit, data = plot_data,
             legend.title = "Cryptic Detected",
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
  ) + 
    ggtitle(plot_name)
  
  print(survival_plots[[event]])
  } else {
    message(glue::glue("Skipping event: {event} - {gene_name} due to insufficient data."))
  }
}

dev.off()

# which cancers have the most cryptic events? -----------------------------

cryptic_events <- unique(tcga_cryptics_metatable$coords_cryptic)
cryptic_plots <- list()
cryptic_tables <- tibble::tibble()

pdf("cryptic_metatable_plots.pdf", width = 10, height = 10)

for(event in cryptic_events) {
  
  gene_name <- tcga_cryptics_metatable |> 
    filter(coords_cryptic == event) |> 
    pull(gene_name) |> 
    unique()
  plot_name <- glue::glue("{event} - {gene_name}")
  
  cryptic_plots[[event]] <- tcga_cryptics_metatable |> 
    filter(coords_cryptic == event) |> 
    filter(cryptic_count >= 2) |> 
    drop_na() |> 
    ggplot(aes(x = cryptic_count, 
               y = fct_reorder(cancer_abbrev, cryptic_count, median))) +
    geom_boxplot(aes(fill = cancer_abbrev)) +
    scale_x_continuous(trans = scales::pseudo_log_trans()) +
    labs(
      x = "Cryptic Count",
      y = "Cancer Type",
    ) +
    ggtitle(plot_name) +
    theme(plot.title = element_text(size=10),
          legend.position = "none")
  
  print(cryptic_plots[[event]])
}

dev.off()


cryptic_tables <- list()

for(event in cryptic_events) {
  tmp <- tcga_cryptics_metatable |> 
    filter(coords_cryptic == event) |> 
    filter(cryptic_count >= 2) |> 
    drop_na(cancer_abbrev) |> 
    janitor::tabyl(cancer_abbrev) |> 
    arrange(desc(percent)) |> 
    mutate(event = event, gene_name = gene_name)
  
  cryptic_tables <- rbind(cryptic_tables, tmp)
}

# mutational burden -------------------------------------------------------

mutational_burden_metatable_plots <- list()

pdf("mutational_burden_metatable_plots.pdf", width = 10, height = 12)

for(event in cryptic_events) {
  
  gene_name <- tcga_cryptics_metatable |> 
    filter(coords_cryptic == event) |> 
    pull(gene_name) |> 
    unique()
  plot_name <- glue::glue("{event} - {gene_name}")
  
  plot_data <- tcga_cryptics_metatable |> 
    filter(coords_cryptic == event) |> 
    separate(study_id,into = ('study_start'),remove = FALSE) |> 
    unique() |> 
    mutate(mutation_count = as.numeric(mutation_count)) |> 
    filter(!is.na(mutation_count) & mutation_count != "NA") |> 
    filter(!is.na(cryptic_detected)) |> 
    filter(n_detected_cryptic > 2) 
  
  if (nrow(plot_data) > 0) {
    mutational_burden_metatable_plots[[event]] <- plot_data |>  
    ggplot(aes(x = cryptic_detected,
               y = log10(mutation_count))) + 
    geom_boxplot(aes(fill = cryptic_detected)) + 
    labs(x = "Cancer Type",
         y = "Log10 Mutation Count") +
    ggtitle(plot_name) +
    facet_wrap(~cancer_abbrev) +
    scale_y_continuous(trans = scales::pseudo_log_trans()) +
    stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.format", hide.ns = TRUE, size = 3) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(hjust = 1.5),  # Adjust the y-axis label spacing
          plot.title = element_text(size = 12),  # Adjust title font size
          strip.text = element_text(size = 10),
          plot.margin = unit(c(1, 2, 1, 2), "cm")) +
    guides(fill = guide_legend(title = "Cryptic detected"))
    
  print(mutational_burden_metatable_plots[[event]])
  } else {
  message(glue::glue("No data for {event} - {gene_name}."))
  }
}

 # for the printing pdf loop
dev.off()
