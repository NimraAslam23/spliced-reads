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
library(GenomicDataCommons)

if (!require("BiocManager", quietly = TRUE, force=TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks", force=TRUE)

# pull back all mutation data from TCGA -----------------------------------

    #project_list <- GDCquery_projects()
project_list <- c("TCGA-BRCA", "TCGA-THCA", "TCGA-UCEC", "TCGA-ACC", 
                  "TCGA-KICH", "TCGA-HNSC", "TCGA-LIHC", "TCGA-MESO", 
                  "TCGA-LAML", "TCGA-KIRP", "TCGA-KIRC", "TCGA-GBM", 
                  "TCGA-LGG", "TCGA-SARC", "TCGA-PCPG", "TCGA-READ", 
                  "TCGA-PAAD", "TCGA-LUAD", "TCGA-PRAD", "TCGA-OV", 
                  "TCGA-LUSC", "TCGA-TGCT", "TCGA-THYM", "TCGA-UVM", 
                  "TCGA-SKCM", "TCGA-UCS", "TCGA-STAD")

for (project in project_list) {
  file_output = glue::glue("data/{project}.csv")
  if(!file.exists(file_output)){
    query <- GDCquery(
      project = project, 
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      access = "open"
    )
    GDCdownload(query)  
    maf <- GDCprepare(query)
    if (!is.null(maf$data) && nrow(maf$data) > 0) {
      write.csv(maf)
    }
    
  }else{
    print('file exists')
    print(file_output)
  }

}

write.table(maf, file='\\Users\\nimraaslam\\Documents\\GitHub\\spliced-reads\\lcm\\all_mutation_data.txt')
write.table(maf, file='all_mutation_data.txt')


all_mutation_data_orig <- data.frame(maf)
all_mutation_data <- data.frame(maf)

# left join clinical df to mutation df ------------------------------------

all_mutation_data <- all_mutation_data |> 
  mutate(case_submitter_id = str_extract(Tumor_Sample_Barcode, "([^-]+-[^-]+-[^-]+)")) |> 
  select(Hugo_Symbol, Entrez_Gene_Id, Chromosome, Start_Position, End_Position, Strand, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele1,
         Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode, Mutation_Status, HGVSc, HGVSp, HGVSp_Short, Transcript_ID, Exon_Number, all_effects,
         Allele, Gene, One_Consequence, Consequence, Amino_acids, SYMBOL, BIOTYPE, SIFT, PolyPhen, IMPACT, VARIANT_CLASS, case_submitter_id)

mutation_clinical_data <- all_mutation_data |> 
  left_join(cBio_clinical, by="case_submitter_id") |> 
  janitor::clean_names()



