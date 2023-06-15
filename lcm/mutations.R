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

all_mutation_data <- data.frame(maf)








