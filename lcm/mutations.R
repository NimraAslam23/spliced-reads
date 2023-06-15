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


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snapcount", force = TRUE)

# pull back all mutation data from TCGA -----------------------------------

# List of TCGA projects
    #project_list <- GDCquery_projects()
project_list_1 <- c("TCGA-BRCA", "TCGA-THCA", "TCGA-UCEC", "TCGA-ACC")
project_list_2 <- c("TCGA-KICH", "TCGA-HNSC", "TCGA-LIHC", "TCGA-MESO")
project_list_3 <- c("TCGA-LAML", "TCGA-KIRP", "TCGA-KIRC", "TCGA-GBM")
project_list_4 <- c("TCGA-LGG", "TCGA-SARC", "TCGA-PCPG", "TCGA-READ")
project_list_5 <- c("TCGA-PAAD", "TCGA-LUAD", "TCGA-PRAD", "TCGA-OV")
project_list_6 <- c("TCGA-LUSC", "TCGA-TGCT", "TCGA-THYM", "TCGA-UVM")
project_list_7 <- c("TCGA-SKCM", "TCGA-UCS", "TCGA-STAD")

all_mutation_data <- data.frame()

for (project in project_list_1) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
    GDCdownload(query)
    maf_1 <- GDCprepare(query)
    if (!is.null(maf_1$data) && nrow(maf_1$data) > 0) {
      all_mutation_data <- rbind(all_mutation_data, maf_1$data)
    }
}


for (project in project_list_2) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query)
  maf_2 <- GDCprepare(query)
  if (!is.null(maf_2$data) && nrow(maf_2$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf_2$data)
  }
}


for (project in project_list_3) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query)
  maf_3 <- GDCprepare(query)
  if (!is.null(maf_3$data) && nrow(maf_3$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf_3$data)
  }
}


for (project in project_list_4) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query)
  maf_4 <- GDCprepare(query)
  if (!is.null(maf_4$data) && nrow(maf_4$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf_4$data)
  }
}

for (project in project_list_5) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query)
  maf_5 <- GDCprepare(query)
  if (!is.null(maf_5$data) && nrow(maf_5$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf_5$data)
  }
}

for (project in project_list_6) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query)
  maf_6 <- GDCprepare(query)
  if (!is.null(maf_6$data) && nrow(maf_6$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf_6$data)
  }
}

for (project in project_list_7) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query)
  maf_7 <- GDCprepare(query)
  if (!is.null(maf_7$data) && nrow(maf_7$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf_7$data)
  }
}


all_m




# still need to go through this code --------------------------------------

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
    GDCdownload(query) #What is this line doing? Does it actually assign to somethign? 
    maf <- GDCprepare(query)
    if (!is.null(maf$data) && nrow(maf$data) > 0) {
      write.csv(maf)
    }
    
  }else{
    print('file exists')
    print(file_output)
  }
  
}




all_mutation_data <- data.frame()

for (project in project_list) {
  query <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  
  GDCdownload(query)
  
  maf <- GDCprepare(query)
  
  if (!is.null(maf$data) && nrow(maf$data) > 0) {
    all_mutation_data <- rbind(all_mutation_data, maf$data)
  }
}







