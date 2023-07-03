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
} 

is_null = purrr::map(tmp, function(df){is.null(dim(df))}) #anonymous function function(df)
tcga_cryptics_metatable = tmp[which(is_null == FALSE)] |> rbindlist() |> 
  mutate(coords_cryptic = two_junc[id,]$incl)

write.table(tcga_cryptics_metatable, file="tcga_cryptics_metatable.txt", sep=",")


