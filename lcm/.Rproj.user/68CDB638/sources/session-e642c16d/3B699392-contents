library(tidyverse)
library(repurrrsive)
library(jsonlite)
library(clipr)

# TCGA Data and Plots -----------------------------------------------------

jir_new$gdc_cases.submitter_id
write_clip(jir_new$gdc_cases.submitter_id)

STMN2_events_mutations <- read.table("STMN2_events_mutations.tsv", sep = "\t", header = TRUE, na.strings="")

# Import data from TCGA ---------------------------------------------------

STMN2_events_genes <- read_json("STMN2_events_genes.json")

STMN2_events_genes <- tibble(STMN2_events_genes = parse_json(STMN2_events_genes))


