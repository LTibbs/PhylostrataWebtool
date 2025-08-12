library(tidyverse)
library(data.table)
library(phylostratr)
library(ape)

## set command line arguments ----
Args <- commandArgs(T)
arg.1 <- as.numeric(Args[1])
arg.2 <- as.character(Args[2])

print(arg.1)
print(arg.2)

setwd(paste0("phylostratr_main/phylostratr_",arg.2,"/"))

# set current proteome
maize.proteomes.dt <- fread(paste0("master_proteomes.csv")) 

current.proteome=maize.proteomes.dt$strata.name[arg.1]

# load tree
load(paste0("tree_filled_", arg.2, ".Rdata"))

setwd(paste0(current.proteome))
getwd()

# set focal organism
strata.cellular.organisms@focal_species=current.proteome

# run phylostratr and save results at each step
results <- strata_blast(strata.cellular.organisms, blast_args = list(nthreads = 8)) 
save(results, file = paste0("results.",current.proteome, ".1.Rdata"))

results <- results%>%
  strata_besthits() 
save(results, file = paste0("results.",current.proteome, ".2.Rdata"))

results <- results%>%  
  merge_besthits()
save(results, file = paste0("results.",current.proteome, ".3.Rdata"))

phylostrata <- stratify(results)
save(phylostrata, file =paste0("phylostrata.",current.proteome, ".Rdata"))