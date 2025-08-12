# Code used to set up the phylostratr analysis of maize
# The results of these analyses form the basis for 
# the Maize Phylostrata Tool at https://phylostrata.maizegdb.org/

# load required packages
library(tidyverse)
library(data.table)
library(ape)
library(xml2)

library(devtools)
install_github('ltibbs/phylostratr')
# install_github('arendsee/phylostratr', force=T)
library(phylostratr)


# Download proteome quality and taxonomy information from UniProt ---------------------------------

# First time only: download info from UniProt and process
# Proteome quality: downloaded xml on Dec 23, 2024 from:
# https://www.uniprot.org/proteomes?query=%28taxonomy_id%3A131567%29+AND+%28proteome_type%3A1%29

# read in xml downloaded from uniprot
prot.xml <- xml2::read_xml("proteomes_taxonomy_id_131567_AND_proteo_2024_12_23.xml")
char.prot <- as.character(prot.xml)
split.again <- str_split(char.prot, pattern="\n")

info.list <- vector(mode="list", 
                    length=12427 # NOTE: I manually set length to current # of reference proteomes on UniProt
) 
j <- 1 # use to iterate through the reference proteomes
for(i in seq_along(split.again[[1]])) {
  current.line <- split.again[[1]][i]
  # print(current.line)
  if(str_detect(current.line, "<redundantProteome")) {next} # this is the most common line, skip it
  if(str_detect(current.line, "upid")) { # pull UPID
    current.upid <- str_match(current.line, "<upid>\\s*(.*?)\\s*</upid>")[,2]
    
    info.list[[j]] <- tibble(UPID=current.upid, taxonomy=NA, isRef=NA,
                             isRep=NA, annScore=NA, buscoScore=NA, cpd=NA,)
  }
  if(str_detect(current.line, "taxonomy")) { # pull taxonomy id
    current.tax <- str_match(current.line, "<taxonomy>\\s*(.*?)\\s*</taxonomy>")[,2]
    
    info.list[[j]]$taxonomy <- current.tax
  }
  if(str_detect(current.line, "isReferenceProteome")) { # pull if it is a reference proteome
    current.isRef <- str_match(current.line, "<isReferenceProteome>\\s*(.*?)\\s*</isReferenceProteome>")[,2]
    
    info.list[[j]]$isRef <- current.isRef
  }
  if(str_detect(current.line, "isRepresentativeProteome")) { # pull if it is a representative proteome
    current.isRep <- str_match(current.line, "<isRepresentativeProteome>\\s*(.*?)\\s*</isRepresentativeProteome>")[,2]
    
    info.list[[j]]$isRep <- current.isRep
  }
  if(str_detect(current.line, "annotationScore")) { # pull annotation score
    current.annScore <- parse_number(current.line)
    
    info.list[[j]]$annScore <- current.annScore
  }
  
  if(str_detect(current.line, "<scores name=\"busco\">")) { # indicator variable to show if currently reading the busco section of the xml
    busco.ready <- T
    # print(current.line)
    # break
  }
  
  if(str_detect(current.line,"property name=\"completed\"")) { # pull completed # of proteins
    if(busco.ready == T) {
      current.completed <- parse_number(current.line)
      
      info.list[[j]]$buscoCompleted <- current.completed
    }
  }
  
  if(str_detect(current.line,"property name=\"completedSingle\"")) { # pull completed # of proteins
    if(busco.ready == T) {
      current.completedSingle <- parse_number(current.line)
      
      info.list[[j]]$buscoCompletedSingle <- current.completedSingle
    }
  }
  
  if(str_detect(current.line,"property name=\"completedDuplicated\"")) { # pull completed # of proteins
    if(busco.ready == T) {
      current.completedDuplicated <- parse_number(current.line)
      
      info.list[[j]]$buscoCompletedDuplicated <- current.completedDuplicated
    }
  }
  
  if(str_detect(current.line,"property name=\"total\"")) { # pull completed # of proteins
    if(busco.ready == T) {
      current.total <- parse_number(current.line)
      
      info.list[[j]]$buscoTotal <- current.total
    }
  }
  
  if(str_detect(current.line,"property name=\"score\"")) { # pull busco score
    if(busco.ready == T) {
      current.buscoScore <- parse_number(current.line)
      
      info.list[[j]]$buscoScore <- current.buscoScore
      busco.ready <- F
    }
  }
  
  if(str_detect(current.line, "<scores name=\"cpd\">")) { # indicator variable to show if currently reading the cpd section of the xml
    cpd.ready <- T
    # print(current.line)
  }
  if(str_detect(current.line,"property name=\"status\"")) { # pull cpd status
    if(cpd.ready == T) {
      current.cpd <- str_match(current.line, "<property name=\"status\" value=\"\\s*(.*?)\\s*\"/>")[,2]
      
      info.list[[j]]$cpd <- current.cpd
      cpd.ready <- F
    }
  }
  
  if(str_detect(current.line, "</proteome")) { # detect the end of this proteome and move on
    j <- j+1
    cpd.ready <- F
    busco.ready <- F
    rm(current.line)
  }
  
}
# save results
info.dt <- rbindlist(info.list, fill=T)
fwrite(file = "proteome_quality.csv", info.dt)

# read results:
info.dt <- fread("proteome_quality.csv")

# combine with taxonomy info
# to get a sense of typical scores for a given group
tax.info <- c(uniprot_downstream_ids("2"), # bacteria
              uniprot_downstream_ids("2157"), # archaea
              uniprot_downstream_ids("2759") # eukarya
) %>%
  taxizedb::classification() %>%
  Filter(f=is.data.frame) 

tax.wide.list <- vector("list", length=length(tax.info))
for(i in seq_along(tax.info)) {
  tax.wide.list[[i]]<- pivot_wider(tax.info[[i]] %>% select(-id), 
                                   values_from=name, names_from=rank) %>% # pivot wider
    mutate(Taxid=names(tax.info)[[i]]) # add taxon id
  
  # find out what taxon level is missing (if any) and add it as NA
  diffs <- setdiff(c("Taxid", "superkingdom", "kingdom", "phylum", "class",
                     "order", "family", "genus", "species", "subspecies", "strain"),
                   colnames(tax.wide.list[[i]]))
  if(length(diffs)>0) {
    dif.tibble <- matrix(nrow=1, ncol=length(diffs))
    colnames(dif.tibble) <- diffs
    
    tax.wide.list[[i]] <- cbind(tax.wide.list[[i]], dif.tibble)
  }
  
  # order columns as desired
  tax.wide.list[[i]] <- tax.wide.list[[i]]  %>%
    select(Taxid, superkingdom, kingdom, phylum, class, order, family, genus, species, subspecies, strain)
}
tax.wide <- as.data.table(rbindlist(tax.wide.list))
tax.wide$Taxid <- as.numeric(tax.wide$Taxid)
fwrite(file="taxonomy.info.csv", x=tax.wide)

# manually fix some issues: 
# UP000050996 is 1637975
info.dt$taxonomy[which(info.dt$UPID=="UP000050996")] <- 1637975

# read in and combine taxonomy and quality data
tax.wide <- fread("taxonomy.info.csv")
prot.tax.info <- left_join(info.dt, tax.wide, by=c("taxonomy"="Taxid")) %>%
  filter(!is.na(taxonomy)) # manually checked these: 2 are unclassified bacteria, one I fixed above (1637975)

# fix ones with taxonomy missing
sum(is.na(prot.tax.info$superkingdom)) #52 issues
na.tax.info <- prot.tax.info %>%
  filter(is.na(superkingdom))
for(i in seq_along(na.tax.info$taxonomy)) {
  current.taxa <- na.tax.info$taxonomy[i]
  
  current.class <- taxizedb::classification(current.taxa)[[1]]
  
  if(sum(is.na(current.class))>0) {
    print(paste0("Issue with taxid ",current.taxa)) 
    next
  } # some taxa not found
  
  current.taxonomy <- pivot_wider(current.class%>% select(-id), 
                                  values_from=name, names_from=rank) %>% # pivot wider
    mutate(Taxid=current.taxa) # add taxon id
  
  
  # find out what taxon level is missing (if any) and add it as NA
  diffs <- setdiff(c("Taxid", "superkingdom", "kingdom", "phylum", "class",
                     "order", "family", "genus", "species", "subspecies", "strain"),
                   colnames(current.taxonomy))
  if(length(diffs)>0) {
    dif.tibble <- matrix(nrow=1, ncol=length(diffs))
    colnames(dif.tibble) <- diffs
    
    current.taxonomy <- cbind(current.taxonomy, dif.tibble)
  }
  
  # order columns as desired
  current.taxonomy <- current.taxonomy  %>%
    select(Taxid, superkingdom, kingdom, phylum, class,
           order, family, genus, species, subspecies, strain)
  
  # save
  na.tax.info[i,] <- left_join(na.tax.info %>% 
                                 select(-colnames(current.taxonomy)[-1]) %>%
                                 filter(taxonomy==current.taxa),
                               current.taxonomy[1,],
                               by=c("taxonomy"="Taxid"))
  rm(current.taxonomy, current.class, current.taxa)
}
sum(is.na(na.tax.info$superkingdom)) # only 5 issues left, ok

prot.tax.nona <- prot.tax.info %>% filter(!is.na(superkingdom))

prot.tax.info <- rbind(prot.tax.nona,
                       na.tax.info) 
stopifnot(sum(is.na(prot.tax.info$superkingdom))==5) # 5 known issues

# Write out combined proteome quality and taxonomy info from UniProt
fwrite("proteome_quality_taxonomy.info.csv", 
       x=prot.tax.info%>%
         filter(!is.na(superkingdom)))

# Calculate custom weights -------------------------------------------------------
# Calculating custom weights from UniProt information
# Include information from BUSCO and UniProt Annotation Score
# Manually increase weights for species of interest (model and/or crop species)
# Completely exclude proteomes with unusually low quality

# cutoff to exclude low quality:
pct.cutoff <- 10 # NOTE: percentile cutoff to exclude taxon completely, provide as # between 0 and 100

# read in the quality info from above
quality.info <- read.csv("proteome_quality_taxonomy.info.csv", na.strings = "") 
# quality.info <- quality.info %>%
#   mutate(pct.buscoCompleted=buscoCompleted/buscoTotal,
#          pct.buscoCompletedSingle=buscoCompletedSingle/buscoTotal,
#          pct.buscoCompletedDuplicated=buscoCompletedDuplicated/buscoTotal)

# calculate scores at which taxa will be completely excluded, 
# even if they are the only ones available for a group
exclude.busco <- quantile(quality.info$buscoScore, na.rm=T, probs=pct.cutoff/100)
exclude.ann <- quantile(quality.info$annScore, na.rm=T, probs=pct.cutoff/100)

taxid.0 <- unique(c(quality.info$taxonomy[quality.info$buscoScore<exclude.busco],
                    quality.info$taxonomy[quality.info$annScore<exclude.ann],
                    quality.info$taxonomy[is.na(quality.info$buscoScore)]))

min.weight <- exclude.busco-((5-exclude.ann)*5) # get the minimum of the continuous weights

# calculate and scale the weights
quality.info <- quality.info %>%
  mutate(score=ifelse(taxonomy %in% taxid.0, 0, # set zero weight to exclude it
                      1.1+(buscoScore-((5-annScore)*5)-min.weight)/((100-min.weight)/.9) )) # scale results to 1.1-2 range: subtract minimum weight and divide by (100-min.weight)/.9, then add 1.1, in order to scale results to the 1.1-2 range 

# manually increase weights by 50% for species with genome databases,
# to link to in the phylostratr output visualization (so e.g. fewer bacteria included)
# based on: https://guides.lib.berkeley.edu/bioinformatics/species-specific,
# https://www.alliancegenome.org/,
# https://www.agbiodata.org/databases

database.vec <- c("Escherichia coli", # EcoCyc
                  "Bacillus subtilis", # subtiwiki
                  "Mycoplasma genitalium", # 243273, mycoplasma/mycoplasmoides genitalium
                  "Zea mays",
                  "Arabidopsis thaliana",
                  "Drosophila melanogaster",
                  "Caenorhabditis elegans",
                  "Saccharomyces cerevisiae",
                  "Glycine max",
                  "Ananas comosus",
                  "Solanum lycopersicum",
                  "Medicago truncatula",
                  "Setaria italica",
                  "Gossypium hirsutum",
                  "Brassica rapa", 
                  "Brachypodium distachyon",
                  "Hordeum vulgare",
                  'Miscanthus giganteus',
                  # "Oryza sativa", # given as UPID instead to specify japonica
                  "Triticum aestivum"
)
up.database <- "UP000059680" # want japonica rice specifically

quality.info <- quality.info %>%
  mutate(score=ifelse(species %in% database.vec, score*1.5, score)) %>% # increase scores for species with databases
  mutate(score=ifelse(UPID %in% up.database, score*1.5, score))  %>%
  mutate(sum.tax.na=is.na(kingdom)+is.na(phylum)+ is.na(class)+ is.na(order)+ is.na(family))  %>%
  mutate(score=ifelse(sum.tax.na>=3, 0, score)) # exclude those that are not well classified (have lots of missing taxonomy information)

zero.weights <- c( # setting these to 0 (to exclude entirely) because proteomes not actually available at UniProt FTP as of Jan 31 2025
  "Candidatus Entotheonella factor",           "Candidatus Methanohalarchaeum thermophilum",
  "Candidatus Termititenax aidoneus",          "Candidatus Mcinerneyibacterium aminivorans",
  "Candidatus Nephthysia bennettiae",          "Candidatus Naiadarchaeum limnaeum",
  "Candidatus Termititenax persephonae"
) 
zero.wt.class <- c("Candidatus Entotheonellia",
                   "Candidatus Termititenacia Candida",
                   "Candidatus Termititenaci")

quality.info$score[which(quality.info$species %in% zero.weights)] <- 0
quality.info$score[which(quality.info$class %in% zero.wt.class)] <- 0

my.weights <- quality.info$score
names(my.weights) <- quality.info$taxonomy


# Make tree using custom weights ------------------------------------------

# Set focal species Zea mays
focal_taxid <- '4577' 

# build prokaryotic tree
pro.tree <- uniprot_sample_prokaryotes(downto="class", weights=my.weights,
                                       drop.names=names(my.weights)[my.weights==0])
save(pro.tree, file = "prokaryote_tree.Rdata")

load("prokaryote_tree.Rdata")

# Make focused tree around taxon of interest (Zea mays)
strata.cellular.organisms <- uniprot_strata(focal_taxid,  # get stratified relatives in UniProt
                                            from=1,  # stratum to begin from, where 1 is 'cellular organisms' (2 is default)
                                            drop.names=names(my.weights)[my.weights==0] # want to completely exclude some taxa
) %>% 
  strata_apply(f=diverse_subtree, 
               n=5, 
               weights=my.weights) %>%
  use_custom_prokaryote_tree(prokaryote.phylo=pro.tree)

# Add additional desired taxa: 
# they may be related to some already included, but are still of interest
# due to good annotations or databases we want to link to
strata.cellular.organisms <- strata.cellular.organisms %>%
  add_taxa(c('4615', # pineapple
             '3635', # cotton
             '6239',# c elegans
             '9606',# human
             '10090', # mouse
             '7955', # zebrafish
             '3311' # ginkgo
  ))

save(strata.cellular.organisms, file = paste0("zea_tree", Sys.Date(),".Rdata"))

load("zea_tree2025-02-03.Rdata")

# look at quality of chosen proteomes
chosen.quality<- quality.info %>%
  mutate(chosen=taxonomy %in% strata.cellular.organisms@tree$tip.label) %>%
  filter(chosen==T)
plyr::count(chosen.quality$superkingdom)

# Download UniProt proteomes ------------------------------------------------------

# This code is run on command line (not R) to download proteomes from UniProt:
mkdir uniprot-seqs
while IFS=',' read field1 field2 field3 field4; do
  echo "Field 1: $field1, Field 2: $field2, Field 3: $field3, Field 4: $field4, Field 5: $field5"

  field4="${field4%"${field4##*[![:space:]]}"}"   # remove trailing whitespace
  
  # download proteome files
  url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$field4/$field1/${field1}_${field2}.fasta.gz"
  
  wget -c "$url"
  
  # unzip and move files
  gunzip ${field1}_${field2}.fasta.gz
  mv ${field1}_${field2}.fasta  ./uniprot-seqs/$field2.faa
  
done < proteome_download_info.csv 
# NOTE: in this file "proteome_download_info.csv", 
# each row is one proteome to download. 
# So, each taxon in the tree for phylostratr will be one row.
# There are 4 comma-separated columns for each proteome: 
# UPID, taxonomy #, species name (as binomial with underscore instead of space), 
# and superkingdom.
  
# Now back in R again: 
# storing the filenames of downloaded proteomes
strata.cellular.organisms <- strata.cellular.organisms %>%
  uniprot_fill_strata
save(strata.cellular.organisms, file = paste0("zea_tree_filled", Sys.Date(),".Rdata"))
load("zea_tree_filled2025-02-14.Rdata")

# for Arabidopsis, ginkgo, and maize, I want to use proteomes from other sources rather than default UniProt,
# so download those:

# Arabidopsis: downloaded from https://datacommons.cyverse.org/browse/iplant/home/araport/public_data/Araport11_Release_201606/annotation/Araport11_genes.201606.pep.fasta.gz and unpacked
# and name as 3702.faa

# Ginkgo: download SGTW.faa.tar.bz2 from https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies/SGTW-Ginkgo_biloba
# then combine using cat SGTW.*FAA > 3311.faa

# Maize: will download proteomes from MaizeGDB next. 
# For now, remove the UniProt maize fasta so it is not accidentally used
unlink("uniprot-seqs/4577.faa")

# Add MaizeGDB proteomes --------------------------------------------------

# This R code will print out the code that needs to be submitted on the command line
# in order to download the MaizeGDB proteomes.
# NOTE: The file master_proteomes.csv contains 1 row per proteome
# and 6 columns: assembly, identifier, strata.name (a short name to be used for this taxon),
# haploid (T or F is this proteome haploid), annotation (NAM, CAAS, or PanAnd),
# and group (temperate, tropical, popcorn, sweetcorn, or mixed):
my.proteomes <- fread("master_proteomes.csv")
for (i in 1:nrow(my.proteomes)) {
  current.tax <- my.proteomes$assembly[i]
  current.id <- my.proteomes$identifier[i]
  
  # get the proteins, download into folder MaizeGDB_proteomes 
  if(my.proteomes$haploid[i]){
    print(paste0("wget https://download.maizegdb.org/", current.tax, "/",
                 current.tax, "a/",
                 current.tax, "a_", current.id, ".protein.fa.gz"))
    print(paste0("wget https://download.maizegdb.org/", current.tax, "/",
                 current.tax, "b/",
                 current.tax, "b_", current.id, ".protein.fa.gz"))
    
    # and index:
    print(paste0("samtools faidx ",current.tax, "a_", current.id, ".protein.fa.gz"))
    print(paste0("samtools faidx ",current.tax, "b_", current.id, ".protein.fa.gz"))
    
  } else {
    print(paste0("wget https://download.maizegdb.org/", current.tax, "/",
                 current.tax, "_", current.id, ".protein.fa.gz"))
    print(paste0("samtools faidx ",current.tax, "_", current.id, ".protein.fa.gz"))
    
  }
  
}
# fix Huangzaosi in command line, it has some proteins with unexpected characters
module load seqkit
seqkit grep --invert-match -f rm.seqs.txt Zm-Huangzaosi-REFERENCE-CAAS_FIL-1.0_Zm00095aa.1.protein.fa.gz -o temp.faa
mv temp.faa Zm-Huangzaosi-REFERENCE-CAAS_FIL-1.0_Zm00095aa.1.protein.fa.gz

# pull longest isoforms: in MaizeGDB_proteomes folder
file.list <- list.files(pattern="*.protein.fa.gz$") 
for(i in seq_along(file.list)){
  temp <- fread(paste0(file.list[i], ".fai")) %>%
    mutate(geneid=str_split_i(string = V1, pattern="_",i=1)) 
  geneids=sort(unique(temp$geneid))
  temp <- temp%>%
    group_by(geneid) %>%
    slice_max(order_by=V2, n=1, with_ties = F) %>% # V2 contains the protein length in AAs
    ungroup()
  new.geneids=sort(unique(temp$geneid))
  stopifnot(all.equal(geneids, new.geneids))
  fwrite(temp %>% select(V1), 
         paste0(gsub(file.list[i],pattern=".protein.fa.gz", replacement=""),
                ".long.seqs.txt"), col.names = F)
}
# This prints out code to submit on the command line:
for(i in seq_along(file.list)){
  current.name <- gsub(file.list[i],pattern=".protein.fa.gz", replacement="")
  print(paste0("seqkit grep -f ",paste0(current.name,
                                        ".long.seqs.txt "),
               current.name, ".protein.fa.gz -o ",
               current.name, "_longest.protein.faa"))
}
print("cp *longest.protein.faa ../uniprot-seqs") # move the contents of MaizeGDB_proteomes into the uniprot-seqs directory for analysis

# Add Zea and Andropogoneae to the tree:

load("zea_tree_filled2025-02-14.Rdata")
# add species:
strata.cellular.organisms <- add_taxa(strata.cellular.organisms, "1293079") # Zea nicaraguensis
strata.cellular.organisms <- add_taxa(strata.cellular.organisms, "4576") #  Zea diploperennis
strata.cellular.organisms <- add_taxa(strata.cellular.organisms, "4563") #  Tripsacum dactyloides
strata.cellular.organisms <- add_taxa(strata.cellular.organisms, "1079125") #  Anatherum virginicum

# Make andropogoneae tree
and.tree <- Strata(
  tree=ape::read.tree(file="Andropogoneae_tree.txt"), # NOTE: Andropogoneae_tree.txt contains the text tree of the Andropogoneae to be added to the existing phylostratr tree
  data = list(faa=list(
    Anatherum_virginicum = 'uniprot-seqs/Anatherum_virginicum.faa',
    `4558`='uniprot-seqs/4558.faa', # sorghum
    `422564`="uniprot-seqs/422564.faa", # miscanthus
    Tripsacum_dactyloides_FL = 'uniprot-seqs/Tripsacum_dactyloides_FL.faa',
    Tripsacum_dactyloides_KS = 'uniprot-seqs/Tripsacum_dactyloides_KS.faa',
    Zea_nicaraguensis = 'uniprot-seqs/Zea_nicaraguensis.faa',
    Zea_diploperennis_Momo = 'uniprot-seqs/Zea_diploperennis_Momo.faa',
    Zea_diploperennis_Gigi = 'uniprot-seqs/Zea_diploperennis_Gigi.faa',
    Zea_mays_parviglumis_TIL01 = 'uniprot-seqs/Zea_mays_parviglumis_TIL01.faa',
    Zea_mays_parviglumis_TIL11 = 'uniprot-seqs/Zea_mays_parviglumis_TIL11.faa',
    Zea_mays_mexicana_TIL18 = 'uniprot-seqs/Zea_mays_mexicana_TIL18.faa',
    Zea_mays_mexicana_TIL25 = 'uniprot-seqs/Zea_mays_mexicana_TIL25.faa',
    Zea_mays_huehuetenagensis = 'uniprot-seqs/Zea_mays_huehuetenagensis.faa',
    Tx303 = 'uniprot-seqs/Tx303.faa',
    Mo18W = 'uniprot-seqs/Mo18W.faa',
    M37W = 'uniprot-seqs/M37W.faa',
    B73 = 'uniprot-seqs/B73.faa',
    Oh7B = 'uniprot-seqs/Oh7B.faa',
    Oh43 = 'uniprot-seqs/Oh43.faa',
    MS71 = 'uniprot-seqs/MS71.faa',
    M162W = 'uniprot-seqs/M162W.faa',
    Ky21 = 'uniprot-seqs/Ky21.faa',
    B97 = 'uniprot-seqs/B97.faa',
    A632 = 'uniprot-seqs/A632.faa',
    Chang7_2 = 'uniprot-seqs/Chang7_2.faa',
    Dan340 = 'uniprot-seqs/Dan340.faa',
    Huangzaosi = 'uniprot-seqs/Huangzaosi.faa',
    Jing724 = 'uniprot-seqs/Jing724.faa',
    Jing92 = 'uniprot-seqs/Jing92.faa',
    PH207 = 'uniprot-seqs/PH207.faa',
    Xu178 = 'uniprot-seqs/Xu178.faa',
    Ye478 = 'uniprot-seqs/Ye478.faa',
    Zheng58 = 'uniprot-seqs/Zheng58.faa',
    P39 = 'uniprot-seqs/P39.faa',
    Il14H = 'uniprot-seqs/Il14H.faa',
    Tzi8 = 'uniprot-seqs/Tzi8.faa',
    NC358 = 'uniprot-seqs/NC358.faa',
    NC350 = 'uniprot-seqs/NC350.faa',
    Ki3 = 'uniprot-seqs/Ki3.faa',
    Ki11 = 'uniprot-seqs/Ki11.faa',
    CML69 = 'uniprot-seqs/CML69.faa',
    CML52 = 'uniprot-seqs/CML52.faa',
    CML333 = 'uniprot-seqs/CML333.faa',
    CML322 = 'uniprot-seqs/CML322.faa',
    CML277 = 'uniprot-seqs/CML277.faa',
    CML247 = 'uniprot-seqs/CML247.faa',
    CML228 = 'uniprot-seqs/CML228.faa',
    CML103 = 'uniprot-seqs/CML103.faa',
    S37 = 'uniprot-seqs/S37.faa',
    Hp301 = 'uniprot-seqs/Hp301.faa'
  )),
  focal_species='4577') 

strata.cellular.organisms <- replace_branch(strata.cellular.organisms, y=and.tree, node='147429') 
save(file="tree_filled_maize.Rdata", strata.cellular.organisms)
load("tree_filled_maize.Rdata")

