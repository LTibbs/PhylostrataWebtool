# Use this code to process Rdata results from phylostratr
# into json files for webtool visualization
# Please see places with "NOTE" for places that users might want or need to customize

library(tidyverse)
library(data.table)
library(jsonlite)

# load phylostrata results
load(paste0("phylostrata.", current.taxa, ".Rdata"))

# NOTE: if you want to combine any mrca levels 
# into a single phylostratum (ps) (e.g, family Poaceae and PACMAD), 
# do that here: 
# for example, for maize: 
results <- results %>%
  mutate(new.ps=ifelse(mrca=='131567',1, # cellular oganisms
                       ifelse(mrca=='2759',2, # domain eukarya
                              ifelse(mrca=='33090',3, # kingdom viridiplantae
                                     ifelse(mrca %in% c('35493','131221'),4, # phylum streptophyta + subphylum Streptophytina
                                            ifelse(mrca%in% c('3193','58023','78536','58024'),5, # embryophyta (land plants) +tracheophyta (vascular plants), euphyllophyta, and spermatophyta (seed plants)
                                                   ifelse(mrca %in% c('3398','1437183'),6, # class Magnoliopsida + mesangiospermae
                                                          ifelse(mrca %in% c('4447',"1437197","4734"),7,# monocots, petrosaviidae, commelinid
                                                                 ifelse(mrca=="38820",8, # Order Poales
                                                                        ifelse(mrca %in% c("4479","147370"),9, # family Poaceae and PACMAD
                                                                               ifelse(mrca=="147369",10, # subfamily Panicoideae
                                                                                      ifelse(mrca %in% c("Andropogoneae", "Tripsacinae"),11,
                                                                                             ifelse(mrca=="Zea",12,
                                                                                                    ifelse(mrca=="Zea_mays",13,
                                                                                                           ifelse(mrca %in% c("Zea_mays_mays", "temperate", "B73"),14, NA)))))))))))))))

stopifnot(sum(is.na(results$new.ps))==0)

# otherwise, if you want to keep the default strata, use:
# results <- results %>% mutate(new.ps=ps)

# get the new.ps for each phylostrata result
new.phylostrata <- left_join(phylostrata, 
                             results %>% 
                               select(mrca,ps,new.ps) %>%
                               distinct)


n.strata <- length(unique(new.phylostrata$new.ps)) # how many strata do you have?

# transform into percent and fill percent for the webtool
new.phylostrata <- phylostrata %>% 
  mutate(pct=(100/n.strata)*((n.strata-new.ps)+.9)) %>%
  mutate(fill_pct=ifelse(new.ps==1,100,pct))

# remove transcript id for json reduction
new.phylostrata <- new.phylostrata %>%
  mutate(id=str_split_i(qseqid, "_",i=1))

# output for downloads
fwrite(new.phylostrata %>% select(id, phylostratum=new.ps), 
       paste0("phylostrata_results.csv.gz"))

# make this into json and output
write(toJSON(setNames(as.data.frame(new.phylostrata %>% 
                                      select(id,strata=new.ps, fill_pct)),
                      c("id","strata","fill_pct")), pretty=T),
      "gene_fill.json")
# NOTE: will need to manually edit this json so that it starts:
# {"pages": [ { "id": 
# by adding the {"pages": part at the beginning and } at the end of the file

# Get the number of taxa in which it was found, out of the total
# get total taxa for each new level:
taxa.count.denominator <- results %>%
  select(staxid, new.ps) %>%
  distinct %>%
  count(new.ps) %>%
  rename(n.denominator=n)
# for each gene, how many taxa is it found in?
taxa.count.numerator <- results %>%
  group_by(qseqid, new.ps) %>%
  count(new.ps)%>%
  rename(n.numerator=n)%>%
  mutate(id_ps=paste0(qseqid, "_", new.ps))
# and make missing ones 0
taxa.count.zeros <- expand.grid(unique(results$qseqid),c(1:n.strata))
colnames(taxa.count.zeros) <- c("id", "new.ps")
taxa.count.zeros <- taxa.count.zeros %>%
  mutate(id_ps=paste0(id, "_", new.ps)) %>%
  filter(!(id_ps %in% unique(taxa.count.numerator$id_ps))) %>%
  mutate(n.numerator=0) %>%
  select(qseqid=id,new.ps,n.numerator,id_ps)
taxa.count.numerator <- rbind(taxa.count.numerator, taxa.count.zeros)
temp <- plyr::count(taxa.count.numerator$qseqid)
stopifnot(sum(temp$freq!=n.strata)==0) # check that all gene x strata levels are represented
rm(temp)

# combine numerator and denominator
taxa.count.table <- left_join(taxa.count.numerator %>% 
                                select(-id_ps), 
                              taxa.count.denominator)
taxa.count.table  <- taxa.count.table %>%
  mutate(detect.frac=paste0(n.numerator, "/", n.denominator)) %>%
  rename(id=qseqid, strata=new.ps)

# get best hits for each gene and stratum:
best.hits <- results%>%
  group_by(qseqid, new.ps) %>%
  slice_min(evalue,with_ties=F)

best.hits2 <- best.hits %>% 
  mutate(best.hit=ifelse(staxid %in% c("3311", "3702"), sseqid, # for ginkgo and Arabidopsis just keep ids as is
                         ifelse(grepl("^[0-9]", staxid), #Uniprot has all-numeric taxids
                                str_split_i(sseqid,"\\|",2), # pull hit id from UniProt taxids 
                                NA ))) %>% # NOTE: If you have additional non-UniProt formatted taxids, add additional conditions here to process them
  select(id=qseqid, row=new.ps, BestOrgID=staxid, BestE=evalue, best.hit)
stopifnot(sum(is.na(best.hits2$best.hit))==0)

# get species names from the ids
# NOTE: use the IDs from your own tree; these are provided as examples
tax.key <- rbindlist(list(fread("eukaryote_20240228_UPIDs.csv"),
                          fread("prokaryote_20240228_UPIDs.csv")))
best.hits3 <- left_join(best.hits2, 
                        tax.key %>% 
                          mutate(BestOrgID=as.character(taxonomy))%>%
                          select(BestOrgID,species_binomial))
# fix Ginkgo
best.hits <- best.hits3 %>%
  mutate(ifelse(BestOrgID=="3311", "Ginkgo_biloba",species_binomial))

rm(best.hits2, best.hits3)

# get proteomes of interest:
# these are for the example species for each stratum in the webtool
# NOTE: set your example species names here:
ex.species <- c("Escherichia_coli",
                "Caenorhabditis_elegans",
                "Chlamydomonas_reinhardtii",
                "Chara_braunii",
                "Ceratopteris_richardii",
                "Musa_acuminata",
                "Ananas_comosus",
                "Oryza_sativa",
                "Panicum_virgatum",
                "Sorghum_bicolor",
                "Arabidopsis_thaliana" 
)
proteome.info <- tax.key %>%
  filter(species_binomial %in% ex.species)

# read in proteomes for example species
# NOTE: these can be downloaded from uniprot via:
# https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgo%2Ccc_subcellular_location%2Corganism_id&format=tsv&query=%28proteome%3Aproteome.info${UPID}%29
uniprot.list <- vector("list", length=nrow(proteome.info))
for(j in 1:nrow(proteome.info)) {
  if(proteome.info$species_binomial[j]=="Arabidopsis_thaliana") {next} # will do arabidopsis separately because via TAIR
  uniprot.list[[j]] <- fread(paste0(proteome.info$UPID[j],".tsv.gz")) 
}
uniprot.dt <- rbindlist(uniprot.list)

# NOTE: the following code was specifically for 
# converting the TAIR proteome IDs of Arabidopsis thaliana to UniProt IDs
# in order to use the UniProt GO term and subcellular location info. 
# So, it is only needed if you want to do this.
# Issue: TAIR conversion includes genes that are not in the reference.
# So, want to prioritize reference proteome members 1st, 
# then reviewed proteins, then go from there as needed

# all proteins
# NOTE: this file is provided, but an updated version can be downloaded via 
# https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cxref_proteomes&format=tsv&query=%28%28taxonomy_id%3A3702%29%29
arab.full.uniprot <- fread("uniprotkb_taxonomy_id_3702_20250327.tsv.gz") 
arab.full.uniprot <- arab.full.uniprot %>%
  mutate(first.proteome=str_split_i(Proteomes, ":", i=1)) %>%
  mutate(first.chr=str_split_i(Proteomes, "[;:]", i=2))
# reviewed
# NOTE: this file is provided, but an updated version can be downloaded via:
# https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cxref_proteomes&format=tsv&query=%28%28taxonomy_id%3A3702%29+AND+%28reviewed%3Atrue%29%29
arab.reviewed.uniprot <- fread("uniprotkb_organism_id_3702_AND_reviewed_20250328.tsv.gz")
arab.full.uniprot <- arab.full.uniprot %>%
  mutate(reviewed=Entry %in% arab.reviewed.uniprot$Entry)
arab.full.uniprot %>% select(first.proteome, reviewed) %>% plyr::count()
# So, the majority of reviewed proteins are from the reference proteome
# But other reviewed proteins are in:
# "" proteome, UP000078284

# proteome prioritization:
#1: UP000006548 = reference proteome

# UP000516314 has assembled chromosomes, chloroplast, and mitochondrion but 89.1% BUSCO. Submitted Sep 2020

# UP000078284 has assembled chromosomes plus 410 unassembled, 95.2% BUSCO. PNAS 2016.

# UP000426265 has no assembled chromosomes, but BUSCO 99.7%. Submitted Nov 2019
# UP000434276 similar. Submitted Dec 2019

# So, rank as following:
#  1 = UP000006548
#  2 = reviewed == T
#  3 = UP000078284 or UP000516314
#  4 = UP000426265 or UP000434276
#  5 = proteome "" (blank)

arab.full.uniprot <- arab.full.uniprot %>%
  mutate(protein.rank=ifelse(first.proteome=="UP000006548", 1, 
                             ifelse(reviewed==T,2,
                                    ifelse(first.proteome %in% c("UP000078284","UP000516314"),3,
                                           ifelse(first.proteome %in% c("UP000426265","UP000434276"),4,
                                                  ifelse(first.proteome=="",5,6))))))

# for arabidopsis, need to convert between TAIR and uniprot IDs
# NOTE: this is downloaded from https://www.arabidopsis.org/download/list?dir=Proteins%2FId_conversions
arab.key <- fread("AGI2uniprot.txt", header=F) %>%
  rename(tair.name=V1,uniprot.trans.id=V2) %>%
  mutate(uniprot.gene.id=str_split_i(uniprot.trans.id, "-",1))
arab.rank.key <- left_join(arab.key, arab.full.uniprot, by=c("uniprot.gene.id" = "Entry")) %>%
  group_by(tair.name) %>%
  slice_min(protein.rank, with_ties = F) # pull a single UniProt ID for each TAIR ID
arab.results <- results %>%
  filter(staxid=="3702") 
arab.results <- left_join(arab.results, arab.rank.key, by=c("sseqid"="tair.name"))
plyr::count(arab.results$protein.rank)
# NOTE: these details are specifically from processing the maize B73 results
# and are provided as an example, but you will likely need to inspect your own "problem genes" yourself:
# So, most (~30k) are from the reference proteome
# others (<100 total) in ranks 3,4,5,
# and finally 1770 not found
sum(is.na(arab.results$`Gene Ontology (GO)`))
# OK looks like these 170 have had names updated on Uniprot. Fix them:

# some are due to having no transcript #, so look at only 1st part of gene name
problem.genes <- arab.results %>%
  filter(is.na(`Gene Ontology (GO)`)) %>%
  mutate(sseqid.front=str_split_i(sseqid, "\\.",1))

fix.arab.key <- left_join(arab.key, arab.full.uniprot, by=c("uniprot.gene.id" = "Entry"))  %>%
  mutate(tair.front=str_split_i(tair.name, "\\.",1))%>%
  group_by(tair.front) %>%
  slice_min(protein.rank, with_ties = F)

fix.problem.genes <- left_join(problem.genes %>%
                                 select(colnames(results),sseqid.front),
                               fix.arab.key, 
                               by=c("sseqid.front"="tair.front")) %>%
  filter(!is.na(Organism))%>%
  select(colnames(arab.results))
nrow(fix.problem.genes) # ok fixed 59

problem.genes <- problem.genes %>%
  filter(!(sseqid %in% fix.problem.genes$sseqid))
# 111 left

# some have capitalization issues
fix.problem.genes2 <- left_join(problem.genes %>%
                                  mutate(sseqid.front=toupper(sseqid.front))%>%
                                  select(colnames(results),sseqid.front),
                                fix.arab.key %>%
                                  mutate(tair.front=toupper(tair.front)),
                                by=c("sseqid.front"="tair.front")) %>%
  filter(!is.na(Organism)) %>%
  select(colnames(arab.results))
nrow(fix.problem.genes2) # fixed 109 more

# remaining issues:
problem.genes <- problem.genes %>%
  filter(!(sseqid %in% fix.problem.genes2$sseqid))
nrow(problem.genes) # only 2 left!

# and for these two (Zm00001eb412040_P004 and Zm00001eb325140_P001),
# both are qseqid AT4G35335.1,
# and from AGI2TAIR, At4g35335 = F4JNE0 OR F4JN00
# so just pick one -- in this case, F4JN00 is reviewed and other is not, so I'll use it

problem.genes <- problem.genes %>%
  mutate(uniprot.trans.id="F4JN00")

fix.problem.genes3 <- left_join(problem.genes %>%
                                  select(colnames(results),uniprot.trans.id),
                                fix.arab.key %>% 
                                  filter(uniprot.trans.id%in%problem.genes$uniprot.trans.id) %>%
                                  ungroup %>%
                                  select(-tair.name,-tair.front) %>%
                                  distinct,
                                by=c("uniprot.trans.id")) %>%
  filter(!is.na(Organism)) %>%
  select(colnames(arab.results))
nrow(fix.problem.genes3) # fixed 2 more

# no remaining issues:
problem.genes <- problem.genes %>%
  filter(!(sseqid %in% fix.problem.genes3$sseqid))
stopifnot(nrow(problem.genes)==0)

# add these into results
fix.arab <- rbindlist(list(fix.problem.genes, fix.problem.genes2, fix.problem.genes3))
arab.results <- arab.results %>%
  filter(!sseqid %in% fix.arab$sseqid)
arab.results <- rbind(fix.arab, arab.results)
stopifnot(sum(is.na(arab.results$uniprot.gene.id))==0)

#### Finished with Arabidopsis result formatting

# pull example taxa results for table:
ex.results <- results %>%
  filter(staxid %in% as.character(proteome.info$taxonomy)) %>%
  filter(staxid != "3702") %>% # exclude Arabidopsis because will re-join them later
  mutate(uniprot.id=str_split_i(sseqid,"\\|",2))
# add the uniprot info for them
uniprot.ex.results <- left_join(ex.results,
                                uniprot.dt %>%
                                  mutate(staxid=as.character(`Organism (ID)`)) %>%
                                  rename(uniprot.id=Entry))
stopifnot(sum(is.na(uniprot.ex.results$`Protein names`))==0)
stopifnot(sum(uniprot.ex.results$`Protein names`=="")==0)

# combine arabidopsis back in:
full.ex.results <- rbind(uniprot.ex.results, 
                         arab.results %>% 
                           rename(uniprot.id=uniprot.gene.id) %>% 
                           mutate(`Organism (ID)`="3702") %>%
                           select(colnames(uniprot.ex.results)))
# add organism info
proteome.info$Organism <- gsub("_"," ", proteome.info$species_binomial)
full.ex.results <- left_join(full.ex.results, 
                             proteome.info %>%
                               mutate(`Organism (ID)`=as.character(taxonomy)) %>%
                               select(-UPID))
full.ex.results$`Subcellular location [CC]` <- gsub("SUBCELLULAR LOCATION: ","", full.ex.results$`Subcellular location [CC]`)

# look for issues
thing=full.ex.results %>% 
  filter(is.na(`Entry Name`))
stopifnot(nrow(thing)==0)

# and make into json
uniprot.json <- left_join(best.hits %>%
                            select(id, row, BestOrg=species_binomial,
                                   BestE, BestHit=best.hit),
                          full.ex.results %>% 
                            select(id=qseqid, row=new.ps, ExOrg=Organism, ExE=evalue, 
                                   `Protein Name`=`Protein names`,
                                   `GO terms`=`Gene Ontology (GO)`, 
                                   `Subcellular Loc`=`Subcellular location [CC]`,
                                   ExHit=uniprot.id)
)
uniprot.json$BestOrg <- gsub("_", " ", uniprot.json$BestOrg)

# to save space, simplify the ginkgo ids
# best.hits %>% 
#   ungroup %>%
#   filter(BestOrgID=="3311") %>% 
#   mutate(best.hit.end=str_split_i(best.hit, "-", i=4), 
#          best.hit.start=str_split_i(best.hit, "-", i=1)) %>%
#   select(best.hit.end, best.hit.start) %>%
#   distinct
uniprot.json$BestHit <- gsub(pattern = "-Ginkgo_biloba",replacement = "",
                             uniprot.json$BestHit)

# combine with the detection fraction:
json.pct <- taxa.count.table %>%
  select(id, row=strata, detect.frac) 

json.full <- left_join(json.pct,uniprot.zea.json)

# add example organism to all rows and set row names/order
organism.key <- json.full %>%
  ungroup %>%
  select(ExOrg, row) %>%
  filter(!is.na(ExOrg)) %>%
  distinct
json.full <- left_join(json.full %>% select(-ExOrg),organism.key) %>%
  select(c(id,row,detect.frac,BestOrg, BestE, BestHit, ExOrg, ExE,`Protein Name`,
           `GO terms`, `Subcellular Loc`,ExHit))
json.full <- json.full %>%
  mutate(BestHit=ifelse(BestOrg=="Arabidopsis thaliana", ExHit, BestHit)) # align Arabidopsis TAIR ids to UniProt IDs

# format organism etc names
json.full$ExOrg <- gsub("_"," ",json.full$ExOrg)
json.full$BestOrg <- gsub("_"," ",json.full$BestOrg)
json.full$id <- str_split_i(json.full$id, pattern="_", i=1)
json.full$ExHit <- str_split_i(json.full$ExHit, pattern="_", i=1)

# output and save data as Rdata, json, and csv.gz
save(file=paste0("json.full.loc.Rdata"), json.full)
load(paste0("json.full.loc.Rdata"))

# output as csv.gz
json.full <- json.full %>%
  arrange(id,row)
fwrite(file="phylostratr_results.csv.gz", json.full)

# and json:
json.long <- json.full %>%
  arrange(id, row)%>%
  group_by(id) %>%
  nest() %>%
  mutate(
    rows = map(data, ~ .x %>%
                 select(row, detect.frac,ExOrg, `Protein Name`, `GO terms`, `Subcellular Loc`,
                        BestOrg, BestE, BestHit, ExE, ExHit) %>%
                 pmap(function(row, detect.frac, ExOrg, `Protein Name`, `GO terms`, `Subcellular Loc`,
                               BestOrg, BestE, BestHit, ExE, ExHit) {
                   list(
                     row = row,
                     detect.frac = detect.frac,
                     hits = list(
                       ExOrg = ExOrg,
                       `Protein Name` = `Protein Name`, 
                       `GO terms` = `GO terms`, 
                       `Subcellular Loc` = `Subcellular Loc`,
                       BestOrg=BestOrg,
                       BestE=BestE,
                       BestHit=BestHit,
                       ExE=ExE,
                       ExHit = ExHit
                     )
                   )
                 })
    ))
json.long.cleaned <- json.long %>%
  select(id, rows)

write(toJSON(json.long.cleaned, pretty=T),
      paste0("cleaned_uniprot_table_loc.json"))   
# NOTE: will need to manually edit this json so that it starts:
# {"pages": [ { "id": 
# by adding the {"pages": part at the beginning and } at the end of the file
# see uniprot_table_loc_mini.json in folder for example
