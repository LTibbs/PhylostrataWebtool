current_dataset=maize # provide a name for the current dataset

cd phylostratr_main/uniprot-seqs # move to uniprot-seqs directory

# 1. make diamond databases 
# need diamond installed and loaded: see https://github.com/bbuchfink/diamond_docs/blob/master/2%20Installation.MD
# here I have it installed in a conda environment 
module load miniconda
source activate conda_envs/diamond/
# make diamond databases: from within the uniprot_seqs directory:
for f in *.faa; do
    conda_envs/diamond/diamond makedb --in ${f} -d phylostratr_$current_dataset/diamond_db/${f%.faa} 
done

# copy uniprot sequences into the working directory for current dataset
cp -rs *faa ../phylostratr_${current_dataset}/uniprot-seqs

# move into directory for current dataset
cd ../phylostratr_$current_dataset

# loop through query proteomes:
# this loop is based on the master_proteomes.csv file which contains all potential query proteomes
# available from MaizeGDB (NAM, CAAS, and PanAnd projects).
# The master_proteomes.csv file has 6 columns: 
# assembly: assembly name from MaizeGDB
# identifier: identifier from MaizeGDB
# strata.name: a short name to be used for phylostratr output
# haploid: whether this accession is provided as two haploid genomes (T) or not (F)
# annotation: which annotation project this proteome comes from (NAM, CAAS, PanAnd)
# group: for NAM only, whether the accession is tropical, temperate, mixed, sweet corn, or popcorn
# If you would like to test with the small example data, replace master_proteomes.csv with test_proteome.csv
while IFS=',' read field1 field2 field3 field4 field5 field6; do
    echo "Field 1: $field1, Field 2: $field2, Field 3: $field3" 

    current_short_name=$field3

    # 2. Make file structures and copy in files for each genome
    mkdir $current_short_name
    cp -rs phylostratr_${current_dataset}/diamond_db phylostratr_${current_dataset}/${current_short_name}/
    cp -rs phylostratr_${current_dataset}/uniprot-seqs phylostratr_${current_dataset}/${current_short_name}/

    # 3. run diamond blastp with current proteome as query
    cd phylostratr_${current_dataset}/${current_short_name}/uniprot-seqs/
    sbatch diamond_blast.sh ${current_short_name}

    # return to main directory for this dataset
    cd ../../phylostratr_$current_dataset
done < master_proteomes.csv # replace master_proteomes.csv with test_proteome.csv for small example run

# 4. after blast is done above, rename and move files
# Again, to test with the small example data, replace master_proteomes.csv with test_proteome.csv
cd ../phylostratr_main
while IFS=',' read field1 field2 field3 field4 field5 field6; do
    echo "Field 1: $field1, Field 2: $field2, Field 3: $field3" 

    current_short_name=$field3

    # copy blast results so they're visible for R, then submit R script
    cp phylostratr_${current_dataset}/${current_short_name}/uniprot-seqs/*.tab phylostratr_${current_dataset}/${current_short_name}/

done < master_proteomes.csv # replace master_proteomes.csv with test_proteome.csv for small example run

# submit as array, one for each proteome in the master_proteomes.csv file, e.g. 1-48
sbatch strata_blast_array.sh $current_dataset




