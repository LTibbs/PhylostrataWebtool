# Maize Phylostratr Analysis
This folder contains the scripts and input data files used to set up and run the phylostratr analysis of maize. The results of these analyses form the basis for the Maize Phylostrata Tool at https://phylostrata.maizegdb.org/.

1. [maize_phylostratr_setup.R](https://github.com/LTibbs/PhylostrataWebtool/blob/main/maize_phylostratr_analysis/maize_phylostratr_setup.R):
This R script is used to set up for the phylostratr analysis.

2. [phylostratr_pipeline.sh](https://github.com/LTibbs/PhylostrataWebtool/blob/main/maize_phylostratr_analysis/phylostratr_pipeline.sh):
This bash script is used to organize and submit the sub-scripts needed to actually run phylostratr.

3. [diamond_blast.sh](https://github.com/LTibbs/PhylostrataWebtool/blob/main/maize_phylostratr_analysis/diamond_blast.sh)
This bash script is referenced by phylostratr_pipeline.sh above and is used to run Diamond Blastp for each query proteome against all other proteomes. This is used to speed up phylostratr's strata_blast().

4. [strata_blast_array.sh](https://github.com/LTibbs/PhylostrataWebtool/blob/main/maize_phylostratr_analysis/strata_blast_array.sh)
This bash script is used to submit the phylostratr jobs as an array, here --array=1-48, where each represents a single proteome used as the query.

5. [strata_blast_array.R](https://github.com/LTibbs/PhylostrataWebtool/blob/main/maize_phylostratr_analysis/strata_blast_array.R)
This R script is submitted by strata_blast_array.sh above, and runs the phylostratr analysis for each query proteome. 
