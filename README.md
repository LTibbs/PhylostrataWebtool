# PhylostrataWebtool

This repository contains the code for the Maize Phylostrata Tool at [https://phylostrata.maizegdb.org/](https://phylostrata.maizegdb.org/). It includes code for the analysis underlying the tool, processing the results of that analysis for visualization, and the code for the tool itself.

## maize_phylostratr_analysis
Contains the code used for phylostratr analysis of maize. Based on an updated version of the phylostratr package, see [https://github.com/LTibbs/phylostratr](https://github.com/LTibbs/phylostratr) and [https://doi.org/10.1093/bioinformatics/btz171](https://doi.org/10.1093/bioinformatics/btz171) for more information about this package.

## data_processing
Contains the code used to process phylostratr output to format it for display in the web tool. For users who would like to implement this tool for their own species or use case, sections to customize are annotated with "NOTE" in the comments.

## webtool
In progress! Code used to make the Maize Phylostrata Tool. Again, sections to customize for other species or uses are annotated with "NOTE" in the comments.
