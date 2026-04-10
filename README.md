# PhylostrataWebtool

This repository contains the code for the MaizeGDB Phylostrata Tool at [https://phylostrata.maizegdb.org/](https://phylostrata.maizegdb.org/). It includes code for the analysis underlying the tool, processing the results of that analysis for visualization, and the code for the tool itself. This code can also be used to create a similar web tool for other species.

Read the pre-print at https://doi.org/10.64898/2025.12.19.695500.

## maize_phylostratr_analysis
Contains the code used for phylostratr analysis of maize, as well as description of input and output files for an example run with test data. Analysis based on an updated version of the phylostratr package, see [https://github.com/arendsee/phylostratr](https://github.com/arendsee/phylostratr) or [https://github.com/LTibbs/phylostratr](https://github.com/LTibbs/phylostratr) and [https://doi.org/10.1093/bioinformatics/btz171](https://doi.org/10.1093/bioinformatics/btz171) for more information about this package.  

## data_processing
Contains the code used to process phylostratr output to format it for display in the web tool. For users who would like to implement this tool for their own species or use case, sections to customize are annotated with "NOTE" in the comments.

## webtool
Code used to make the MaizeGDB Phylostrata Tool. Again, sections to customize for users who would like to make a version of this web tool for other species are annotated with "NOTE" in the comments.
