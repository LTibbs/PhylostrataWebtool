# Data Processing

The goal of the code in this subfolder is to process the `Rdata` outupt from Phylostratr and combine it with UniProt data (subcellular localization and GO terms) to create the `json` files used for the webtool visualization. The goal is to end with two json files:

1. The `uniprot_table_loc.json` file is used to build the details page for each protein. It has the best hit for each protein within each stratum, as well as other details including subcellular localization and GO terms in the chosen example species for each stratum. For an example of this file format, see `uniprot_table_loc_mini.json`.
2. The `gene_fill.json` file is used to build the visual search result image showing the stratum level of each protein. It has the stratum of each protein as well as the fill percentage for this protein, which is calculated based on the stratum number. For an example of this file format, see `gene_fill_mini.json`.

In this code, I have used `NOTE` to denote places where the user may want or need to make changes to customize the results. For example, the user can specify the example species for each phylostratum or choose to combine multiple strata into one. 

