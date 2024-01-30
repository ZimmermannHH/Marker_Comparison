# Marker_Comparison

Scripts for in silico PCRs and output created for the comparison of short markers for DNA metabarcoding of sedimentary ancient DNA (sedaDNA).

Here, I provide the code that was used for the processing, filtering, and analysis of the markers for the paper by [Zimmermann et al.](https://www.nature.com/articles/s41467-023-36845-x) Marine ecosystem shifts with deglacial sea-ice loss inferred from ancient DNA shotgun sequencing. Nat Commun 14, 1650 (2023). https://doi.org/10.1038/s41467-023-36845-x

The code to process the data is explained exemplary in this README and bash scripts are stored in the `scripts` folder of this repository. Further data processing, filtering and the analysis of the data was carried out in `R`. The R-scripts are part of this repository and the README guides through the order of running them.

For each marker, the folder `species_lists` contains a list of species that can be resolved on species level. 

To run the scripts, the file `marker.lineages.df.txt` needs to be downloaded from GEUS Dataverse https://doi.org/10.22008/FK2/A8PSA1

Run scripts in the following order:
1. marker_analysis.R
2. resolution_coverage.R
3. barcode_resolution_coverage_supplement.R 
4. amplicon_size.R
5. phylo_tree.R

