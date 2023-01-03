# FEP+ benchmark and experimental reproduceability survey
A data set of proteins and congeneric series of ligands for the validation of structure-based binding free energy 
prediction methods and a survey of experimental reproduceability. The data in this repository was used in the 
[accompanying manuscript](https://doi.org/10.26434/chemrxiv-2022-p2vpg).

## Contents

### `fep_benchmark_inputs/`
Input structures and files for the data set for binding free energy calculations. The input files have been set up for 
two use cases: for running with Schrodinger's FEP+ program (`fep_benchmark_inputs/fep_plus_inputs/`) and for running
with any other software (`fep_benchmark_inputs/structure_inputs/`). 

Metadata on the origin of each data set are included in each subdirectory. 

### `experimental_survey_data/`
The experimental binding data extracted from studies that used at least two different assays to measure the binding 
affinity of the same set of compounds to the same protein. These data were used to estimate to experimental 
reproducibility of binding and functional assays.

### `21_4_results/`
Output files from running Schrodingers 21-4 release of FEP+ on the files in `fep_benchmark_inputs/fep_plus_inputs/`. The data is 
included as Schrodinger compatible files as well as `CSV` files.

### `analysis_scripts/`
Python tools and command-line scripts that were used to generate the results in the [accompanying manuscript](https://doi.org/10.26434/chemrxiv-2022-p2vpg).