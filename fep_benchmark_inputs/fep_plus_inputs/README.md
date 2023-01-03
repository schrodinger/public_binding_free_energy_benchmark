# Benchmark inputs for FEP+

This directory contains the FEP+ input files (`.fmp`s) and submission scripts for running the extended benchmark data
set. There are 103 `fmp` files for a wide range of targets and congeneric series. The input 
files are separated into 14 separate directories that reflect the origin of the data sets. Every `fmp` file in this 
directory forms is part of the benchmark. 

## Contents
* `apply_all_corrections.sh` Apply post simulation corrections to all output fmp files that require them. This 
script also collects the names and locations of the output fmp files into text files which facilitate analysis.
* `post_simulation_corrections.md` A description of the post simulation corrections. These corrections are 
applied with `apply_all_corrections.sh`.
* `binding_mode_correction.py` The script to automatically apply binding mode corrections. This is needed for some
systems.
* `README.md`. This file. 

## Subset directories
The names of each directory and the corresponding subset title in the [accompanying manuscript](https://doi.org/10.26434/chemrxiv-2022-p2vpg).
* `bayer_macrocycles/`: Bayer Macrocycles
* `charge_annhil/`: FEP+ charge change set
* `fragments/`: FEP+ fragment data set
* `gpcrs/`: GPCRs
* `jacs_set/`: FEP+ R-group set
* `janssen_bace/`: Janssen BACE1 data sets
* `macrocycles/`: FEP+ macrocycles
* `mcs_docking/`: MCS docking sets
* `merck/`: The public Merck data set
* `misc/`: Miscellaneous data sets
* `opls_stress/`: OPLS stress set
* `opls_ddag`: OPLS private drug discovery data sets
* `scaffold_hoping/`: FEP+ scaffold hopping
* `waterset/`: FEP+ buried water set

### Subset directory contents
Each directory has the following files
* `subset_metadata.csv`. Metadata on each subsystem, the naming scheme and the PDB each model was based on.  
* `*.fmp` files. These are FEP+ format maps/perturbation graphs
* `setup_and_run.sh`. A submission script for all FEP+ maps in the that uses the settings described in the accompanying paper. 
* ` run_all_corrections.sh`. A script to run all corrections on the fmp files in that directory.
* Additional files related to post-processing corrections to the results.


