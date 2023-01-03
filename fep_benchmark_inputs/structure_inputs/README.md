# Benchmark structural inputs
The protein and ligand structures of the benchmark. In total, there are 103 protein structures 
associated ligand binding poses that have been seperated in 14 sub-directories that indicate the were the input data was
adapted from. 

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
* `scaffold_hoping/`: FEP+ scaffold hopping
* `waterset/`: FEP+ buried water set

### Subset directory contents
Each directory contains 
* `subset_metadata.csv`. Metadata on each subsystem, the file naming scheme and the PDB each model was based on.  
* `*_ligands.sdf` The ligands along with additional rotomers and protomer/tautomer states. Each ligand is tagged with the
experimental binding free energy (in kcal/mol). 
* `*_protien.pdb` The target structure.
* `*_edges.csv` The edges of the perturbation graph along with experimental relative binding free energy.

