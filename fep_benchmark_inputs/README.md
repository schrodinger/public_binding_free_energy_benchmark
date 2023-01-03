# Inputs for running the binding free energy benchmark
The inputs for the free energy benchmark are divided into two sections:
* `fep_plus_inputs/` which contains the inputs for the FEP+ software (e.g. `fmp` files) and the scripts required to
run the benchmark.
* `structure_inputs/` which contain structural inputs in a format can be read by many computational chemistry platforms
(i.e. `PDB` and `sdf` files).

This benchmark data set was constructed and validated using FEP+ and the data in the `structure_inputs` was derived from
the FEP+ `fmp` graph files..

## Arrangment of the input files
The benchmark data set consists of a number of congeneric series that have been grouped based on the publication where 
the inputs were based on. As documented in the supporting information of the [accompanying manuscript](https://doi.org/10.26434/chemrxiv-2022-p2vpg),
many systems were modified from their original source. 

## Metadata
Details on the congeneric series subset, file naming scheme, protein name, original PDB reference, experimental binding
affinity dynamic range, number of input structures, and the DOI of the source on the source where the structures were 
taken from, can be found in the file below:
* `benchmark_metadata.csv`

The `Number of nodes` column in `benchmark_metadata.csv` refers to the number of ligand structural inputs. Some ligands 
had multiple structural inputs as they were represented by multiple rotamers or protomer/tautomers.