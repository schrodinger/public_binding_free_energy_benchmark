# Contents
### `edge_predictions/`
The collection of CSV files that contain the raw outputs from FEP+ (21-4) using the structures from 'structural_inputs'. 
The columns correspond to
* The names of the ligands involved in the perturbation
* The experimental relative binding free energy
* The raw prediction from FEP+ ("Bennett ddG") and it's estimated error ("Bennett std. error")
* The cycle closure corrected predicted relative binding free energies ("CCC ddG") and their estimated error 
("CCC std. error")

The topology of the perturbation graph is defined by the edges in each CSV file. The direction of each edge is Lig1 -> Lig2.

### `processed_output_fmps/`
The post-calculation and post correction FEP+ `fmp` files that contain the final set of results. These are used to 
calculate the accuracy statistics of FEP+.

Many of the edges in "edge_predictions" are between conformers, protomers, and tautomers of the same molecule. These 
must to processed to produce a single relative binding free energy for each molecule. Each `fmp` file contained in each 
subdirectory should be read with Schrodinger's FEP+ in order to properly account for multiple rotamers and protomers.

### `ligand_predictions/`
The data in these `csv` files have been derived from the `fmp` files in `processed_output_fmps/` and be used to calculate 
the overall accuracy of FEP+ on this benchmark without Schrodinger's tools.

The `csv` files in each subdirectory contain predicted binding free energies after the application of conformer, rotamer symmetry,
and pka/tautomer corrections. After correction, different conformers and tautomers/protomers of the same molecules have 
the same predicted binding free energy. Rather than relative binding free energies, the CSV files contained here have
predicted _absolute_ binding free energies that are known up to an unknown additive constant. The additive constant was
set such that the mean of the predicted absolute binding free energies equals the mean of experimental binding free 
energies.

In each CSV file, the columns refer to
* The name of the ligand
* The experimental binding free energy
* The predicted binding free energy after corrections as well as the estimated standard error.

From these absolute predictions, the pairwise errors and correlation statistics can be calculated. Alternate conformers 
have been removed from the set of the predictions in each `csv`, but alternate protomers and tautomers are still listed. 
If these `csv`s files are used to calculate pairwise errors and correlation statistics without futher processing, the 
ligands with multiple protomers and tautomers will be treated as _seperate_ ligands and thus will err the accuracy metrics.

### `summary_statistics/`
The pairiwise, edgewise, and absolute correlation statistics for each subgroup as `csv` files along with 
`group_summaries.csv`, which contains the aggregated accuracy staitics for each group. 

These `csv` files were created with `../analysis_scripts/write_group_summary_tables.py` using the `fmp` files in 
`processed_output_fmps/` as input. Assuming one has an FEP+ license and the `SCHRODINGER` variable set, one can create
these files with
```python
$SCHRODINGER/run python3 ../analysis_scripts/write_group_summary_tables.py processed_output_fmps -e fmp
```

##### Note
The statistics in `group_summaries.csv` includes a row with the Name `opls_ddag`, which refers to Schrodinger's internal
drug discovery data set. 