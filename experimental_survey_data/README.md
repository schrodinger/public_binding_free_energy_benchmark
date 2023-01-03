# Experimental survey data
Comparative assay data where at least two different assays were used to measure the binding affinities of the same set 
of compounds to the same protein. 
The focus of this survey is on the _reproducibility_ of the measured binding free energy, as opposed to the 
_repeatability_ of a particular assay. 
To be included in this survey, affinities had to be reported as dissociation constants 
(Kd), inhibition constants, (Ki) or the ligand concentration that achieve fifty percent inhibition (IC50).

## Contents
Binding affinity comparisons are stored as CSV files in the directory 
* `publicly_accessible_survey_data/`: the collection of CSV files from publicly accessible sources. Each CSV file contains two sets of 
measurements for the same chemical series; the columns correspond to the measurements from a different assay and the 
rows correspond to different ligands. 

The metadata of the CSV files in `publicly_accessible_survey_data/` is contained in
* `publicly_accessible_survey_metadata.csv`, which lists information such as the protein name, the assays types, the DOI
of the where the data was taken from as well as where in the source document where the affinities can be found.


### Note
The experimental data used in this [manuscript](https://doi.org/10.26434/chemrxiv-2022-p2vpg) included data from 
Schrodinger's drug discovery programs which is not included here.