## Benchmark subsets with systems that have post simulation corrections
Here is a breakdown of the corrections applied to the maps. This can be applied automatically wtih `apply_all_corrections.sh`.

Prior to the 22-3 release, pKa/tauotmer and binding mode corrections had to be run seperately from the command line and
custom scripts. The corrections outlined in this section are for releases 22-2 and earlier and can be streamlined for 
more recent releases. 

FEP maps and directories that are _not_ listed below require no post processing.
 
###  `charge_annhil/` FEP+ charge-change 
Every `fmp` file, except `thrombin_whole_map.fmp` requires a pKa correction. The solvent pKas can be found in the accompanying 
`.txt` files. The are 2 sets of solvent pKas: manually chosen (`*_pka.txt`)or automatically entered with epiK 
(`*_epikpka.txt`). These pKa files can be used as input to the `pka_tautomer_correction.py` script.

### `fragments/` FEP+ fragements
There are 2 jak2 maps that require pka/tautomer corrections as well as binding mode corrections. These correction 
scripts are:
* `run_jak2_set1_corrs.sh`, which uses `jak2_set1_pka.txt` as input.
* `run_jak2_set2_corrs.sh`
 
### `gpcrs/` GPCRs
* A2A: pka, binding mode, and symmetry corrections can be run with `run_a2a_corrs.sh` which uses `a2a_epik_pka.txt` and 
`a2a_symbmcorr.py`.
* The output 2 P2Y1 maps can be merged using: 
    - `$SCHRODINGER/run -FROM scisol merge_graph.py p2y1_meta_sub_out.fmp p2y1_ortho_sub_out.fmp -o p2y1_merged.fmp -sc`
    

### `jacs_set/` FEP+ R-group set
* The JNK1 output requires binding mode corrections: 
    -  `$SCHRODINGER/run -FROM scisol binding_mode_correction.py jnk1_manual_flips_out.fmp -o jnk1_manual_flips_out_bmcorr_out.fmp`
* The MCL1 map also requires binding mode corrections:
    - `$SCHRODINGER/run -FROM scisol binding_mode_correction.py mcl1_extra_flips_out.fmp -o mcl1_extra_flips_out_bmcorr_out.fmp`
    
    
### `janssen_bace/` Janssen BACE1 data sets
* One molecule from the retrospective set from Ciordai et al. requires a pKa correction for one molecule. This can be 
run using `run_bace_ciordia_retro_corrs.sh`
* The P2 pocket subset from Keranen et al. requires binding mode corrections. These can be run using 
`run_bace_keranen_p2_corrs.sh`

### `mcs_docking` MCS docking sets
* The HNE map has some molecules with multiple protonation states. These can be corrected for using
    - `$SCHRODINGER/run -FROM scisol pka_tautomer_correction.py hne_500mM_out.fmp -pka-file hne_epik_pka.txt -o hne_500mM_pkacorr_out.fmp`

### `merck/` The public Merck data sets
* The CDK8 map requires symmetry and binding mode corrections. These can be implemented with `cdk8_symbmcorr.py` script
like so:
    -  `$SCHRODINGER/run cdk8_symbmcorr.py cdk8_5cei_new_helix_loop_extra_out.fmp -o cdk8_symbmcorr_out.fmp`
* The Eg5 map has multiple protomers that can be handled usig `run_eg5_corrections.sh`.
* Some ligands in the hif2-alpha map have 2 rotamer states whereas other require symmetry corrections. These can be 
handled with the python script:
    - `$SCHRODINGER/run hif2a_symbmcorr.py hif2a_automap_out.fmp -o hif2a_automap_symbmcorr.fmp`
* Most ligands in the PFKFB3 have multiple rotamers whereas some require a symmetry correction, which require the 
application of the script:
    - `$SCHRODINGER/run pfkfb3_symbmcorr.py pfkfb3_automap.fmp -o pfkfb3_automap_symbmcorr_out_fmp`
* Some ligands in the SYK map have different protomers, tautomers, and rotamers, which can be addressed using `run_syk_corrs.sh`.
* The TNKS2 requires correcting for multiple protomers and one symmetric rotamer state. These can be corrected using
` run_tnks2_corrs.sh`.

### `misc/` Miscellaneous
* Some of the ligands in the BTK map have 2 rotamer states that have been added to the map. These can be accounted for
with 
    - `$SCHRODINGER/run binding_mode_correction.py btk_extra_flip_out.fmp -o btk_extra_flip_bmcorr_out.fmp`

### `waterset/` FEP+ buried water set
* The scytalone dehydratase map has a bespoke pKa correction script that was written for the validation of the FEP+'s
GCMC method. The corrections can be applied using
    - `$SCHRODINGER/run sd_analysis_tools.py -i scyt_dehyd_out.fmp -o scyt_dehyd_pkacorr_out.fmp -d macro`

