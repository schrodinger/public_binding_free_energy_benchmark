$SCHRODINGER/run cdk8_symbmcorr.py cdk8_5cei_new_helix_loop_extra/cdk8_5cei_new_helix_loop_extra_out.fmp -o cdk8_5cei_new_helix_loop_extra/cdk8_symbmcorr_out.fmp
./run_eg5_corrections.sh
$SCHRODINGER/run hif2a_symbmcorr.py hif2a_automap/hif2a_automap_out.fmp -o hif2a_automap/hif2a_automap_symbmcorr_out.fmp
#$SCHRODINGER/run pfkfb3_symbmcorr.py pfkfb3_automap/pfkfb3_automap.fmp -o pfkfb3_automap/pfkfb3_automap_symbmcorr_out_fmp
./run_syk_corrs.sh 
./run_tnks2_corrs.sh
