# 1. Apply the symmetry correction to all ligands except 7.
$SCHRODINGER/run tnks2_symcorr.py tnks2_fullmap/tnks2_fullmap_out.fmp -o tnks2_fullmap/tnks2_symcorr_out.fmp
# 2. Apply the pKa correction
$SCHRODINGER/run -FROM scisol pka_tautomer_correction.py tnks2_fullmap/tnks2_symcorr_out.fmp -pka-file tnks2_macropka.txt -o tnks2_fullmap/tnks2_symcorr_pkacorr_out.fmp -s _
