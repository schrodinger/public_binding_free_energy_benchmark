$SCHRODINGER/run -FROM scisol pka_tautomer_correction.py jak2_set1_out.fmp -p jak2_set1_pka.txt -o jak2_set1_pkacorr_out.fmp
$SCHRODINGER/run -FROM scisol binding_mode_correction.py jak2_set1_pkacorr_out.fmp -o jak2_set1_pkacorr_bmcorr_out.fmp

