$SCHRODINGER/run -FROM scisol binding_mode_correction.py syk_4puz_fullmap/syk_4puz_fullmap_out.fmp -o syk_4puz_fullmap/syk_bmcorr_out.fmp
$SCHRODINGER/run -FROM scisol pka_tautomer_correction.py syk_4puz_fullmap/syk_bmcorr_out.fmp -p syk_epik_pops.txt -s- -o syk_4puz_fullmap/syk_pkacorr_out.fmp
