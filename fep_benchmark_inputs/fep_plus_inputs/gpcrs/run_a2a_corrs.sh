$SCHRODINGER/run a2a_symbmcorr.py a2a_hip278/a2a_hip278_out.fmp -o a2a_hip278/a2a_hip278_symbmcorr_out.fmp
$SCHRODINGER/run -FROM scisol pka_tautomer_correction.py a2a_hip278/a2a_hip278_symbmcorr_out.fmp -o a2a_hip278/a2a_hip278_sbpkacorr_out.fmp -s_ -pka-file a2a_epik_pka.txt

