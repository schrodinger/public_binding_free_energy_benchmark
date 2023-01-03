$SCHRODINGER/run sd_analysis_tools.py -i scyt_dehyd/scyt_dehyd_out.fmp -o scyt_dehyd/scyt_dehyd_pkacorr_out.fmp -d macro
$SCHRODINGER/run -FROM scisol binding_mode_correction.py throm_nozob_hip75/throm_nozob_hip75_out.fmp -o throm_nozob_hip75/throm_nozob_hip75_bmcorr_out.fmp
