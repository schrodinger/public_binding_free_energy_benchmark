
mkdir renin_customcore
cp renin_customcore.fmp renin_customcore
cd renin_customcore
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ffbuilder -ff-host bolt_cpu:50 -ppj 2 -time 20000.0 -ensemble muVT -seed 1 -vacuum -lambda_windows 12 -core_hopping_lambda_windows 16 -charged_lambda_windows 24 -JOBNAME renin_customcore renin_customcore.fmp -TMPLAUNCHDIR -QARG "-P dev_GPU" -modify-dihe
cd ../

mkdir hne
cp hne.fmp hne
cd hne
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ffbuilder -ff-host bolt_cpu:50 -ppj 2 -time 20000.0 -ensemble muVT -salt 0.5 -seed 1 -vacuum -lambda_windows 12 -core_hopping_lambda_windows 16 -charged_lambda_windows 24 -JOBNAME hne_500mM hne.fmp -TMPLAUNCHDIR  -QARG "-P dev_GPU"
cd ../ 
