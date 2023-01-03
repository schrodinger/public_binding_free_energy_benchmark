mkdir wagner_brd4
cp wagner_brd4.fmp wagner_brd4
cd wagner_brd4
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME wagner_brd4 wagner_brd4.fmp -QARG "-P dev_GPU"
cd ../

mkdir ftase_extraligs_custcore_stereo
cp ftase_extraligs_custcore_stereo.fmp ftase_extraligs_custcore_stereo
cd ftase_extraligs_custcore_stereo
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 40000.0 -vacuum -ensemble muVT -core_hopping_lambda_windows 48 -seed 1 -JOBNAME ftase_extraligs_custcore_stereo ftase_extraligs_custcore_stereo.fmp -QARG "-P dev_GPU"
cd ../

