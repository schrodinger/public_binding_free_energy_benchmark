PROTEINS=(hsp90_kung brd41_ASH106 taf12 scyt_dehyd chk1 throm_nozob_hip75 urokinase hsp90_woodhead)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:30 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done
