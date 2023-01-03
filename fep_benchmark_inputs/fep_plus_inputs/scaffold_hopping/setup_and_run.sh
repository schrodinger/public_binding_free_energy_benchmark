PROTEINS=(Bace1_4zsp CHK1_3u9n_corehop Era_2q70 Fxa_2ei8 TPSB2_3v7t)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

