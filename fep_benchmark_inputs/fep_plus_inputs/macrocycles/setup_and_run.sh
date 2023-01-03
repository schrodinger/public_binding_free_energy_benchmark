PROTEINS=(2B8V_lig24and25_alpha05 2E9P_lig4to7_alpha05 2Q15_lig17to21_alpha05 3RKZ_lig62to70_alpha05 MHT1_lig3_alpha05)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

PROTEINS=(ck2_custcore_hotlys hsp90_3hvd_custcore)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -modify_dihe -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done
