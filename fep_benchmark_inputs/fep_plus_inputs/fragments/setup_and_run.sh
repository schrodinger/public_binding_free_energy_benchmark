PROTEINS=(frag_liga_auto frag_mcl1_noweak frag_mup1 frag_p38 hsp90_frag_2rings hsp90_frag_single_ring jak2_set1 jak2_set2_extra t4lysozyme_uvt)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

