PROTEINS=(cdk8_5cei_new_helix_loop_extra cmet eg5_extraprotomers hif2a_automap pfkfb3_automap shp2 syk_4puz_fullmap tnks2_fullmap)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done
