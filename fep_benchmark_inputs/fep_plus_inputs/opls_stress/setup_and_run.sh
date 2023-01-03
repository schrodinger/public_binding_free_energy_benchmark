PROTEINS=(chk1_set1 chk1_set2 chk1_set3 chk1_set4 chk1_set5 chk1_set6 chk1_set7 cr_bace1 cr_bace2 fxa_set3 fxa_set4 fxa_set5 fxa_set6 hc_bace1 hc_bace2 pb_bace3)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

mkdir fxa_yoshikawa_set
cp fxa_yoshikawa_set.fmp fxa_yoshikawa_set
cd fxa_yoshikawa_set
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -modify_dihe -ensemble muVT -seed 1 -JOBNAME fxa_yoshikawa_set fxa_yoshikawa_set.fmp -QARG "-P dev_GPU"
cd ../
