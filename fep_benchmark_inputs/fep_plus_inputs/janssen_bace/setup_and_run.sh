PROTEINS=(bace_ciordia_prospective bace_ciordia_retro bace_keranen_p2)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

mkdir bace_p3_arg368_in
cp bace_p3_arg368_in.fmp bace_p3_arg368_in
cd bace_p3_arg368_in
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -core-hopping-lambda-windows 32 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME bace_p3 bace_p3_arg368_in.fmp -QARG "-P dev_GPU"
cd ..
