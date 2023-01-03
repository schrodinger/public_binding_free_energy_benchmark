PROTEINS=(a2a_hip278 p2y1_meta_sub p2y1_ortho_sub)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

mkdir ox2_hip_custcore
cp ox2_hip_custcore.fmp ox2_hip_custcore
cd ox2_hip_custcore
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -modify_dihe -ensemble muVT -seed 1 -JOBNAME ox2 ox2_hip_custcore.fmp -QARG "-P dev_GPU"
cd ../
