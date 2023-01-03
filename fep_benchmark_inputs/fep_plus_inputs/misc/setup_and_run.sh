PROTEINS=(btk_extra_flip cdk8_koehler hfaah hiv_prot_ekegren)
for P in ${PROTEINS[*]}
do
    mkdir "${P}"
    cp ${P}.fmp ${P}
    cd "${P}"
    "${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -ffbuilder -ff-host bolt_cpu:50 -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME ${P} ${P}.fmp -QARG "-P dev_GPU"
    cd ../
done

# Error with forcefield builder in 22-3, so using a custom OPLSDIR:
mkdir galectin3_extra
cp -r galectin3_extra.fmp galectin3_ffb_oplsdir galectin3_extra
cd galectin3_extra
"${SCHRODINGER}/fep_plus" -HOST bolt_cpu -SUBHOST bolt_gpu -ppj 2 -OPLSDIR galectin3_ffb_oplsdir -time 20000.0 -vacuum -ensemble muVT -seed 1 -JOBNAME galectin3_extra galectin3_extra.fmp -QARG "-P dev_GPU"
cd ../
