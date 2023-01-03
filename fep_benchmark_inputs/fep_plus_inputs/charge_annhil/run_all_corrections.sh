PROTEINS=(cdk2 dlk egfr ephx2 irak4_s2 irak4_s3 itk jak1 jnk1 ptp1b tyk2)
for P in ${PROTEINS[*]}
do
    cd "${P}"
    $SCHRODINGER/run -FROM scisol pka_tautomer_correction.py "${P}"_out.fmp -pka-file ../"${P}"_pka.txt -o "${P}"_pkacorr_out.fmp    
    cd ../
done

