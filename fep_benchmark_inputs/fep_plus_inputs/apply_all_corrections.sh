###################
# Bayer macrocycles
###################
echo 'Collecting Bayer macrocycle data'

## Name the outfmp files
echo  "bayer_macrocycles/ftase_extraligs_custcore_stereo/ftase_extraligs_custcore_stereo_out.fmp" > bayer_macrocycles_outs.txt
echo "bayer_macrocycles/wagner_brd4/wagner_brd4_out.fmp" >>  bayer_macrocycles_outs.txt 

#######################
#FEP+ charge change set
#######################
echo 'Collecting FEP+ charge set data'

## Apply the corrections
cd charge_annhil
./run_all_corrections.sh
cd ../

## Name the output fmp file locations
PROTEINS=(cdk2 dlk egfr ephx2 irak4_s2 irak4_s3 itk jak1 jnk1 ptp1b tyk2)
for P in ${PROTEINS[*]}
do
  echo "charge_annhil/${P}/${P}_pkacorr_out.fmp" >> charge_annhil_outs.txt
done
echo "charge_annhil/thrombin_whole_map/thrombin_whole_map_out.fmp" >> charge_annhil_outs.txt

############
# Fragements
############
echo 'Collecting fragement set data'

cd fragments
./run_all_corrections.sh
cd ../

PROTEINS=(frag_liga_auto frag_mcl1_noweak frag_mup1 frag_p38 hsp90_frag_2rings hsp90_frag_single_ring t4lysozyme_uvt)
for P in ${PROTEINS[*]}
do
    echo "fragments/${P}/${P}_out.fmp" >> fragments_outs.txt
done

echo "fragments/jak2_set1/jak2_set1_pkacorr_bmcorr_out.fmp" >> fragments_outs.txt
echo "fragments/jak2_set2_extra/jak2_set2_extra_bmcorr_out.fmp" >> fragments_outs.txt

#######
# GPCRs
#######
echo 'Collecting GPCR data'
# Apply the corrections
cd gpcrs
./run_all_corrections.sh
cd ../

## Name the files
echo "gpcrs/p2y1_merged_out.fmp" > gpcrs_outs.txt
echo "gpcrs/a2a_hip278/a2a_hip278_sbpkacorr_out.fmp" >> gpcrs_outs.txt
echo "gpcrs/ox2_hip_custcore/ox2_out.fmp" >> gpcrs_outs.txt

#####################################
# FEP+ R-group perutbation (JACS set)
#####################################
echo 'Collecting FEP+ R-group perutbation (JACS set) data'
cd jacs_set/
./run_all_corrections.sh
cd ../

echo "jacs_set/jnk1_manual_flips/jnk1_manual_flips_out_bmcorr_out.fmp" > jacs_set_outs.txt
echo "jacs_set/mcl1_extra_flips/mcl1_extra_flips_bmcorr_out.fmp" >> jacs_set_outs.txt

PROTEINS=(bace cdk2 jnk1_manual_flips mcl1_extra_flips p38 ptp1b thrombin_core tyk2)
for P in ${PROTEINS[*]}
do
    echo "jacs_set/${P}/${P}_out.fmp" >> jacs_set_outs.txt
done

################
# Janssen BACE1
###############
echo 'Collecting the Jansseen BACE1 data'

cd janssen_bace 
./run_all_corrections.sh
cd ../

echo "janssen_bace/bace_ciordia_retro/bace_ciordia_retro_pkacorr_out.fmp" > janssen_bace_outs.txt
echo "janssen_bace/bace_keranen_p2/bace_keranen_p2_bmcorr_out.fmp" >> janssen_bace_outs.txt
echo "janssen_bace/bace_ciordia_prospective/bace_ciordia_prospective_out.fmp" >> janssen_bace_outs.txt 
echo "janssen_bace/bace_p3_arg368_in/bace_p3_arg368_in_out.fmp" >> janssen_bace_outs.txt

##################
# FEP+ Macrocycles
##################
echo 'Collecting FEP+ macrocycle data'
PROTEINS=(2B8V_lig24and25_alpha05 2E9P_lig4to7_alpha05 2Q15_lig17to21_alpha05 3RKZ_lig62to70_alpha05 MHT1_lig3_alpha05 ck2_custcore_hotlys hsp90_3hvd_custcore)
for P in ${PROTEINS[*]}
do
    echo "macrocycles/${P}/${P}_out.fmp" >> macrocycle_outs.txt
done

#############
# MCS docking
#############
echo 'Collecting MCS docking data'
cd mcs_docking
./run_all_corrections.sh
cd ../

echo "mcs_docking/hne/hne_500mM_pkacorr_out.fmp" > mcs_docking_outs.txt
echo "mcs_docking/renin_customcore/renin_customcore_out.fmp" >> mcs_docking_outs.txt

#######
# Merck
######
echo 'Collecting Merck data'
cd  merck
./run_all_corrections.sh
cd ../

echo 'merck/cdk8_5cei_new_helix_loop_extra/cdk8_symbmcorr_out.fmp' > merck_outs.txt
echo 'merck/cmet/cmet_out.fmp' >> merck_outs.txt
echo 'merck/eg5_extraprotomers/eg5_extraprotomers_pkacorr_out.fmp' >> merck_outs.txt
echo 'merck/hif2a_automap/hif2a_automap_symbmcorr_out.fmp' >> merck_outs.txt
#echo 'merck/pfkfb3_automap/pfkfb3_automap_symbmcorr_out_fmp' >> merck_outs.txt
echo 'merck/shp2/shp2_out.fmp' >> merck_outs.txt
echo 'merck/syk_4puz_fullmap/syk_pkacorr_out.fmp'  >> merck_outs.txt
echo 'merck/tnks2_fullmap/tnks2_symcorr_pkacorr_out.fmp' >> merck_outs.txt

###############
# Miscellaneous
###############
echo 'Collecting the miscellaneous data sets'
cd misc
./run_all_corrections.sh
cd ../

PROTEINS=(cdk8_koehler galectin3_extra hfaah hiv_prot_ekegren)
for P in ${PROTEINS[*]}
do
    echo "misc/${P}/${P}_out.fmp" >> misc_outs.txt
done
echo 'misc/btk_extra_flip/btk_extra_flip_bmcorr_out.fmp' >> misc_outs.txt

######################
# OPLS drug discovery
#####################
echo 'Collecting the OPLS drug discovery data sets'
PROTEINS=(bathonP_ethers bathonP_thq bathonP_thq_ring hero0 hero1 hero3 hero5 iris lak1 lak2 lak3 orion)
for P in ${PROTEINS[*]}
do
    echo "opls_ddag/${P}/${P}_out.fmp" >> opls_ddag_outs.txt
done

#############
# OPLS Stress
#############
echo 'Collecting the OPLS stress set data'

PROTEINS=(chk1_set1 chk1_set2 chk1_set3 chk1_set4 chk1_set5 chk1_set6 chk1_set7 cr_bace1 cr_bace2 fxa_set3 fxa_set4 fxa_set5 fxa_set6 hc_bace1 hc_bace2 pb_bace3 fxa_yoshikawa_set)
for P in ${PROTEINS[*]}
do
    echo "opls_stress/${P}/${P}_out.fmp" >> opls_stress_outs.txt
done

###################
# Scaffold hoppping
###################
echo 'Collecting the scaffold hopping data set'

PROTEINS=(Bace1_4zsp CHK1_3u9n_corehop Era_2q70 Fxa_2ei8 TPSB2_3v7t)
for P in ${PROTEINS[*]}
do
    echo "scaffold_hopping/${P}/${P}_out.fmp" >> scaffold_hopping_outs.txt
done

##########
# Waterset
##########
echo 'Collecting the water displacement data set'

cd waterset
./run_all_corrections.sh
cd ../

PROTEINS=(hsp90_kung brd41_ASH106 taf12 urokinase hsp90_woodhead)
for P in ${PROTEINS[*]}
do
   echo "waterset/${P}/${P}_out.fmp" >> waterset_outs.txt
done
echo 'waterset/scyt_dehyd/scyt_dehyd_pkacorr_out.fmp' >> waterset_outs.txt
echo 'waterset/throm_nozob_hip75/throm_nozob_hip75_bmcorr_out.fmp' >> waterset_outs.txt

