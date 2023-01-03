cp jak2_set1_pka.txt run_jak2_set1_corrs.sh jak2_set1
cd jak2_set1
./run_jak2_set1_corrs.sh
cd ../

cd jak2_set2_extra
cp ../run_jak2_set2_corrs.sh .
./run_jak2_set2_corrs.sh
rm run_jak2_set2_corrs.sh
cd ../

