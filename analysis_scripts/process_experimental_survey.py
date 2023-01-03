import analysis_functions as af
import argparse

def main(argv=None):
    usage = """
        As input, this script requires the directory that contains the CSV files of the individual experimental binding
        free energy comparisons. This data has been taken from publicly assessible sources and the CSV filenames have
        been hard-coded in this script. 
         
        > python process_experimental_survey.py directory1
            
        Optionally, a second directory can also be supplied. This second directory is expected to contain the comparative
        binding free energy data from Schrodinger's drug discovery projects. 
        
        > python process_experimental_survey.py directory1 --drug_discovery_dir directory2
        """
    description = """
        Estimate the degree of reproducibility of experimental measured binding free energies using comparative assay
        data. This script process chemical series that have had their binding free energies measured at least two 
        different methods. 
        
        Each assay comparison has been seperated into 3 catogories: those where two binding assays are compared, where
        a binding assay and a functional assay are compared, and where two functional assays are compared. This script 
        produces reproducibility estimates for each category. 
        """
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument(
        'dir',
        type=str,
        help="The directory that contains the publicly accessible experimental survey data.")
    parser.add_argument(
        '-d',
        '--drug_discovery_dir',
        type=str,
        help="The directory that contains Schrodinger's drug discovery comparative assay data, default=None", default=None)
    args = parser.parse_args(argv)

    dir1 = args.dir
    dir2 = args.drug_discovery_dir
    # Splitting the survey into the different catagories of comparison
    binding_comparisons = [f'{dir1}/aaron2010_spr_itc.csv',
                           f'{dir1}/hang2009_hcvpoly_con1_spr_fluor.csv',
                           f'{dir1}/jecklin2009_ms_itc.csv',
                           f'{dir1}/jecklin2009_ms_spr.csv',
                           f'{dir1}/jecklin2009_spr_itc.csv',
                           f'{dir1}/mandine2001_sh2_spa_spr.csv',
                           f'{dir1}/mason2012_spr_lan.csv',
                           f'{dir1}/murphy2006_binding.csv',
                           f'{dir1}/murthy2007_itc_spr.csv',
                           f'{dir1}/navratilova2007_caII_itc_spr.csv',
                           f'{dir1}/newman2012_spr_itc.csv',
                           f'{dir1}/peterson2018_itc_fp.csv',
                           f'{dir1}/rogez2013_bovine_itc_tsa.csv',
                           f'{dir1}/rogez2013_bovine_spr_itc.csv',
                           f'{dir1}/rogez2013_bovine_spr_tsa.csv',
                           f'{dir1}/rogez2013_human_itc_tsa.csv',
                           f'{dir1}/rogez2013_human_spr_itc.csv',
                           f'{dir1}/rogez2013_human_spr_tsa.csv',
                           f'{dir1}/schnapp2016_dppiv_itc_spr.csv',
                           f'{dir1}/ycas2020_bptf_labeled_spr_alphascreen.csv',
                           f'{dir1}/ycas2020_bptf_nmr_alphascreen.csv',
                           f'{dir1}/ycas2020_bptf_nmr_labeled_spr.csv',
                           f'{dir1}/ycas2020_bptf_nmr_spr.csv',
                           f'{dir1}/ycas2020_bptf_spr_alphascreen.csv']

    if dir2 is not None:
        binding_comparisons.extend([f'{dir2}/projectD_lantha_discover.csv', f'{dir2}/projectE_lantha_discover.csv'])

    binding_vs_functional_comparisons = [f'{dir1}/baum2009_thrombin_ki_itc.csv',
                                         f'{dir1}/hang2009_hcvpoly_bk_ic50_fluor.csv',
                                         f'{dir1}/hang2009_hcvpoly_con1_ic50_fluor.csv',
                                         f'{dir1}/kung2011_hsp90_spa_itc.csv',
                                         f'{dir1}/li2016_dppiv_spr_ic50.csv',
                                         f'{dir1}/mason2012_lan_ic50.csv',
                                         f'{dir1}/mason2012_spr_ic50.csv',
                                         f'{dir1}/nikolovska2004_xiap_ki_kd.csv',
                                         f'{dir1}/crawford2016_binding_BRD41.csv',
                                         f'{dir1}/crawford2016_binding_BRD42.csv',
                                         f'{dir1}/crawford2016_binding_BRD9.csv',
                                         f'{dir1}/crawford2016_binding_BRPF1.csv',
                                         f'{dir1}/crawford2016_binding_CECR2.csv',
                                         f'{dir1}/crawford2016_binding_CREBBP.csv',
                                         f'{dir1}/crawford2016_binding_TAF12.csv',
                                         f'{dir1}/patil2018_ic50_spr_nottopbottom.csv',
                                         f'{dir1}/schnapp2016_dppiv_ic50_spr.csv',
                                         f'{dir1}/schnapp2016_dppiv_itc_ic50.csv',
                                         f'{dir1}/schindler2020_functional_spr.csv']

    if dir2 is not None:
        binding_vs_functional_comparisons.extend([ f'{dir2}/projectA_spr_biochem.csv',
                                         f'{dir2}/projectB_biochem_phospho.csv',
                                         f'{dir2}/projectC_spr_biochem.csv',
                                         f'{dir2}/projectD_discover_trfret.csv',
                                         f'{dir2}/projectD_discover_atpkm.csv',
                                         f'{dir2}/projectD_lantha_atpkm.csv',
                                         f'{dir2}/projectD_lantha_trfret.csv',
                                         f'{dir2}/projectE_lantha_atpkm.csv',
                                         f'{dir2}/projectE_lantha_trfret.csv',
                                         f'{dir2}/projectE_discover_trfret.csv',
                                         f'{dir2}/projectE_discover_atpkm.csv'])

    functional_comparisons = [f'{dir1}/chen2013_angiotensin_inhibition.csv',
                              f'{dir1}/jia2006_cot_inhibition_nottopbottom.csv',
                              f'{dir1}/katz2003_upa_competition.csv',
                              f'{dir1}/moonshot2020_covid_protease_inhibition_nottopbottom.csv']
    if dir2 is not None:
        functional_comparisons.extend([f'{dir2}/projectD_atpkm_trfret.csv',
                              f'{dir2}/projectE_atpkm_trfret.csv'])

    all_comparitive_files = binding_comparisons + binding_vs_functional_comparisons + functional_comparisons

    print('The overall experimental error in the survey')
    print('--------------------------------------------')
    af.summarize_experimental_error(all_comparitive_files)

    print('Biophysical vs biophysical error')
    print('--------------------------------')
    print(f'Number of comparisons = {len(binding_comparisons)}')
    af.summarize_experimental_error(binding_comparisons)

    print('Biophysical vs biochemical error')
    print('---------------------------------')
    print(f'Number of comparisons = {len(binding_vs_functional_comparisons)}')
    af.summarize_experimental_error(binding_vs_functional_comparisons)

    print('Biochemical vs biochemical error')
    print('---------------------------------')
    print(f'Number of comparisons = {len(functional_comparisons)}')
    af.summarize_experimental_error(functional_comparisons)

if __name__== '__main__':
    main()



