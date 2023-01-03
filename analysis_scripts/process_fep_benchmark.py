import analysis_functions as af
import argparse
import importlib
from glob import glob
import os
def main(argv=None):
    usage = """
    The script takes output FMP or CSV files and returns aggregate accuracy statistics. The FMP or CSV files must be 
    grouped in subdirectories of the input leading directory. For example, if there are 2 groups with 3 output files 
    with the following directory structure:
    
        upper_dir
            ├── group_1 
                  ├── system_1_out.fmp (or .csv)
                  ├── system_2_out.fmp (or .csv)
                  ├── system_3_out.fmp (or .csv)
            ├── group_2 
                  ├── system_A_out.fmp (or .csv)
                  ├── system_B_out.fmp (or .csv)
                  ├── system_C_out.fmp (or .csv)
            
    then this script can be run with
    
        > $SCHRODINGER/run python3 process_fep_benchmark.py upper_dir -e fmp
        
    The aggregate statistics for the FEP+ manuscript were calculated using the following command
        
        > $SCHRODINGER/run python3 process_fep_benchmark.py ../21_4_results/processed_output_fmps -e fmp 
    
    Alternatively, the aggregate statistics can be estimated with the CSVs in ../21_4_results/ligand_predictions
        
        > python process_fep_benchmark.py ../21_4_results/ligand_predictions -e csv
    
    Using the CSVs files in ../21_4_results/ligand_predictions only provides approximately accurate statistics as 
    ligands with multiple protomers or tautomers are over-counted. Ligands with multiple protomers or tautomers are 
    properly accounted for when using FMP files.
    """
    description = """
    Assess the accuracy of the FEP+ benchmark using the same metrics as discussed in the maximal and current accuracy of 
    rigorous protein-ligand binding free energy calculations".
    
    The overall accuracy is calculated as well as the break down per group of fmps with the --latex_tables flag.
    """
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument(
        'upper_dir',
        type=str,
        help="The upper directory where all the results are contained. FEP accuracy can be broken by "
             "group. Each group is expected to be a subdirectory of 'upper-dir' and each subdirectory should contain "
             "either FMP or CSV files.")
    parser.add_argument(
        '-e',
        '--ext',
        type=str,
        choices=['fmp', 'csv'],
        help="The file extension of the results. Results can be either FMP files or CSVs.")
    args = parser.parse_args(argv)

    if args.ext == 'fmp':
        if importlib.util.find_spec('schrodinger') is None:
            raise Exception('Schrodinger must be installed to use FMP files as input. Use CSV files instead.')

    files = []
    for entry in glob(f'{args.upper_dir}/*'):
        if os.path.isdir(entry):
            files.extend(glob(f'{entry}/*{args.ext}'))

    # Get the analysis metrics for each map
    if args.ext == 'fmp':
        results, diffs = af.parse_fep_data(files)
    elif args.ext == 'csv':
        results, diffs = af.parse_fep_data_from_csv(files)
    else:
        raise Exception(f'Only "fmp" and "csv" are accessible file extenstions. You have entered {args.ext}.')

    # Correcting for the ligands in the thrombin water displacement set map that overlap with thrombin JACS set map.
    throm_name_waterset = 'throm_nozob_hip75_sbmcorr_out'
    throm_name_jacsset = 'thrombin_core_out'
    for i, name in enumerate(results['entries']):
        if name == throm_name_waterset:
            throm_name_waterset_ind = i
        if name == throm_name_jacsset:
            throm_name_jacsset_ind = i

    if throm_name_waterset_ind is not None and throm_name_jacsset_ind is not None:
        results['number of compounds'][throm_name_waterset_ind] -= results['number of compounds'][throm_name_jacsset_ind]

    print('FEP+ benchmark summary')
    print('-----------------------')
    af.summarize_fep_error(results)
    print()


if __name__== '__main__':
    main()
