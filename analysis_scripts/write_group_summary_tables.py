import analysis_functions as af
import argparse
import importlib
from glob import glob
import os
import pandas as pd

def main(argv=None):
    usage = """
    The script expects that all output fmp or csv files are located in subdirectories of the supplied main directory. 
    For example, if there are 2 groups with 3 output files with the following directory structure:
    
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
    
        > $SCHRODINGER/run python3 write_group_summary_tables.py upper_dir -e fmp

    If you would like to generate the grouped statistics for the output FEP files in this repository, please type
        
        > $SCHRODINGER/run python3 write_group_summary_tables.py ../21_4_results/processed_output_fmps -e fmp
    
    or, if you do not have a Schrodinger installation, you can use CSV files with:
    
        > python write_group_summary_tables.py ../21_4_results/ligand_predictions -e csv
    
    NOTE: using the CSVs files in ../21_4_results/ligand_predictions only provides approximately accurate statistics as 
    ligands with multple protomers or tautomers are over-counted.
    
    The output CSVs are written to the working directory.
    """
    description = """
    Write out the per-group FEP+ accuracy results to CSV files as well as the summary accuracy for each group.
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

    group_summaries = []
    for entry in glob(f'{args.upper_dir}/*'):
        if os.path.isdir(entry):
            group_name = (entry.split('/')[-1])
            files = glob(f'{entry}/*{args.ext}')
            if args.ext == 'fmp':
                results, diffs = af.parse_fep_data(files)
            elif args.ext == 'csv':
                results, diffs = af.parse_fep_data_from_csv(files)
            else:
                raise Exception(f'Only "fmp" and "csv" are accessible file extenstions. You have entered {args.ext}.')
            df = pd.DataFrame(results)
            df.to_csv(f'{group_name}_results.csv', index=False, float_format='%.2f')
            if args.ext == 'fmp':
                summary = [group_name,
                           results['number of compounds'].sum(),
                           results['number of edges'].sum(),
                           af.summarize_fep_error(results, verbose=False)]
            else:
                summary = [group_name,
                           results['number of compounds'].sum(),
                           af.summarize_fep_error(results, verbose=False)]
            group_summaries.append(summary)

    # Now write out the summary table for each group. Every stat has confidence intervals calculated by boostrap
    # sampling.
    if args.ext == 'fmp':
        group_results = {'Name':[], 'No. compounds':[], 'No. edges':[],
                         'Pairwise MUE':[], 'Pairwise MUE, lower 95%':[], 'Pairwise MUE, upper 95%':[],
                         'Pairwise RMSE':[],'Pairwise RMSE, lower 95%':[], 'Pairwise RMSE, upper 95%':[],
                         'Edgewise MUE':[], 'Edgewise MUE, lower 95%':[],'Edgewise MUE, upper 95%':[],
                         'Edgewise RMSE':[], 'Edgewise RMSE, lower 95%':[],'Edgewise RMSE, upper 95%':[],
                         'R-squared': [], 'R-squared, lower 95%':[], 'R-squared, upper 95%':[],
                         'Kendall tau':[], 'Kendall tau, lower 95%':[], 'Kendall tau, upper 95%':[]}
    else:
        group_results = {'Name':[], 'No. compounds':[],
                         'Pairwise MUE':[], 'Pairwise MUE, lower 95%':[], 'Pairwise MUE, upper 95%':[],
                         'Pairwise RMSE':[],'Pairwise RMSE, lower 95%':[], 'Pairwise RMSE, upper 95%':[],
                         'R-squared': [], 'R-squared, lower 95%':[], 'R-squared, upper 95%':[],
                         'Kendall tau':[], 'Kendall tau, lower 95%':[], 'Kendall tau, upper 95%':[]}

    for summary in group_summaries:
        group_results['Name'].append(summary[0])
        group_results['No. compounds'].append(summary[1])
        if args.ext == 'fmp':
            group_results['No. edges'].append(summary[2])
            pair_rmse, pair_mue, edge_rmse, edge_mue, r2, tau = summary[3]
        else:
            pair_rmse, pair_mue, r2, tau = summary[2]
        group_results['Pairwise MUE'].append(pair_mue[0])
        group_results['Pairwise MUE, lower 95%'].append(pair_mue[1])
        group_results['Pairwise MUE, upper 95%'].append(pair_mue[2])
        group_results['Pairwise RMSE'].append(pair_rmse[0])
        group_results['Pairwise RMSE, lower 95%'].append(pair_rmse[1])
        group_results['Pairwise RMSE, upper 95%'].append(pair_rmse[2])
        if args.ext == 'fmp':
            group_results['Edgewise MUE'].append(edge_mue[0])
            group_results['Edgewise MUE, lower 95%'].append(edge_mue[1])
            group_results['Edgewise MUE, upper 95%'].append(edge_mue[2])
            group_results['Edgewise RMSE'].append(edge_rmse[0])
            group_results['Edgewise RMSE, lower 95%'].append(edge_rmse[1])
            group_results['Edgewise RMSE, upper 95%'].append(edge_rmse[2])
        group_results['R-squared'].append(r2[0])
        group_results['R-squared, lower 95%'].append(r2[1])
        group_results['R-squared, upper 95%'].append(r2[2])
        group_results['Kendall tau'].append(tau[0])
        group_results['Kendall tau, lower 95%'].append(tau[1])
        group_results['Kendall tau, upper 95%'].append(tau[2])

    df = pd.DataFrame(group_results)
    df.to_csv('group_summaries.csv', index=False, float_format='%.2f')


if __name__== '__main__':
    main()
