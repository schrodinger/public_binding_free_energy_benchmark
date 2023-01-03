import analysis_functions as af
import argparse
import pandas as pd
from glob import glob
import os


def main(argv=None):
    usage = """
    The script takes output FMP (not CSV files) and prints formatted latex tables that show some accuracy statistics of 
    each system. These latex tables were used in the supporting information of the accompanying manuscript.
 

        > $SCHRODINGER/run python3 print_latex_tables.py ../21_4_results/processed_output_fmps

    Optionally, output metadata can be supplied that provides nicer system names:
        
        > $SCHRODINGER/run python3 print_latex_tables.py ../21_4_results/processed_output_fmps -m ../21_4_results/benchmark_output_metadata.csv
    
    The CSV file with the metadata should be a table with columns 'Input file naming scheme' ,
    'Output file naming scheme', and 'Reference PDB'.
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
        '-m',
        '--metadata',
        type=str,
        help="The CSV file with the output metadata, default=None")
    args = parser.parse_args(argv)

    files = []
    for entry in glob(f'{args.upper_dir}/*'):
        if os.path.isdir(entry):
            files.extend(glob(f'{entry}/*fmp'))

    if  args.metadata is not None:
        df_meta = pd.read_csv(args.metadata)
    else:
        df_meta = None

    for entry in glob(f'{args.upper_dir}/*'):
        if os.path.isdir(entry):
            files = glob(f'{entry}/*fmp')
            if len(files) >= 1:
                group_abbreviation = entry.split('/')[-1]
                # TODO: GPCRs currently requires manually intervention.
                if group_abbreviation == 'gpcrs':
                    continue
                print(group_abbreviation)
                if df_meta is not None:
                    # Use the metadata to create dictionaries
                    out2pdb = {}
                    out2protein = {}
                    df_group = df_meta.loc[df_meta['Group abbreviation'] == group_abbreviation]
                    for i, row in df_group.iterrows():
                        out2pdb[row['Output file naming scheme']] = row['Reference PDB']
                        out2protein[row['Output file naming scheme']] = row['Protein']
                    af.print_latex_table(files, out2pdb=out2pdb, out2protein=out2protein)
                else:
                    af.print_latex_table(files)
                    print()
                    print()


if __name__ == '__main__':
    main()
