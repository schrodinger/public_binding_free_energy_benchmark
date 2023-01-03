import analysis_functions as af
import importlib
import matplotlib.pylab as plt
import numpy as np
from scipy import stats
import argparse
from glob import glob
import os

def _summarize_differences_stats(results, diffs):
    """
    Print a summary of the statistics of the supplied array - which is assumed to be a set of differences. Nothing is
    returned and the results are printed out.

    Parameters
    ----------
    results: dict
        The output from analysis_functions.parse_fep_data, or analysis_functions.parse_fep_data_from_csv or
        analysis_functions.parse_experimental_data which contains the summary statistics from the FEP benchmark or
        experimental survey.
    diffs: numpy.ndarray
        An array pairwise differences between experimental measurements or FEP predictions against experiemntal
        measurements.
    """
    af.error_diff_stats(diffs)
    print()
    print('The Shapiro-Wilk test tests for the null hypothesis that the differences are drawn from a normal distribution:')
    print(stats.shapiro(diffs))
    print()
    print('Median pairwise RMSE from assays: {:.2f} kcal/mol'.format(np.percentile(results['Pairwise RMSE'],50)))
    print()


def main(argv=None):
    usage = """
        As input, this script takes the experimental survey and the FEP benchmark data as input. The FEP data can be in
        CSV format (e.g from ../21_4_results/ligand_predictions) or as FMP format (e.g. from 
        ../21_4_results/processed_output_fmps). The FMP or CSV files must be grouped in subdirectories of the input 
        leading directory. 
        
        The plot that was used in the publication using the output FEP+ FMP files was created using 
        
            > $SCHRODINGER/run python3 generate_boxplots_and_histograms.py --fep_dir ../21_4_results/processed_output_fmps -e fmp --exp_files ../experimental_survey_data/*/*csv -o scatterplot_exp_fep.png
        
        Alternatively, the CSV output files can also be be processed by this script. 
        
            > python generate_boxplots_and_histograms.py --fep_dir ../21_4_results/ligand_predictions/ -e csv --exp_files ../experimental_survey_data/*/*csv -o scatterplot_exp_fep_from_csv.png
        """
    description = """
        The left panel shows boxplots comparing the root-mean-square error (RMSE) between relative binding free energies
        from different experimental assays (left) and the FEP+ predictions against experimental data (right). The size 
        of each data point is proportional to the number of ligands in the series in either an assay comparison or 
        perturbation graph. The two largest data points in the experimental survey are from the COVID moonshot project 
        and project A from table S4. The median RMSE in the experimental survey is 0.85 kcal/mol and the median in the 
        FEP+ benchmark is 1.08 kcal/mol. The right plot shows the all pairwise relative binding free energy differences 
        from the experimental survey and all pairwise FEP+ errors. The error distributions are bell-shaped and can be 
        approximated by t-distributions.
        """
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument(
        '--fep_dir',
        type=str,
        help="The upper directory that contains all of the FMP or CSV outputs files in subdirectories.")
    parser.add_argument(
        '-e',
        '--ext',
        type=str,
        choices=['fmp', 'csv'],
        help="The file extension of the FEP results. Results can be either FMP files or CSVs.")
    parser.add_argument(
        '--exp_files',
        type=str,
        nargs='+',
        help="A list of CSV files that contain experimental data comparisons.")
    parser.add_argument(
        '-o',
        '--outname',
        type=str,
        help="The name of the png file that is produced.")


    args = parser.parse_args(argv)

    if args.ext == 'fmp':
        if importlib.util.find_spec('schrodinger') is None:
            raise Exception('Schrodinger must be installed to use FMP files as input. Use CSV files instead.')

    ##################################
    ### Load the experimental data ###
    ##################################

    # Excluding these experimental comparisons as they are comparing assays that are too similiar.
    notlist = ['bocquet2015_spr_detergent_nanoc9',
               'bocquet2015_spr_detergent_nanohis',
               'bocquet2015_spr_nanoC9_nanohis',
               'jia2006_cot_inhibition',
               'patil2018_ic50_spr',
               'moonshot2020_covid_protease_inhibition',
               'ycas2020_bptf_spr_labeled_spr',
               'murphy2006_binding']

    exp_results, exp_diffs = af.parse_experimental_data(args.exp_files, notlist)

    print('Experimental pairwise error distribution stats:')
    print('----------------------------------------------')
    _summarize_differences_stats(exp_results, exp_diffs)

    ##################################
    ######## Load the FEP data #######
    ##################################

    files = []
    for entry in glob(f'{args.fep_dir}/*'):
        if os.path.isdir(entry):
            files.extend(glob(f'{entry}/*{args.ext}'))

    if len(files) == 0:
        raise Exception(f'No FEP results files found in the subdirectories of {args.fep_dir}. '
                        f'Check the directory name and its contents.')

    extensions = [f.split('.')[-1] for f in files]
    if all([e == 'csv' for e in extensions]):
        fep_results, fep_diffs = af.parse_fep_data_from_csv(files)
    elif all([e == 'fmp' for e in extensions]):
        fep_results, fep_diffs  = af.parse_fep_data(files)

    else:
        raise Exception(f'The supplied files must all be FMP files or CSV files. The following file extensions have '
                        f'been supplied: {set(extensions)}')

    print()
    print('FEP pairwise error distribution stats:')
    print('--------------------------------------')
    _summarize_differences_stats(fep_results, fep_diffs)

    ##################################
    ########## The plotting ##########
    ##################################

    def num2sizes(number_in_set):
        return number_in_set * 5 + 20

    fig, axis = plt.subplots(1, 2, figsize=(13, 6))
    exp_color = 'C4'
    fep_color = 'C2'
    ax = axis[1]
    fep_diffs = np.array(fep_diffs)
    sym_fep_diffs = np.hstack((fep_diffs, -fep_diffs))

    exp_diffs = np.array(exp_diffs)
    sym_exp_diffs = np.hstack((exp_diffs, -exp_diffs))

    plt.ioff()
    fep_heights, fep_bins, patches = plt.hist(sym_fep_diffs, bins=30, density=True)
    exp_heights, exp_bins, patches = plt.hist(sym_exp_diffs, bins=70, density=True)
    # Taking the bin centers, not the edges
    fep_bins = fep_bins[0:-1] + (fep_bins[1] - fep_bins[0]) / 2
    exp_bins = exp_bins[0:-1] + (exp_bins[1] - exp_bins[0]) / 2
    plt.close()

    #################################
    fig, axis = plt.subplots(1, 2, figsize=(13, 6))

    ############ BOXPLOT ############
    ax = axis[0]

    ax.boxplot((exp_results['Pairwise RMSE'], fep_results['Pairwise RMSE']), widths=0.5, showfliers=False,
               patch_artist=True, vert=True,
               boxprops={'facecolor': 'white', 'color': 'k', 'alpha': 1}, medianprops=dict(color=exp_color),
               zorder=1)
    bplot = ax.boxplot((exp_results['Pairwise RMSE'], fep_results['Pairwise RMSE']), widths=0.5, showfliers=False,
                       patch_artist=True, vert=True,
                       boxprops={'facecolor': exp_color, 'color': 'k', 'alpha': 0.5}, medianprops=dict(color=exp_color),
                       zorder=1)

    for patch, color in zip(bplot['boxes'], [exp_color, fep_color]):
        patch.set_facecolor(color)
    for patch, color in zip(bplot['medians'], [exp_color, fep_color]):
        patch.set_color(color)
        patch.set_linewidth(2)

    # The experimental data
    inds = np.flip(np.argsort(exp_results['number']))
    sorted_rmsds = exp_results['Pairwise RMSE'][inds]
    sorted_num = exp_results['number'][inds]
    sizes = num2sizes(sorted_num)

    xpos = np.repeat(1, len(exp_results['Pairwise RMSE']))
    ax.scatter(xpos, sorted_rmsds, s=sizes, alpha=1, color='white', edgecolor='None', linewidths=2, zorder=1)
    ax.scatter(xpos, sorted_rmsds, s=sizes, alpha=0.4, color=exp_color, edgecolor=None, zorder=2)
    ax.scatter(xpos, sorted_rmsds, s=sizes, color='None', edgecolor='black', linewidths=0.5, zorder=3)

    # The FEP+ data
    inds = np.flip(np.argsort(fep_results['number of compounds']))
    sorted_rmsds = fep_results['Pairwise RMSE'][inds]
    sorted_num = fep_results['number of compounds'][inds]

    sizes = num2sizes(sorted_num)

    xpos = np.repeat(2, len(fep_results['Pairwise RMSE']))
    ax.scatter(xpos, sorted_rmsds, s=sizes, alpha=1, color='white', edgecolor='None', linewidths=2, zorder=1)
    ax.scatter(xpos, sorted_rmsds, s=sizes, alpha=0.3, color=fep_color, edgecolor=None, zorder=2)
    ax.scatter(xpos, sorted_rmsds, s=sizes, color='None', edgecolor='black', linewidths=0.4, zorder=3)

    # markerfacecolor="None"

    for i in (0, 0.5, 1., 1.5, 2., 2.5, 3.):
        ax.axhline(i, ls='--', lw=1, color='grey', zorder=0)

    ax.set_ylabel('Pairwise RMSE (kcal/mol)', fontsize=18)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xticks((1, 2))
    ax.set_xticklabels(labels=('Experimental survey', 'FEP+ benchmark'), fontdict={'fontsize': 16})

    ############ Histograms ############
    ax = axis[1]

    ax.grid(ls='--', zorder=0)
    ax.fill_between(fep_bins, 0, fep_heights / fep_heights.max(), color='white', alpha=1, zorder=2)
    ax.fill_between(fep_bins, 0, fep_heights / fep_heights.max(), color='C2', alpha=0.5, label='FEP+', zorder=2)
    ax.plot(fep_bins, fep_heights / fep_heights.max(), color='C2', lw=4, zorder=2)

    ax.fill_between(exp_bins, 0, exp_heights / exp_heights.max(), color='white', alpha=1, zorder=2)
    ax.plot(exp_bins, exp_heights / exp_heights.max(), color='C4', lw=4, zorder=3)
    ax.fill_between(exp_bins, 0, exp_heights / exp_heights.max(), color='C4', alpha=0.5, label='Exp.', zorder=3)

    ax.set_ylim(0, 1.06)
    ax.set_yticklabels(())

    ax.tick_params(axis='y', labelsize=16, color='white')
    ax.tick_params(axis='x', labelsize=14)
    ax.legend(fontsize=16)
    ax.set_xlabel(r'$\Delta\Delta G$ error (kcal/mol)', fontsize=16)
    ax.set_ylabel('Scaled density', fontsize=16)

    plt.tight_layout()
    plt.savefig(args.outname, dpi=200)

if __name__== '__main__':
    main()
