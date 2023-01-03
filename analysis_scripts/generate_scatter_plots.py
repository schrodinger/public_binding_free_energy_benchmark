import numpy as np
import pandas as pd
import argparse

import helper_functions as hf
import matplotlib.pylab as plt

# Scatter plot parameters
MARKEREDGEWIDTH = 2
MARKERSIZE = 7
CAPSIZE=5
FMT = 'o'
ZORDER = 2


def gen_subplot_scatter(ax, dgx, dgy, dgx_err=None, dgy_err=None, scatterplot=True):
    """
    Generate a scatter plot of binding free energies that shows some statistics. This is a helper function that is
    written for only one particular plot.

    Parameters
    ----------
    ax: matplotlib.axes._subplots.AxesSubplot
        The subplot where the axis will go.
    dgx: list-like
        The binding free energies for the x-axis.
    dgy: list-like
        The binding free energies for the y-axis
    dgx_error: list-like
        The uncertainty on the x-axis binding free energies values. Assumed to represent standard errors.
    dgy_error: list-like
        The uncertainty on the y-axis binding free energies values. Assumed to represent standard errors.
    scatterplot: bool
        Whether or not to add the scatter plot.

    """
    # Scatter plot with a shaded areas that denote 1 kcal/mol and 2 kcal/mol
    nudge = 0.2
    data = np.hstack((dgx, dgy))
    ymin = xmin = np.min(data) - nudge
    ymax = xmax = np.max(data) + nudge

    xy_line = np.linspace(ymin, ymax)

    # Plot the shaded regions showing 1 and 2 kcal/mol error regions
    ax.fill_between(xy_line, xy_line - 2.0, xy_line + 2.0, color='grey', alpha=0.2, lw=0)
    ax.fill_between(xy_line, xy_line - 1.0, xy_line + 1.0, color='grey', alpha=0.4, lw=0)
    # Plot the x=y line.
    ax.plot(xy_line, xy_line, color='black', alpha=0.8)

    # The scatter plot
    if scatterplot:
        ax.errorbar(dgx, dgy, xerr=dgx_err, yerr=dgy_err,
                    markersize=MARKERSIZE,
                    fmt=FMT,
                    zorder=ZORDER,
                    capsize=CAPSIZE,
                    markeredgewidth=MARKEREDGEWIDTH)

    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))

    # Add correlation and pairwise error information.
    diffs = hf.get_pairwise_diffs(dgx, dgy)
    corr_stats = hf.get_absolute_stats(dgx, dgy)
    boostrap_mues, bootstrap_rmsds = hf.bootstrap_pairwise_error(dgx, dgy, nboots=5000)
    brmds, bmues, br2, bootstrap_taus = hf.bootstrap_absolute_stats(dgx, dgy, nboots=5000)

    xpos = xmin + (xmax - xmin) * 0.05
    ypos = ymin + (ymax - ymin) * 0.93
    s = r"$\tau_{{\Delta G}}$ = {0:.2f} [{1:.2f}, {2:.2f}]".format(corr_stats[3],
                                                                   np.percentile(bootstrap_taus, 2.5),
                                                                   np.percentile(bootstrap_taus, 97.5))
    ax.text(x=xpos, y=ypos, s=s, fontsize=15)
    ypos = ymin + (ymax - ymin) * 0.85
    rmsd = np.sqrt(np.mean(diffs ** 2))
    s = r'$RMSE_{{\Delta\Delta G}}$ = {0:.1f} [{1:.1f}, {2:.1f}] kcal/mol'.format(rmsd,
                                                                                  np.percentile(bootstrap_rmsds, 2.5),
                                                                                  np.percentile(bootstrap_rmsds, 97.5))
    ax.text(x=xpos, y=ypos, s=s, fontsize=15)

def main(argv=None):
    usage = """
        As input, this script requires the directory that contains the CSV files of the individual experimental binding
        free energy comparisons and the FEP+ data. The data sets that this script expects are hard-coded.
        
        To generate the plot used in the publication, the following command was used: 
        
            > python generate_scatter_plots.py scatterplot_data -o scatterplot.png
        """
    description = """
        Scatter plots showing the range of agreement of dGs between different experimental assays (top row) and 
        the agreement between FEP+ predictions and experiment (bottom row). The leftmost column shows examples where the 
        pairwise RMSE of relative binding free energies was much better than average, the middle column shows examples 
        where the RMSE was close to the average, and the rightmost column shows examples where the RMSE was worse than 
        average. The top left of each plot shows Kendall tau and pairwise RMSE for each data set. Points in the dark 
        gray area are measurements or predictions that are within 1 kcal/mol of each other, and points in the light gray 
        area agree within 2 kcal/mol. The top left plot shows that isothermal titration calorimetry (ITC) and 
        fluorescence polarization (FP) binding free energy measurements of galectin ligands are offset by around 
        1 kcal/mol - this offset is irrelevant for rank ordering and does not affect the correlation or the pairwise 
        RMSE metric. The offset of the FEP+ dG predictions was determined by ensuring the mean of the 
        dGs was equal to the experimental mean.
        """
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument(
        'dirname',
        type=str,
        help="The directory that contains the experimental and FEP data that will be plotted.")
    parser.add_argument(
        '-o',
        '--outname',
        type=str,
        help="The name of the scatter plot (a png file) that is produced.")

    args = parser.parse_args(argv)

    # Load the experimental comparison data
    df_exp_galectin = pd.read_csv(f'{args.dirname}/peterson2018_itc_fp_exp_unc.csv')
    df_exp_galectin['dG FP error (kcal/mol)'].loc[pd.isna(df_exp_galectin['dG FP error (kcal/mol)'])] = 0.0
    df_exp_fak = pd.read_csv(f'{args.dirname}/mason2012_spr_ic50.csv')
    df_exp_dppiv = pd.read_csv(f'{args.dirname}/schnapp2016_dppiv_itc_spr_exp_unc.csv')

    # Load the FEP+ predictions
    df_fep_galectin = pd.read_csv(f'{args.dirname}/galectin3_extra_out.csv')
    df_fep_galectin['Exp. dG error (kcal/mol)'].loc[pd.isna(df_fep_galectin['Exp. dG error (kcal/mol)'])] = 0.0
    df_fep_mcl1 = pd.read_csv(f'{args.dirname}/mcl1_extra_flips_bmcorr_out.csv')
    df_fep_fxa = pd.read_csv(f'{args.dirname}/fxa_yoshikawa_set_out.csv')

    # Plotting everything
    fig, axis = plt.subplots(2, 3, figsize=(15, 10))

    ### The experimental comparison boxplots
    print('Galectin')
    print('#########')
    ax = axis[0, 0]

    gen_subplot_scatter(ax, df_exp_galectin['dG ITC (kcal/mol)'], df_exp_galectin['dG FP (kcal/mol)'],
                        scatterplot=False)
    ax.errorbar(df_exp_galectin['dG ITC (kcal/mol)'],
                df_exp_galectin['dG FP (kcal/mol)'],
                xerr=[df_exp_galectin['dG ITC (kcal/mol)'] - df_exp_galectin['dG ITC lower 68%'],
                      df_exp_galectin['dG ITC upper 68%'] - df_exp_galectin['dG ITC (kcal/mol)']],
                yerr=df_exp_galectin['dG FP error (kcal/mol)'],
                markersize=MARKERSIZE,
                fmt=FMT,
                zorder=ZORDER,
                capsize=CAPSIZE,
                markeredgewidth=MARKEREDGEWIDTH)
    ax.set_title('Experimental galectin comparison', fontsize=16)
    ax.set_xlabel(r'$\Delta G$ from ITC (kcal/mol)', fontsize=14)
    ax.set_ylabel(r'$\Delta G$ from biophysical FP (kcal/mol)', fontsize=14)
    print()

    print('FAK')
    print('#########')
    ax = axis[0, 1]
    gen_subplot_scatter(ax, df_exp_fak['SPR - Kd (kcal/mol)'], df_exp_fak['Functional - IC50 (kcal/mol)'])
    ax.set_title('Experimental FAK comparison', fontsize=16)
    ax.set_xlabel(r'$\Delta G$ from SPR (kcal/mol)', fontsize=14)
    ax.set_ylabel(r'$\Delta G$ from biochemical IC50 (kcal/mol)', fontsize=14)
    print()

    print('DPPIV')
    print('#########')
    ax = axis[0, 2]
    gen_subplot_scatter(ax,
                        df_exp_dppiv['dG - ITC (kcal/mol)'],
                        df_exp_dppiv['dG SPR (kcal/mol)'],
                        dgx_err=df_exp_dppiv['dG error - ITC (kcal/mol)'],
                        dgy_err=df_exp_dppiv['dG error SPR (kcal/mol)'])
    ax.set_title('Experimental DPPIV comparison', fontsize=16)
    ax.set_xlabel(r'$\Delta G$ from ITC (kcal/mol)', fontsize=14)
    ax.set_ylabel(r'$\Delta G$ from SPR (kcal/mol)', fontsize=14)
    print()

    print('Galectin')
    print('#########')
    ax = axis[1, 0]
    gen_subplot_scatter(ax,
                        df_fep_galectin['Exp. dG (kcal/mol)'],
                        df_fep_galectin['Pred. dG (kcal/mol)'],
                        dgx_err=df_fep_galectin['Exp. dG error (kcal/mol)'],
                        dgy_err=df_fep_galectin['Pred. dG std. error (kcal/mol)'])
    ax.set_title('FEP+ galectin comparison', fontsize=16)
    ax.set_xlabel(r'$\Delta G$ from biophysical FP (kcal/mol)', fontsize=14)
    ax.set_ylabel(r'$\Delta G$ predicted by FEP+ (kcal/mol)', fontsize=14)
    print()

    print('MCL1')
    print('#########')
    ax = axis[1, 1]
    gen_subplot_scatter(ax,
                        df_fep_mcl1['Exp. dG (kcal/mol)'],
                        df_fep_mcl1['Pred. dG (kcal/mol)'],
                        dgx_err=df_fep_mcl1['Exp. dG error (kcal/mol)'],
                        dgy_err=df_fep_mcl1['Pred. dG std. error (kcal/mol)'])
    ax.set_title('FEP+ MCL1 comparison', fontsize=16)
    ax.set_xlabel('$\Delta G$ from biochemical inhibition (kcal/mol)', fontsize=14)
    ax.set_ylabel(r'$\Delta G$ predicted by FEP+ (kcal/mol)', fontsize=14)
    print()

    print('FXa')
    print('#########')
    ax = axis[1, 2]
    gen_subplot_scatter(ax,
                        df_fep_fxa['Exp. dG (kcal/mol)'],
                        df_fep_fxa['Pred. dG (kcal/mol)'],
                        dgy_err=df_fep_fxa['Pred. dG std. error (kcal/mol)'])
    ax.set_title('FEP+ FXa comparison', fontsize=16)
    ax.set_xlabel(r'$\Delta G$ from biochemical inhibition (kcal/mol)', fontsize=14)
    ax.set_ylabel(r'$\Delta G$ predicted by FEP+ (kcal/mol)', fontsize=14)
    print()

    plt.tight_layout()
    plt.savefig(args.outname, dpi=200)

if __name__== '__main__':
    main()