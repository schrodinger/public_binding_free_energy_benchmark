import numpy as np
import matplotlib.pylab as plt
from itertools import combinations
from scipy import stats

def pretty_scatter(xdata, ydata, figsize=(5,5), nudge=0.2):

    fig, ax = plt.subplots(1, 1,figsize=figsize)
    data = np.hstack((xdata, ydata))
    ymin = xmin = np.min(data) - nudge
    ymax = xmax = np.max(data) + nudge

    xy_line = np.linspace(ymin, ymax)

    ax.fill_between(xy_line, xy_line - 2.0, xy_line + 2.0, color='grey', alpha=0.2, lw=0)
    ax.fill_between(xy_line, xy_line - 1.0, xy_line + 1.0, color='grey', alpha=0.4, lw=0)
    ax.plot(xy_line, xy_line, color='black', alpha=0.8)

    ax.scatter(xdata, ydata, s=55, alpha=1, zorder=2)

    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))
    
    ax.tick_params(axis='x',labelsize=12)
    ax.tick_params(axis='y',labelsize=12)
    
    return fig, ax


def get_pairwise_diffs(data1, data2, inds=None, verbose=True):
    """
    Calculate the pairwise differences between two data series.

    Parameters
    ----------
    data1: list-like of floats
        The data series for one variable.
    data2: List-like of floats
        The data series for the second variable.
    inds: list-like of ints
        The indices of the elements in the data series that will be used to calculate the pairwise error. By default,
        all indices are used.
    verbose: bool
        Whether to print out a summary of the pairwise errors.

    Returns
    -------
    diffs: numpy.ndarray
        The array of all pairwise differences. If the data series have a length of N, then the length of diffs will be
        N * (N - 1) / 2.
    """
    if len(data1) != len(data2):
        raise Exception('Length of inputs do not match')
    
    if inds is None:
        inds = np.arange(len(data1))

    pairwise_diffs1 = np.array([data1[i] - data1[j] for i, j in combinations(inds, 2)])
    pairwise_diffs2 = np.array([data2[i] - data2[j] for i, j in combinations(inds, 2)])
    
    diffs = pairwise_diffs2 - pairwise_diffs1

    if verbose:
        rmsd = np.sqrt(np.mean(diffs ** 2))
        print(len(inds), 'total compounds')
        print(len(pairwise_diffs1), 'pairwise differences to compare')
        print('Pairwise MUE = {:.2f} kcal/mol'.format(np.mean(np.abs(diffs))))
        print('Pairwise RMSD = {:.2f} kcal/mol'.format(rmsd))

    return diffs


def bootstrap_pairwise_error(data1, data2, nboots=5000):
    """
    Generate bootstrap samples of the pairwise errors of 2 data series.

    Parameters
    ----------
    data1: list-like of floats
        The data series for one variable.
    data2: List-like of floats
        The data series for the second variable.
    nboots: int
        The number of bootstrap samples to generate.

    Returns
    -------
    boot_mues: numpy.ndarray
        A 1D array of bootstrap samples of the mean unsigned error (mean absolute error) of the pairwise error of the
        data series.
    boot_rmsds: numpy.ndarray
        A 1D array of bootstrap samples of the root-mean-square difference of the pairwise error of the
        data series.
    """
    boot_mues = np.zeros(nboots)
    boot_rmsds = np.zeros(nboots)
    for i in range(nboots):
        inds = np.random.choice(len(data1), len(data1))
        diffs = get_pairwise_diffs(data1, data2, inds, verbose=False)
        boot_mues[i] = np.mean(np.abs(diffs))
        boot_rmsds[i] = np.sqrt(np.mean(diffs ** 2))

    return boot_mues, boot_rmsds


def get_absolute_stats(dg_set1, dg_set2, verbose=True):
    """
    Get the square of Pearson correlation coefficient and Kendall's tau
    """
    r2 = np.corrcoef(dg_set1, dg_set2)[0,1]**2
    t = stats.kendalltau(dg_set1, dg_set2).correlation
    if verbose:
        print('R-squared = {}'.format(np.round(r2,2)))
        print("Kendaull's tau = {}".format(np.round(t,2)))
    
    # Also through in absolute deviation
    abs_diff = np.abs(dg_set1 - dg_set2)
    rmsd = np.sqrt(np.mean(abs_diff**2))
    if verbose:
        print('Absolute MAE =  {} kcal/mol'.format(np.round(abs_diff.mean(), 2)))
        print('Absolute RMSE = {} kcal/mol'.format(np.round(rmsd,2)))
          
    return rmsd, abs_diff.mean(), r2, t

def bootstrap_absolute_stats(data1, data2, nboots=5000):
    """
    Generate bootstrap samples of the absolute mean error, root-mean-square error, and correlation stats of two data
    series.

    Parameters
    ----------
    data1: list-like of floats
        The data series for one variable.
    data2: List-like of floats
        The data series for the second variable.
    nboots: int
        The number of bootstrap samples to generate.
    """
    boot_mues = np.zeros(nboots)
    boot_rmsds = np.zeros(nboots)
    boot_r2 = np.zeros(nboots)
    boot_t = np.zeros(nboots)
    for i in range(nboots):
        inds = np.random.choice(len(data1), len(data1))
        boot_rmsds[i], boot_mues[i], boot_r2[i], boot_t[i] =  get_absolute_stats(data1[inds], data2[inds], verbose=False)
        #rmsd, mue, r2, t = get_absolute_stats(data1[inds], data2[inds], verbose=False)
        #boot_mues[i] = mue
        #boot_rmsds[i] = rmsd
        #boot_r2[i] = r2

    return boot_mues, boot_rmsds, boot_r2, boot_t

def read_outfiles(txtfilename, prefix=None):
    """
    Append the file paths that are listed in a text file.

    Parameters
    ----------
    txtfilename: str
        The path to a text file that contains a different output FEP+ file in each new line.
    prefix: str
        Additional text to place before the file path, such as the upper directory.
    Returns
    -------
    ext_maps: list of str
        The paths to the FEP maps in the text file.
    """
    if prefix is None:
        prefix = ''
    with open(txtfilename, "r") as file:
        lines = file.readlines()
        fmps = [prefix + l.strip() for l in lines if l.strip() != '']
    return fmps


def print_latex_group_summary_table(csv):
    """
    Print out the FEP+ accuracy statistics for each group in a latex-ready format. This function was used for a table
    in the supporting information.

    Parameters
    ----------
    csv: str
        The CSV file that contains the group aggregate statistics. e.g '../21_4_results/summary_statistics/group_summaries.csv'
    """
    df = pd.read_csv(csv)
    df = df.reindex([4, 2, 12, 5, 6, 9, 11, 13, 7, 1, 0, 3, 10, 8])
    names = ['Public Merck',
             'FEP+ fragments',
             'GPCRs',
             'Janssen BACE1',
             'Bayer macrocycles',
             'MCS docking',
             'FEP+ scaffold-hopping',
             'FEP+ macrocycles',
             'FEP+ buried water',
             'FEP+ charge-change',
             'Miscellaneous',
             'FEP+ R-group',
             'OPLS stress set',
             'OPLS drug discovery']
    df['Nice name'] = names
    print('Group & No. compounds & No. edges & R$^2$ & Edgewise RMSE & Pairwise RMSE \\\\hline')
    for i, r in df.iterrows():
        line = f"{r['Nice name']} & " \
               f"{r['No. compounds']} & " \
               f"{r['No. edges']} & " \
               f"{r['R-squared']} [{r['R-squared, lower 95%']}, {r['R-squared, upper 95%']}] & " \
               f"{r['Edgewise RMSE']} [{r['Edgewise RMSE, lower 95%']}, {r['Edgewise RMSE, upper 95%']}] & " \
               f"{r['Pairwise RMSE']} [{r['Pairwise RMSE, lower 95%']}, {r['Pairwise RMSE, upper 95%']}] \\\ "
        print(line)