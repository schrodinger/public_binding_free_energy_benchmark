import numpy as np
from scipy import stats
import pandas as pd

from schrodinger.application.scisol.packages.fep.graph import Graph
from schrodinger.application.scisol.packages.fep import fep_stats


import helper_functions as hf

def read_exp_csv(filename):
    """
    Read a CSV file from the experimental survey. The CSV file is expected to have 2 columns corresponding to different
    assays. Each row is contains two different binding affinity measurements (expected to be in kcal/mol) for the same
    ligand.

    Parameters
    ----------
    filename: str
        The name of the filename that contains the experimental comparison csv.

    Returns
    -------
    dgs1: numpy.ndarray
        The binding free energy for a series of ligands measured using the first assay.
    dgs2: numpy.ndarray
        The binding free energy for a series of ligands measured using the second assay.
    """
    df = pd.read_csv(filename)
    dgs = [df[d] for d in df]
    if len(dgs) != 2:
        raise Exception(f'Only 2 columns are expected in CSV file. File {filename} contains {len(dgs)} columns.')

    return np.array(dgs[0]), np.array(dgs[1])


def weighted_mean(num_per_set, value_per_set):
    """
    Calculate the weighted mean for a value in the experimental survey or FEP benchmark.

    Parameters
    ----------
    num_per_set: np.array
        The weighting attached to a value
    value_per_set: np.array
        The array that will have its weighted mean calculated

    Returns
    -------
    the weighted mean: float
    """
    return np.sum(num_per_set*value_per_set) / np.sum(num_per_set)


def weighted_rmsd(num_per_set, rmsd_per_set):
    """
    Calculate the overall weighted root-mean-square of a set of root-mean-squares (RMS). The calculation first
    calculates the weighted mean-square before returning the square root.

    Parameters
    ----------
    num_per_set: np.array
        The weighting applied to each RMS value.
    rmsd_per_set: np.array
        The RMS values that will be aggregated.

    Returns
    -------
    the weighted RMS.
    """
    return np.sqrt(np.sum(num_per_set * rmsd_per_set**2) / np.sum(num_per_set))


def get_bootstrap_weighted_value(num_set, value_set, nboots=10000):
    """
    Return the bootstrap mean of an array along with uncertainty.

    Parameters
    ----------
    num_set: np.array
        The weighting applied to each value in the value_set
    value_set: np.array
        The array that will have its bootstrap mean calculated.
    nboots: int
        The number of boostrap samples to take.

    Returns
    -------
    the boostrap mean: float
    the bootstrap standard error: float
    the bootstrap 2.5% confidence limit: float
    the boostrap 97.5% confidence limit: float
    """
    value_samples = np.zeros(nboots)
    for b in range(nboots):
        inds = np.random.choice(len(value_set), len(value_set))
        value_samples[b] = weighted_mean(num_set[inds], value_set[inds])

    return value_samples.mean(), value_samples.std(), np.percentile(value_samples, 2.5), np.percentile(value_samples, 97.5)


def get_bootstrap_weighted_rmsd(num_set, rmsd_set, nboots=10000):
    """
    Calculate the boostrap estimate of the overall weighted root-mean-square of a set of root-mean-squares (RMS)

    Parameters
    ----------
    num_set: np.array
        The weighting applied to each value in the rmsd_set
    rmsd_set: np.array
        The array that contains root-mean-squares that will have their boostrap mean calculated.
    nboots: int
        The number of boostrap samples to take.

    Returns
    -------
    the boostrap mean: float
    the bootstrap standard error: float
    the bootstrap 2.5% confidence limit: float
    the boostrap 97.5% confidence limit: float
    """
    rmsd_samples = np.zeros(nboots)
    for b in range(nboots):
        inds = np.random.choice(len(rmsd_set), len(rmsd_set))
        rmsd_samples[b] = weighted_rmsd(num_set[inds], rmsd_set[inds])

    return rmsd_samples.mean(), rmsd_samples.std(), np.percentile(rmsd_samples, 2.5), np.percentile(rmsd_samples, 97.5)


def calculate_tau_from_fep(graph):
    """
    Calculate Kendall's tau from an FEP graph.

    Parameters
    -----------
    graph: schrodinger.application.scisol.packages.fep.graph.Graph
        The output FEP+ file.

    Returns
    -------
    tau: float
        Kendall's rank correlation coefficient.
    """
    # Extract dG values where both exp and pred dG are present
    exp_dgs = []
    pred_dgs = []
    excluded = 0
    nodes_with_exp = []
    for n in graph.nodes_iter():
        if n.exp_dg is None or n.pred_dg is None:
            continue
        if n.is_ccc_excluded:
            excluded += 1
            continue
        exp_dgs.append(n.exp_dg.val)
        pred_dgs.append(n.pred_dg.val)
        nodes_with_exp.append(n)

    result = stats.kendalltau(exp_dgs, pred_dgs)
    return result.correlation


def collect_fep_pairwise_errors(files):
    """
    Get all the pairwise differences between the FEP+ prediction and the experimental reference values.

    Parameters
    ----------
    files: list-like
        The paths to all the output fmp files.

    Returns
    -------
    results: list
        Every single pairwise difference of the predicted ddGs to experimental ddGs (in kcal/mol).
    """
    pairwise_diffs = []
    for name in files:
        g = Graph.deserialize(name)
        exp_dgs = []
        pred_dgs = []
        for n in g.nodes_iter():
            if n.exp_dg is None or n.pred_dg is None or n.is_ccc_excluded:
                continue
            exp_dgs.append(n.exp_dg.val)
            pred_dgs.append(n.pred_dg.val)
        pairwise_diffs.extend(hf.get_pairwise_diffs(np.array(exp_dgs), np.array(pred_dgs), verbose=False))

    return pairwise_diffs


def parse_fep_data(files):
    """
    Collect the FEP errors from a list of FEP+ fmp files

    TODO: document format of output.

    Parameters
    ----------
    files: list-like
        The paths to all the output fmp files.

    Returns
    -------
    results: dict
        A dictionary containing all numpy arrays of each error metric.
    """
    results = {'entries':[], 'number of compounds':[], 'number of edges':[], 'Pairwise RMSE':[], 'Pairwise MUE':[],
               'Edgewise RMSE':[], 'Edgewise MUE':[], 'R-squared':[], 'Kendall tau':[]}

    pairwise_diffs = []
    for name in files:
        results['entries'].append(name.split('/')[-1].split('.')[0])
        g = Graph.deserialize(name)

        # Collect the aggregate stats
        r = fep_stats.calculate(g)
        results['number of compounds'].append(r['Total compounds'])
        results['number of edges'].append(g.number_of_edges())
        results['Edgewise RMSE'].append(r['RMSE Edgewise'].val)
        results['Edgewise MUE'].append(r['MUE Edgewise'].val)
        results['Pairwise RMSE'].append(r['RMSE Pairwise'].val)
        results['Pairwise MUE'].append(r['MUE Pairwise'].val)
        results['R-squared'].append(r['R^2'])
        results['Kendall tau'].append(calculate_tau_from_fep(g))

        # Collect the pairwise errors
        exp_dgs = []
        pred_dgs = []
        for n in g.nodes_iter():
            if n.exp_dg is None or n.pred_dg is None or n.is_ccc_excluded:
                continue
            exp_dgs.append(n.exp_dg.val)
            pred_dgs.append(n.pred_dg.val)
        pairwise_diffs.extend(hf.get_pairwise_diffs(np.array(exp_dgs), np.array(pred_dgs), verbose=False))

    for key in results:
        results[key] = np.array(results[key])

    return results, np.array(pairwise_diffs)


def parse_fep_data_from_csv(files):
    """
    Collect the FEP errors from a list of FEP+ fmp files

    Parameters
    ----------
    files: list-like
        The paths to all the output csv files.

    Returns
    -------
    results: dict
        A dictionary containing all numpy arrays of each error metric.
    """
    results = {'entries': [],
               'number of compounds': [],
               'Pairwise RMSE': [],
               'Pairwise MUE': [],
               'R-squared': [],
               'Kendall tau': []}

    pairwise_diffs = []
    for name in files:
        results['entries'].append(name.split('/')[-1].split('.')[0])

        df = pd.read_csv(name)
        rmsd, mue, r2, tau = hf.get_absolute_stats(df['Pred. dG (kcal/mol)'], df['Exp. dG (kcal/mol)'], verbose=False)

        diffs = hf.get_pairwise_diffs(df['Pred. dG (kcal/mol)'], df['Exp. dG (kcal/mol)'], verbose=False)
        pairwise_diffs.extend(diffs)

        # Collect the aggregate stats
        results['number of compounds'].append(len(df))
        results['Pairwise RMSE'].append(np.sqrt(np.mean(diffs ** 2)))
        results['Pairwise MUE'].append(np.mean(np.abs(diffs)))
        results['R-squared'].append(r2)
        results['Kendall tau'].append(tau)

    for key in results:
        results[key] = np.array(results[key])

    return results, np.array(pairwise_diffs)


def summarize_fep_error(results, verbose=True):
    """
    Calculate the weighted errors and correlation statistics for the FEP benchmark.

    Parameters
    ----------
    results: dict of numpy.array
        Each key is a error metric that points to an array of these metrics for each map in the benchmark set.
    verbose: bool
        Whether to print out the summary statistics.

    Returns
    -------
    pairwise_rmse: tuple
        The weighted pairwise RMSE along with bootstrapped lower and upper 95% confidence intervals.
    pairwise_mue: tuple
        The weighted pairwise MUE
    edgewise_rmse: tuple
        If edgewise data is present, the edgewise RMSE with bootstrapped lower and upper 95% confidence intervals.
        Skipped if there is not edge data in results.
    edgewise_mue: tuple
        If edgewise data is present, the edgewise MUE with bootstrapped lower and upper 95% confidence intervals.
        Skipped if there is not edge data in results.
    rsquared: tuple
        The weighted R-squared of the absolute binding free energies, along with bootstrapped lower and upper 95%
        confidence intervals.
    kendall_tau: tuple
        The weighted Kendall tau of the absolute binding free energies, along with bootstrapped lower and upper 95%
        confidence intervals.
    """
    # Allowing for cases when edge data is not available:
    if 'number of edges' in results and 'Edgewise RMSE' in results and 'Edgewise MUE' in results:
        edge_data = True
    else:
        edge_data = False

    num_comps = results['number of compounds']
    if edge_data:
        num_edges = results['number of edges']

    # Pairwise RMSEs
    m, std, pair_l, pair_u = get_bootstrap_weighted_rmsd(num_comps, results['Pairwise RMSE'])
    pair_m = weighted_rmsd(num_comps, results['Pairwise RMSE'])

    # Pairwise MUEs
    m, std, pair_mue_l, pair_mue_u = get_bootstrap_weighted_value(num_comps, results['Pairwise MUE'])
    pair_mue_m = weighted_mean(num_comps, results['Pairwise MUE'])

    # Correlation and rank
    m, std, r2_l, r2_u = get_bootstrap_weighted_value(num_comps, results['R-squared'])
    r2_m = weighted_mean(num_comps, results['R-squared'])

    m, std, tau_l, tau_u = get_bootstrap_weighted_value(num_comps, results['Kendall tau'])
    tau_m = weighted_mean(num_comps, results['Kendall tau'])

    if verbose:
        print(f'Pair RMSE = {pair_m:.2f} [{pair_l:.2f}, {pair_u:.2f}] kcal/mol')
        print(f'Pair MUE  = {pair_mue_m:.2f} [{pair_mue_l:.2f}, {pair_mue_u:.2f}] kcal/mol')

    if edge_data:
        m, std, edge_l, edge_u = get_bootstrap_weighted_rmsd(num_edges, results['Edgewise RMSE'])
        edge_m = weighted_rmsd(num_edges, results['Edgewise RMSE'])

        m, std, edge_mue_l, edge_mue_u = get_bootstrap_weighted_value(num_edges, results['Edgewise MUE'])
        edge_mue_m = weighted_mean(num_edges, results['Edgewise MUE'])

        if verbose:
            print(f'Edge RMSE = {edge_m:.2f} [{edge_l:.2f}, {edge_u:.2f}] kcal/mol')
            print(f'Edge MUE  = {edge_mue_m:.2f} [{edge_mue_l:.2f}, {edge_mue_u:.2f}] kcal/mol')

    if verbose:
        print(f'R2        = {r2_m:.2f} [{r2_l:.2f}, {r2_u:.2f}]')
        print(f'Tau       = {tau_m:.2f} [{tau_l:.2f}, {tau_u:.2f}]')

    if edge_data:
        return (pair_m, pair_l, pair_u), \
               (pair_mue_m, pair_mue_l, pair_mue_u), \
               (edge_m, edge_l, edge_u), \
               (edge_mue_m, edge_mue_l, edge_mue_u), \
               (r2_m, r2_l, r2_u), \
               (tau_m, tau_l, tau_u)

    else:
        return (pair_m, pair_l, pair_u), \
               (pair_mue_m, pair_mue_l, pair_mue_u), \
               (r2_m, r2_l, r2_u), \
               (tau_m, tau_l, tau_u)


def collect_experimental_pairwise_diffs(files, notlist=()):
    """
    Get all the pairwise differences from each assay comparison from the experimental survey.

    Parameters
    ----------
    files: list-like
        The paths to all the files that will have their error aggregated.
    notlist: list-like
        The names of the csv files that will be omitted from the analysis.

    Returns
    -------
    results: list
        Every single pairwise difference of the measured ddGs (in kcal/mol).
    """
    pairwise_diffs = []
    for f in files:
        entry = f.split('/')[-1].split('.')[0]
        if entry not in notlist:
            dg1, dg2 = read_exp_csv(f)
            pairwise_diffs.extend(hf.get_pairwise_diffs(dg1, dg2, verbose=False))

    return pairwise_diffs


def parse_experimental_data(files, notlist=()):
    """
    Collect the of experimental errors from a set of csv files that contain all the comparative data.

    Parameters
    ----------
    files: list-like
        The paths to all the files that will have their error aggregated.
    notlist: list-like
        The names of the csv files that will be omitted from the analysis.

    Returns
    -------
    results: dict
        Contains the arrays of the data RMSEs, MUEs, and correlation statistics.
    """
    results = {'entries':[], 'number':[], 'Pairwise RMSE':[], 'Pairwise MUE':[], 'Absolute RMSE':[], 'Absolute MUE':[],
               'R-squared':[], 'Kendall tau':[]}

    pairwise_diffs = []
    for f in files:
        entry = f.split('/')[-1].split('.')[0]
        if entry in notlist:
            #print('Skipping {}'.format(entry))
            pass
        else:
            # Collect the aggregate stats
            results['entries'].append(entry)
            # Correlation stats
            dg1, dg2 = read_exp_csv(f)
            rmsd, mae, c1, c2  = hf.get_absolute_stats(dg1, dg2, verbose=False)
            results['R-squared'].append(c1)
            results['Kendall tau'].append(c2)
            results['number'].append(len(dg1))
            # Abosulute errors
            results['Absolute RMSE'].append(np.sqrt(np.mean((dg1 - dg2)**2)))
            results['Absolute MUE'].append(np.mean(np.abs(dg1 - dg2)))
            # Pairwise errors
            diffs = hf.get_pairwise_diffs(dg1, dg2, verbose=False)
            results['Pairwise RMSE'].append(np.sqrt(np.mean(diffs**2)))
            results['Pairwise MUE'].append(np.mean(np.abs(diffs)))

            # Store the pairwise differences
            pairwise_diffs.extend(diffs)

    for key in results:
        results[key] = np.array(results[key])

    return results, np.array(pairwise_diffs)


def summarize_experimental_error(files, notlist=()):
    """
    Calculate and print the summary statistics from the set of csv files that contail the experimental comparison data.

    Parameters
    ----------
   files: list-like
        The paths to all the files that will have their error aggregated.
    notlist: list-like
        The names of the csv files that will be omitted from the analysis.
    """

    results, diffs = parse_experimental_data(files, notlist)

    print('Total number of comparison data points (including repeated ligands) =', np.sum(results['number']))
    print()

    for key in results:
        if key == 'number' or key == 'entries':
            pass
        elif 'RMSE' in key:
            boot_mean, boot_std, boot_lower, boot_upper = get_bootstrap_weighted_rmsd(results['number'], results[key])
            print('Weighted {} = {:.2} kcal/mol'.format(key, weighted_rmsd(results['number'], results[key])))
            print(
                'Weighted bootstap {} = {:.2} [{:.2f}, {:.2f}] kcal/mol'.format(key, boot_mean, boot_lower, boot_upper))
            print()
        else:
            boot_mean, boot_std, boot_lower, boot_upper = get_bootstrap_weighted_value(results['number'], results[key])
            print('Weighted {} = {:.2} kcal/mol'.format(key, weighted_mean(results['number'], results[key])))
            print(
                'Weighted bootstap {} = {:.2} [{:.2f}, {:.2f}] kcal/mol'.format(key, boot_mean, boot_lower, boot_upper))
            print()


def error_diff_stats(diffs):
    """
    Print out the percentage of differences less than 1 kcal/mol and the fraction of differences greater than 2 kcal/mol.
    """
    abs_diffs = np.abs(diffs)
    total = len(diffs)
    frac_less_1 = np.sum(abs_diffs < 1) / total
    frac_more_2 = np.sum(abs_diffs > 2) / total
    print(f'There are a total of {total} differences')
    print(f'{100 * frac_less_1:.1f}% of the differences are less than 1 kcal/mol' )
    print(f'{100 * frac_more_2:.1f}% of the differences are greater than 2 kcal/mol' )


def print_latex_table(fmpnames, out2pdb=None, out2protein=None):
    """
    Print out the errors for a collection of FEP+ maps in a latex formatted table.

    fmpnames: list of str
        The paths to a every FEP+ output file you want to put in a latex table.
    out2pdb: dict
        A dictionary that links the output filename to a PDB
    out2protein: dict
        A dictionary that links the output filename to a protein name.
    """
    num_compounds = []
    num_edges = []
    pairwise_rmse = []
    edge_rmse = []
    r2 = []
    lines = []
    for name in fmpnames:
        g = Graph.deserialize(name)
        r = fep_stats.calculate(g)
        num_compounds.append(r['Total compounds'])
        num_edges.append(g.number_of_edges())
        edge_rmse.append(r['RMSE Edgewise'].val)
        pairwise_rmse.append(r['RMSE Pairwise'].val)
        r2.append(r['R^2'])

        outname = name.split('/')[-1].split('.')[0]
        if out2pdb is not None:
            pdb = out2pdb[outname]
        else:
            pdb = 'XXXX'
        if out2protein is not None:
            protein = out2protein[outname]
        else:
            protein = outname

        line = r'    {:35} & {} & {} & {} & {:.2f} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} \\'.format(protein,
                                                                                           pdb,
                                                                                           r['Total compounds'],
                                                                                           g.number_of_edges(),
                                                                                           r['R^2'],
                                                                                           r['RMSE Edgewise'].val,
                                                                                           r['RMSE Edgewise'].unc,
                                                                                           r['RMSE Pairwise'].val,
                                                                                           r['RMSE Pairwise'].unc)
        lines.append(line)

    num_compounds = np.array(num_compounds)
    num_edges = np.array(num_edges)
    r2 = np.array(r2)
    edge_rmse = np.array(edge_rmse)
    pairwise_rmse = np.array(pairwise_rmse)

    m, std, pair_l, pair_u  = get_bootstrap_weighted_rmsd(num_compounds, pairwise_rmse)
    pair_m = weighted_rmsd(num_compounds, pairwise_rmse)

    m, std, edge_l, edge_u  = get_bootstrap_weighted_rmsd(num_edges, edge_rmse)
    edge_m = weighted_rmsd(num_edges, edge_rmse)

    m, std, r2_l, r2_u  = get_bootstrap_weighted_value(num_compounds, r2)
    r2_m = weighted_mean(num_compounds, r2)

    last_line = r'\hline & {:35} & {} & {} & {:.2f} [{:.2f}, {:.2f}] & {:.2f} [{:.2f}, {:.2f}] & {:.2f} [{:.2f}, {:.2f}] \\'.format('Total/weighted mean',
                                                                                                 np.sum(num_compounds),
                                                                                                 np.sum(num_edges),
                                                                                                 r2_m,
                                                                                                 r2_l,
                                                                                                 r2_u,
                                                                                                 edge_m,
                                                                                                 edge_l,
                                                                                                 edge_u,
                                                                                                 pair_m,
                                                                                                 pair_l,
                                                                                                 pair_u)
    header = r"""\begin{table}[!ht]
    \centering
    \small
    \caption{}
    \begin{tabular}{c|c|c|c|c|c|c}
    \toprule
    System & PDB & No. compounds & No. edges & R$^2$ & Edgewise RMSE & Pairwise RMSE \\\hline"""
    footer = r"""    \bottomrule
        \end{tabular}
        \label{tab:}
\end{table}
    """

    print(header)
    for l in lines:
        print(l)
    print(last_line)
    print(footer)


