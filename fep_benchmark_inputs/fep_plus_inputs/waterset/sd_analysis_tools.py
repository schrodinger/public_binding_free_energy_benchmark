import numpy as np
from copy import deepcopy
from schrodinger.application.scisol.packages.fep import graph
from schrodinger.application.scisol.packages.fep import fep_stats
from schrodinger.application.desmond.measurement import Measurement
from itertools import combinations


from scipy.optimize import minimize


def selection_penaly(pka, ph=7., protonated=True):
    """
    The free energy to select out the either the protonated or deprotonated state of a molecule in solvent.

    Note
    ----
    The way to think about it is via the grand canonical ensemble. These penalties are the ratios of two grand canonical
    partition functions where one state has been removed.

    Parameters
    ----------
    ph = float
        The pH of the solvent.
    pka: float
        The pKa of the protonated form.
    protonated: bool
        Whether to spit out the penalty to select either the protonated form, or deprotonated form.
    """
    top = 10 ** (pka - ph)
    if protonated:
        return -0.596 * np.log(top / (1 + top))
    else:
        return -0.596 * np.log(1 / (1 + top))


def solvent_penalty(pka1, pka2, ph=7):
    """
    Apply the relative free energy penalty to isoloate the neutral forms of two ligands. Used to for the
    pKa correction of the SD map.
    """
    top1 = 10 ** (pka1 - ph)
    top2 = 10 ** (pka2 - ph)
    return -0.596 * np.log((1 + top1) / (1 + top2))


def apply_solvent_penalty(g_pka, pka_dict):
    """
    Calculate the solvent pKa penalty for the ligands in an fmp map.

    Parameters
    ----------
    g_pka: fmp map
        The map whose entries will be changed in place to reflect the pka correction.
    pka_dict: dict
        The dictionary that stores (by ligand name) the macro pKas of the ligands.
    """
    for e in g_pka.edges_iter():
        l1 = e.short_id_title.split('==>')[0].split(':')[1].strip()
        l2 = e.short_id_title.split('==>')[1].split(':')[1].strip()
        corr = solvent_penalty(pka_dict[l1], pka_dict[l2])
        unc = e.complex_dg.unc
        val = e.complex_dg.val
        e.complex_dg = Measurement(val + corr, unc)
    g_pka.calc_cycle_closure()

def pretty_print_fmp(g):
    """
    Print the statistics from an fmp file
    """
    results = fep_stats.calculate(g)
    for item in results:
        print(item, results[item])

import argparse
def get_arg_parser():
    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program to automatically apply the pKa correction to the scytalone"
                                                 "dehydratase map.")
    parser.add_argument('-i', '--infile', type=str,
                        help="the name of the scytalone out fmp file that.")
    parser.add_argument('-o', '--outfile', type=str,
                        help="the name of the output fmp that has the pKa correction appied.")
    parser.add_argument('-d', '--dict', type=str, choices=['macro', 'fitted'],
                        help="the type of pKa dictionary to use: 'macro' for the pKas predicted by macropka and 'fitted'"
                             " for the pKas that have been optimized to reduce pairwise RMSE, default=macro", default='macro')
    parser.add_argument('--optimize', action='store_true',
                        help='Whether to minimize the pairwise RMSE by altering the pKas of 2d, 3d, 6d, and 8d.', default=False)
    return parser


if __name__ == '__main__':
    args = get_arg_parser().parse_args()

    # The prediction from macropKa on simplified versions of the ligands.
    #PKA_DICT = {'2d': 6.33, '3d': 9.26, '4d': 5.45, '5d': 4.06, '6d': 7.09, '7d': 5.42, '8d': 7.13}
    # 14th January 2019
    # -----------------
    # MacropKa predicts a value of 9.26 for compound 3d. However, there is good experimental evidence that the pKa of
    # 3d is 8.0. So using that value for ligand 3d.
    #PKA_DICT = {'2d': 6.33, '3d': 8.0, '4d': 5.45, '5d': 4.06, '6d': 7.09, '7d': 5.42, '8d': 7.13}
    PKA_DICT = {'2d': 5.60, '3d': 8.89, '4d': 5.23, '5d': 4.06, '6d': 7.09, '7d': 3.61, '8d': 7.10}
    # The result from optimizing the pKas of 2d, 3d, 6d and 8d to minimize the pairwise RMSE of the map.
    OPT_PKA_DICT = {'2d': 6.11923504, '3d': 8.26555333, '4d': 5.45, '5d': 4.06, '6d': 8.0862216, '7d': 5.42, '8d': 7.62419523}

    # Read in the graph:
    g = graph.Graph.deserialize(args.infile)

    if args.dict == 'macro':
        pka_dict = PKA_DICT
    elif args.dict == 'fitted':
        pka_dict = OPT_PKA_DICT
    else:
        raise Exception('Unable to recognize the option for --dict. Please choose a valid option.')

    if args.optimize:
        # Minimize the pKas from the starting dict subject to the constraint that they can be no more than 1 log unit away from their starting value.
        X0 = np.array((pka_dict['2d'], pka_dict['6d'], pka_dict['8d']))
        def pairwise_error_constraint(pkas):
            """
            Defining the function that will be minimized.
            """
            if np.any(np.abs(X0 - pkas) > 1):
                return np.inf
            else:
                g_pka_corr = deepcopy(g)
                pka_dict_mod = {'2d': pkas[0], '3d': 8.0, '4d': 5.45, '5d': 4.06, '6d': pkas[1], '7d': 5.42, '8d': pkas[2]}
                apply_solvent_penalty(g_pka_corr, pka_dict_mod)
                r = fep_stats.calculate(g_pka_corr)
                return r['RMSE Edgewise'].val

        print('\nOptimizing pKas of 2d, 6d and 8d to minimize the edgewise RMSE...')
        result = minimize(fun=pairwise_error_constraint, x0=X0, method='Nelder-Mead', options={'maxiter': 100})
        pka_dict['2d'] = result.x[0]
        pka_dict['6d'] = result.x[1]
        pka_dict['8d'] = result.x[2]
        print('... Optimization compled.')
        print('Fitted values are:')
        print(pka_dict)

    # Apply the pKa correction:
    g_pka_corr = deepcopy(g)
    apply_solvent_penalty(g_pka_corr, pka_dict)
    #print('\nNew results:')
    #pretty_print_fmp(g_pka_corr)
    g_pka_corr.write(args.outfile)
