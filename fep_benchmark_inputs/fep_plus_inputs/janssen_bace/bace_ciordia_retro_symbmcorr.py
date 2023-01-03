from schrodinger.application.scisol.packages.fep import graph
import numpy as np
import argparse
from itertools import combinations

import sys
sys.path.append('/home/ross/Work/schrython/')
import fep_tools as ftools
import binding_mode_correction as bmc


def get_lignames(edge):
    sep = edge.short_id_title.split('==>')
    lig1 = sep[0].split(':')[1].strip()
    lig2 = sep[1].split(':')[1].strip()
    return lig1, lig2


# Add the symmetry corrections. Specific to the map
degree_of_symm = {'17':2}

def get_error_before_corrs(g):
    pred_dgs = []
    exp_dgs = []
    for n in g.nodes_iter():
        ligname = n.short_id_title.split(":")[1].strip()
        if 'flip' not in ligname:
            pred_dgs.append(n.pred_dg.val)
            exp_dgs.append(n.exp_dg.val)

    pred_dgs = np.array(pred_dgs)
    exp_dgs = np.array(exp_dgs)

    pred_dgs += np.mean(exp_dgs) - np.mean(pred_dgs)

    indices = np.arange(len(exp_dgs))
    pred_pairwise = np.array([pred_dgs[i] - pred_dgs[j] for i, j in combinations(indices, 2)])
    exp_pairwise = np.array([exp_dgs[i] - exp_dgs[j] for i, j in combinations(indices, 2)])

    sq_diffs = (exp_pairwise - pred_pairwise)**2
    boot_rmsds = np.zeros(10000)
    for b in range(10000):
        samp_inds = np.random.choice(len(sq_diffs), len(sq_diffs))
        arr_boot = sq_diffs[samp_inds]
        boot_rmsds[b] = np.sqrt(arr_boot.mean())

    pair_rmsd = np.sqrt(np.mean((exp_pairwise - pred_pairwise)**2))
    print('Pairwise RMSD without extra rotamers and symmetry correction = {:.2f} +/- {:.2f} kcal/mol'.format(pair_rmsd, np.std(boot_rmsds)))


def main(argv=None):
    usage="""
    Apply the symmetry and binding mode correction to a particular map.

    $SCHRODINGER/run symbmcorr.py map_out.fmp -o map_bmcorr_out.fmp
    """
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('infile',
                        type=str,
                        help="The name of an out fmp file that contains multiple binding modes.")
    parser.add_argument('-o',
                        dest='outfile',
                        type=str,
                        help="If entered, the name of the corrected map, default=None.",
                        default=None)
    args = parser.parse_args(argv)

    g = graph.Graph.deserialize(args.infile)

    ftools.pretty_print_fmp(g)
    print()
    print(get_error_before_corrs(g))
    print()

    # SYMMETRY CORRECTION
    kT =0.596
    for e in g.edges_iter():
        l1, l2 = get_lignames(e)
        total_corr = 0
        # Symmetry correction
        if l1 in degree_of_symm:
            total_corr += kT * np.log(degree_of_symm[l1])
        if l2 in degree_of_symm:
            total_corr -= kT * np.log(degree_of_symm[l2])
        e.complex_dg += total_corr

    g.calc_cycle_closure()

    # BINDING MODE CORRECTION
    mode_nodes, corrections = bmc.correct_multiple_binding_modes(g)
    bmc.merge_ligand_nodes(g, mode_nodes, corrections)
    g.calc_cycle_closure()

    ftools.pretty_print_fmp(g)

    if args.outfile is not None:
        g.write(args.outfile)

if __name__ == '__main__':
    main()
