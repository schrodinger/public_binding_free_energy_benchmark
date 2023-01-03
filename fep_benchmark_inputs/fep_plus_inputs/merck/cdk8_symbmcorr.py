from schrodinger.application.scisol.packages.fep import graph
import numpy as np
import argparse
from itertools import combinations


import sys
sys.path.append('../')
import binding_mode_correction as bmc


def get_lignames(edge):
    sep = edge.short_id_title.split('==>')
    lig1 = sep[0].split(':')[1].strip()
    lig2 = sep[1].split(':')[1].strip()
    return lig1, lig2

def get_pairwise_diffs(g):
    pred_dg = np.zeros(len(g))
    exp_dg = np.zeros(len(g))
    for i, n in enumerate(g.nodes_iter()):
        pred_dg[i] = n.pred_dg
        exp_dg[i] = n.exp_dg

    #print('Error before', np.mean(np.abs(pred_dg - exp_dg)))
    # Make sure the predictions are centered on the experimental dgs
    pred_dg += np.mean(exp_dg) - np.mean(pred_dg)
    #print('Error after', np.mean(np.abs(pred_dg - exp_dg)))

    indices = np.arange(len(g))
    pred_pairwise = np.array([pred_dg[i] - pred_dg[j] for i, j in combinations(indices, 2)])
    exp_pairwise = np.array([exp_dg[i] - exp_dg[j] for i, j in combinations(indices, 2)])

    return exp_pairwise, pred_pairwise


def get_error_before_corrs(g):
    pred_dgs = []
    exp_dgs = []
    for n in g.nodes_iter():
        ligname = n.short_id_title.split(":")[1].strip()
        if 'flipped' not in ligname:
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
    #print('Pairwise RMSD without extra rotamers and symmetry correction = {:.2f} +/- {:.2f} kcal/mol'.format(pair_rmsd, np.std(boot_rmsds)))


# Add the symmetry corrections. Specific to the map
no_sym_corr = ['13','13 flipped', '14', '15', '16', '16 flipped','17', '17 flipped', '18', '18 flipped', '42','43','43 flipped', '44', '44 flipped', '45','45 flipped']

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

    get_error_before_corrs(g)

    # Remove edges if they don't have predictions:
    remove_edges = []
    for e in g.edges_iter():
        if e.complex_dg is None or e.solvent_dg is None:
            remove_edges.append(e)
    
    if len(remove_edges) > 0:
        #print('Removing {} edges as they do not have results'.format(len(remove_edges))) 
        g.remove_edges_from(remove_edges) 

    remove_nodes = []
    for n in g.nodes_iter():
        if n.pred_dg is None:
            remove_nodes.append(n)
    
    if len(remove_nodes) > 0:
        #print('Removing {} nodes as they do not have results'.format(len(remove_nodes)))
        g.remove_nodes_from(remove_nodes)


    # SYMMETRY CORRECTION
    kT =0.596
    for e in g.edges_iter():
        l1, l2 = get_lignames(e)
        total_corr = 0
        # Symmetry correction
        if l1 not in no_sym_corr:
            total_corr += kT * np.log(2)
        if l2 not in no_sym_corr:
            total_corr -= kT * np.log(2)
        e.complex_dg += total_corr

    g.calc_cycle_closure()

    # BINDING MODE CORRECTION
    mode_nodes, corrections = bmc.correct_multiple_binding_modes(g)
    bmc.merge_ligand_nodes(g, mode_nodes, corrections)
    g.calc_cycle_closure()

    if args.outfile is not None:
        g.write(args.outfile)

if __name__ == '__main__':
    main()
