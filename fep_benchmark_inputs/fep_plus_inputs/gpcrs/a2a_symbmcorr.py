from schrodinger.application.scisol.packages.fep import graph
import numpy as np
import argparse

import sys
sys.path.append('../')
import binding_mode_correction as bmc

# Add the symmetry corrections
degree_of_symm = {'4a':4, '4a flipped':4,
             '4b D':2,
             '4b U':2,
             '4c':4,
             '4d':4,
             '4e D':2, '4eD deprotonated':2,
             '4e U':2, '4eU deprotonated':2,
             '4f':4, '4f deprotonated':4,
             '4g':4,
             '4h':4,
             '4i':4,
             '4jR':2,
             '4jS':2,
             '4kD':2,
             '4kU':2,
             '4lD':2,
             '4lU':2,
             '4m D':2,
             '4m U':2,
             '4n':4,
             '4o':4, '4o deprotonated':4,
             '4q':4,
             '4r D':2,
             '4r U':2}

def get_lignames(edge):
    sep = edge.short_id_title.split('==>')
    lig1 = sep[0].split(':')[1].strip()
    lig2 = sep[1].split(':')[1].strip()
    return lig1, lig2



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

    # Remove edges if they don't have predictions:
    remove_edges = []
    for e in g.edges_iter():
        if e.complex_dg is None or e.solvent_dg is None:
            remove_edges.append(e)

    if len(remove_edges) > 0:
        print('Removing {} edges as they do not have results'.format(len(remove_edges)))
        g.remove_edges_from(remove_edges)

    remove_nodes = []
    for n in g.nodes_iter():
        if n.pred_dg is None:
            remove_nodes.append(n)

    if len(remove_nodes) > 0:
        print('Removing {} nodes as they do not have results'.format(len(remove_nodes)))
        g.remove_nodes_from(remove_nodes)

    kT =0.596
    for e in g.edges_iter():
        l1, l2 = get_lignames(e)
        total_corr = 0
        # Symmetry correction
        total_corr += kT * np.log(degree_of_symm[l1])
        total_corr -= kT * np.log(degree_of_symm[l2])
        e.complex_dg += total_corr

    g.calc_cycle_closure()

    mode_nodes, corrections = bmc.correct_multiple_binding_modes(g)
    bmc.merge_ligand_nodes(g, mode_nodes, corrections)
    g.calc_cycle_closure()
    g.write(args.outfile)

if __name__ == '__main__':
    main()
