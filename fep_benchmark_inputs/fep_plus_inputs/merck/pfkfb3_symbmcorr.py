from schrodinger.application.scisol.packages.fep import graph
import numpy as np
import argparse


import sys
sys.path.append('../')
import binding_mode_correction as bmc


def get_lignames(edge):
    sep = edge.short_id_title.split('==>')
    lig1 = sep[0].split(':')[1].strip()
    lig2 = sep[1].split(':')[1].strip()
    return lig1, lig2


# Add the symmetry corrections. Specific to the map
degree_of_symm = {'44':2, '46':2, '47':2, '48':2, '49':2} 


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

    # SYMMETRY CORRECTION
    kT =0.596
    for e in g.edges_iter():
        if e.complex_dg is not None:
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

    if args.outfile is not None:
        g.write(args.outfile)

if __name__ == '__main__':
    main()
