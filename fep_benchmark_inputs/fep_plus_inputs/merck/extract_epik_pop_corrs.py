import numpy as np
from schrodinger.application.scisol.packages.fep.graph import Graph
import argparse
from collections import defaultdict

epik_penalty_names = ('r_epik_Ionization_Penalty',
                      'r_epik_Ionization_Penalty_Charging',
                      'r_epik_Ionization_Penalty_Neutral',
                      'r_epik_State_Penalty',
                      'r_epik_Charging_Adjusted_Penalty',
                      'i_epik_Tot_Q',
                      'i_epik_Tot_abs_Q',
                      'r_lp_tautomer_probability')
 
def extract_supernodes(g):
    smile2node = defaultdict(list)
    for n in g.nodes_iter():
        if 'r_epik_Ionization_Penalty' in n.struc.property:
            smile2node[n.canonical_smiles].append(n)

    smile2supernode = defaultdict(list)
    for key in smile2node:
        if len(smile2node[key]) > 1:
            smile2supernode[key].extend(smile2node[key])
    
    return smile2supernode


def get_ligname(node):
    return node.short_id_title.split(':')[1].strip()


def write_popfile_from_epik(g, delimiter, outname, kT=0.596):
    """
    Parse the epik ionization penalty for each ligand in the graph that has muliple protomers/tautomers and 
    write out the input for the pka_tautomer_correction.py script.
    """
    # Find the supernodes in the graph:
    smile2supernode = extract_supernodes(g)
    lines = []
    for key in smile2supernode:
        penalties = np.array([n.struc.property['r_epik_Ionization_Penalty'] for n in smile2supernode[key]])
        penalties -=  penalties.min()   # Ionization penalties already reflect normalized weights, but doing this for stability
        populations = np.exp(-penalties / kT ) / np.sum(np.exp(-penalties / kT ))
        pop_str = [str(p) for p in populations]
        lignames = [get_ligname(n) for n in smile2supernode[key]]

        lines.append(delimiter.join(lignames) + delimiter + delimiter.join(pop_str))
            
    with open(outname, 'w') as f:
        f.writelines('\n'.join(lines))


def main(argv=None):
    usage="""
    Extract the Epik protomer/tautomer ionization state penalties and write out a population file for the pKa tautomer 
    correction script.

    $SCHRODINGER/run extract_epik_pop_corrs.py map_out.fmp -o epik_populations.txt
    """
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('infile',
                        type=str,
                        help="The name of an fmp whose ligands have epik state penalties as molecular properties.")
    parser.add_argument('-o',
                        dest='outfile',
                        type=str,
                        help="The name of the population file that will be produced, default=epik_populations.txt.",
                        default='epik_populations.txt')
    parser.add_argument('-d',
                        dest='delimiter',
                        type=str,
                        help="The delimiter with which to space the text in the population file, default=-",
                        default='-')
    args = parser.parse_args(argv)

    g = Graph.deserialize(args.infile)
    write_popfile_from_epik(g, outname=args.outfile, delimiter=args.delimiter)


if __name__ == '__main__':
    main()
