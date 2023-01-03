import argparse
import sys
from collections import defaultdict
from typing import List

import numpy as np

from schrodinger.application.desmond.constants import BOLTZMANN
from schrodinger.application.desmond.measurement import Measurement
from schrodinger.application.scisol.packages.fep import graph
from schrodinger.structutils import analyze


def clean_bennett_from_map(g: graph.Graph):
    """
    Removes the bennett complex and solvent free energy data and replace with the cycle closure correction DDG
    predictions. This helps to simplify the application of the binding mode correction.

    :param g: An FEP+ graph with output data.

    Note
    -----
    The graph is modified in place.
    """
    # TODO: As this creates an inconsistency between the DDGs in the graph and the SID, the SID data is removed.
    # The Jira ticket DESMOND-10777 has been raised to discuss the possibility of updating the SID in the future.
    for e in g.edges_iter():
        e.complex_dg = e.ccc_ddg
        e.solvent_dg = Measurement(0, 0)
        e.vacuum_dg = None
        e.complex_sid = None
        e.solvent_sid = None


def get_binding_mode_nodes(g: graph.Graph) -> List[List[graph.Node]]:
    """
    Extract the nodes that represent different binding modes of the same
    ligand. SMILES are used to identify which ligands are the same.

    :param g: A graph that contains multiple binding modes of the ligands, with
        each binding mode represented as a separate node.

    :return: The nodes of the binding modes as a list of lists. Each entry in
        the list contains all the nodes of the different binding modes of the
        same ligand.
    """
    # Extract the unique SMILES of each node in the graph and place matching SMILES into the same dict entry
    mode_nodes_dict = defaultdict(list)
    for n in g.nodes_iter():
        smile = analyze.generate_smiles(n.struc, unique=True)
        mode_nodes_dict[smile].append(n)

    # Return the nodes whose SMILES appear multiple times
    return [val for val in mode_nodes_dict.values() if len(val) > 1]


def calc_mode_corrections(mode_nodes: List[List[graph.Node]],
                          temperature: float = 300.) -> List[np.ndarray]:
    """
    Calculate the binding mode correction for each mode and for each ligand.
    The presence of other binding modes serves to increase the binding strength
    of a ligands.

    :param mode_nodes: A list of lists that contains the nodes of the different
        binding modes of the same ligand.
    :param temperature: The temperature of the simulation in Kelvin.

    :returns: The binding mode corrections for each ligand. Each correction
        accounts for the presence of the other nodes of the same ligand.
    """
    kT = BOLTZMANN * temperature
    corrections = []
    for nodes in mode_nodes:
        dgs = np.array([n.pred_dg.val for n in nodes])
        corrs = np.zeros(len(dgs))
        for i in range(len(dgs)):
            corrs[i] = -kT * np.log(np.sum(np.exp(-(dgs - dgs[i]) / kT)))
        corrections.append(corrs)

    return corrections


def apply_corrections(mode_nodes: List[List[graph.Node]],
                      corrections: List[np.ndarray]):
    """
    Apply the binding mode correction to each edge that is connected to a node that represents one of at least two
    binding modes of the same ligand. The edge data are modified in place. and nothing is returned.

    :param mode_nodes: A list of lists that contains the nodes of the different binding modes of the same ligand.
    :param corrections: The binding mode corrections for each ligand.
    """
    for nodes, correction in zip(mode_nodes, corrections):
        # There is a different correction for each binding mode of the same ligand.
        for n, corr in zip(nodes, correction):
            # Apply correction to every edge the node is part of.
            for e in n.edges:
                if n == e[0]:
                    # Deduct the binding mode correction
                    e.complex_dg -= Measurement(corr, 0)
                else:
                    # Add the correction
                    e.complex_dg += Measurement(corr, 0)


def correct_multiple_binding_modes(
        g: graph.Graph) -> (List[List[graph.Node]], List[np.ndarray]):
    """
    Automatically detects different modes, calculates the free energy correction, and applies it to the map. The input
    map is modified in place.

    :param g: An FEP+ graph that contains multiple binding modes.

    :returns mode_nodes: A list of lists that contains the nodes of the different binding modes of the same ligand.
    :returns corrections: The binding mode corrections for each ligand.
    """
    # Simplify the map by replacing all the bennett DG information with the cycle closure corrected DGs
    clean_bennett_from_map(g)
    # Find which ligands have multiple binding modes using canonicalized SMILES:
    mode_nodes = get_binding_mode_nodes(g)
    # Calculate the binding mode corrections and apply to the map:
    corrections = calc_mode_corrections(mode_nodes)
    apply_corrections(mode_nodes, corrections)
    g.calc_cycle_closure()

    return mode_nodes, corrections


def merge_ligand_nodes(g: graph.Graph, mode_nodes: List[List[graph.Node]],
                       corrections: List[np.ndarray]):
    """
    Merge all nodes of different binding modes into a single one. The most
    favourable binding mode of each ligand is retained. The graph is modified
    in place.

    :param g: An FEP+ graph that contains multiple binding modes. Cycle closure
        must have called for the map for this function to run.
    :param mode_nodes: A list of lists that contains the nodes of the different
        binding modes of the same ligand.
    :param corrections: The binding mode corrections for each ligand.
    """
    for n in g.nodes_iter():
        if n.pred_dg is None:
            raise ValueError('Cycle closure must be called for the input map.')

    # Determine which nodes to keep and which nodes to discard.
    # The mode that requires the least correction is the most favourable and it will be retained.
    to_merge = []
    to_keep = []
    for c in corrections:
        keep_ind = np.where(c == np.max(c))[0][0]
        to_merge.append(np.where(c != c[keep_ind])[0])
        to_keep.append(keep_ind)

    # Add the edges of the binding modes into the one that will be retained.
    nodes_to_delete = []
    nodes_to_keep = []
    for nodes, merge_inds, keep_inds in zip(mode_nodes, to_merge, to_keep):
        for i in merge_inds:
            # Add the neighbours of the other binding modes to the node that will be kept.
            for n in nodes[i].neighbors:
                g.add_edge(nodes[keep_inds], n)
            nodes_to_delete.append(nodes[i])
            nodes_to_keep.append(nodes[keep_inds])

    # Looping through the nodes that will be kept, going through the edges, and adding the predicted binding data
    for nodes in nodes_to_keep:
        for e in nodes.edges:
            n_initial, n_final = e.direction
            e.complex_dg = n_final.pred_dg - n_initial.pred_dg
            e.solvent_dg = Measurement(0, 0)
            e.vacuum_dg = None

    # Finally deleting the nodes of the duplicate binding modes.
    g.remove_nodes_from(nodes_to_delete)

    return nodes_to_keep


def main(argv=None):
    usage = """
    Let map_out.fmp be an output from FEP+ that has different binding modes of some of the ligands. To
    calculate the relative binding free between ligands in a way that accounts for the presence of multiple binding
    modes, type

    $SCHRODINGER/run -FROM scisol binding_mode_correction.py map_out.fmp -o map_bmcorr_out.fmp

    By default, the output 'map_bmcorr_out.fmp' has nodes of the same ligand merged into one. Edges between merged nodes
    do not have trajectory (fmpdb) information. To not merge the nodes, use the flag '-no-merge', like this

    $SCHRODINGER/run -FROM scisol binding_mode_correction.py map_out.fmp -o map_bmcorr_out.fmp -no-merge
    """
    description = """
    In FEP+, one can run maps that have different binding modes of ligands in different nodes. This script automatically
    detects multiple binding modes in an FEP+ map and applies the binding mode correction to the calculated relative free
    energies. When using this script, there are no limitations to the number of ligands and poses, or on the geometry of the map.

    To aid visualization and MUE and RMSE calculation, nodes of the same ligand are automatically merged into a node
    that represents the most favorable binding mode, although this feature can be switched off.

    When applying the corrections, the Bennett DDGs are replaced with cycle closure corrected DDGs. As this creates an
    inconsistency between the Bennett DDGs in the map and SID data, the SID data is removed from the map.
    """
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument(
        'infile',
        type=str,
        help="The name of an out fmp file that contains multiple binding modes."
    )
    parser.add_argument(
        '-o',
        dest='outfile',
        type=str,
        help=
        "The name of the output fmp file that has the binding mode correction for each edge, default='bmcorr_out.fmp'.",
        default='bmcorr_out.fmp')
    parser.add_argument(
        '-no-merge',
        action='store_true',
        dest='no_merge',
        help=
        'Option to NOT merge nodes of the same binding mode into most favorable, default=False.',
        default=False)
    args = parser.parse_args(argv)

    # Load the map that contains ligands with at least 2 binding modes:
    g = graph.Graph.deserialize(args.infile)

    if g.fep_type != graph.FEP_TYPES.SMALL_MOLECULE:
        sys.exit("Only small molecule graphs are supported")

    # Detect alternate binding modes and apply the correction:
    mode_nodes, corrections = correct_multiple_binding_modes(g)

    if not args.no_merge:
        merge_ligand_nodes(g, mode_nodes, corrections)
        g.calc_cycle_closure()

    g.write(args.outfile)


if __name__ == '__main__':
    main()
