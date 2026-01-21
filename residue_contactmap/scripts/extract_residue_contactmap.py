import argparse
import os

from itertools import chain
from glob import glob

from MDAnalysis import Universe
from MDAnalysis.core.groups import Residue

from process_utils.select import by_pattern_selector, get_sec_str_ca_pattern
from process_utils.analysis import ResidueBinaryContactMapper, AnalyzerWrapper
from process_utils.transform import AssembleQuaternaryStructure, TransformWrapper
from process_utils.batch_process import BatchLoader, BatchCsvWriter, BatchAnalyzer


def residuename_provider(residue: Residue):
    return f"{residue.segid}_{residue.resnum:02d}_{residue.resname}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract NH vectors')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--chain-name-dna', required=True, type=str)
    parser.add_argument('--chain-name-histone', required=True, type=str)
    parser.add_argument('--residue-of-interest-histone', required=True)
    parser.add_argument('--trajectory-start', default=0, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--trajectory-stride', default=1, type=int)
    parser.add_argument('--batch-size', default=100, type=int)
    parser.add_argument('--dt-ns', type=float, default=0.01)
    parser.add_argument('--output-directory', default=".")
    parser.add_argument('--dna-chains', default=["I", "J"], type=list)
    parser.add_argument('--protein-chains', default=["A", "B", "C", "D", "E", "F", "G", "H"], type=list)
    args = parser.parse_args()

    #  load trajectory reference
    trj_reference = Universe(args.path_to_trajectory_reference, topology_format="PDB")

    # set selector to assemble structure that may appear divided at the boundaries due to periodic boundary condition
    sec_str_ca = get_sec_str_ca_pattern(trj_reference, chain_ids=args.protein_chains)
    dna_pattern = f"(name N1 N9)"
    sec_str_ca_and_dna_pattern = f"({dna_pattern}) or ({sec_str_ca})"

    # set path to trajectory files
    nc_files = glob(os.path.join(args.path_to_trajectory, '*.nc'))
    nc_files.sort()

    # set trajectory transforms
    transforms = [
        # 1. assemble the nucleosome particle that may appear divided at the boundaries
        # due to periodic boundary condition in the MD simulation.
        TransformWrapper(transform=AssembleQuaternaryStructure,
                         reference=trj_reference,
                         chain_ids=args.protein_chains + args.dna_chains,
                         atom_selector=sec_str_ca_and_dna_pattern),
    ]

    # set loader for trajectory to batch processing due to file open limits
    batchloader = BatchLoader(reference=trj_reference,
                              trj_list=nc_files[args.trajectory_start:args.trajectory_length],
                              trajectory_stride=args.trajectory_stride,
                              batch_size=args.batch_size,
                              dt_ns=args.dt_ns,
                              transforms=transforms,
                              )

    # set trajectory analyzer to extract coordinates of NH vectors
    # pattern_dna = f"chainID {args.chain_name_dna} and not (name H*)"
    pattern_dna = f"chainID {args.chain_name_dna} and not (name H*)"

    first_rid, last_rid = args.residue_of_interest_histone.split("-")
    resids_of_interest = set(list(range(int(first_rid), int(last_rid) + 1)))
    pattern_h4_tail = f"chainID {args.chain_name_histone} and resid {' '.join(list(map(str, resids_of_interest)))} and not (name H*)"

    analyzer = AnalyzerWrapper(ResidueBinaryContactMapper,
                               selector_partner_1=by_pattern_selector(pattern_h4_tail),
                               selector_partner_2=by_pattern_selector(pattern_dna),
                               distance_cutoff=4.0,
                               residuename_provider=residuename_provider
                               )

    # set writer to save the coordinates of NH vectors
    writer = BatchCsvWriter(output_directory=args.output_directory,
                            header=["time_ns", *[f"{residue.segid}{residue.resnum:03d}-{residue.resname}"
                                                 for residue in trj_reference.select_atoms(pattern_dna).residues]]
                            )

    # process batches of trajectory and save the results
    batchanalyzer = BatchAnalyzer(batchloader=batchloader,
                                  analyzer=analyzer,
                                  writer=writer)
    batchanalyzer.analyse()
