import argparse
import os

from glob import glob

from MDAnalysis import Universe
from MDAnalysis.core.groups import Atom

from process_utils.analysis import AnalyzerWrapper, ExtractAtomCoords
from process_utils.select import get_sec_str_ca_pattern
from process_utils.transform import AssembleQuaternaryStructure, fit_rot_trans_by_pattern, TransformWrapper
from process_utils.batch_process import BatchLoader, BatchCsvWriter, BatchAnalyzer


def atomname_provider(atom: Atom):
    return f"{atom.segid}_{atom.residue.resnum:03d}_{atom.residue.resname}_{atom.name}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract CA coords')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--histone-chains', required=True, type=str)
    parser.add_argument('--residue-of-interest-histone', required=True)
    parser.add_argument('--trajectory-start', default=0, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--trajectory-stride', default=100, type=int)
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
        # 2. overlay all MD frames by superimposing them onto the reference
        # via the secondary-structure Cα atoms from the histone core.
        TransformWrapper(transform=fit_rot_trans_by_pattern,
                         reference=trj_reference,
                         pattern=sec_str_ca)
    ]

    # set loader for trajectory to batch processing due to file open limits
    first_rid, last_rid = args.residue_of_interest_histone.split("-")
    resids_of_interest = set(list(range(int(first_rid), int(last_rid) + 1)))
    pattern_h4_tail = f"chainID {args.histone_chains.replace(',', ' ')} " \
                      f" and resid {' '.join(list(map(str, resids_of_interest)))} and (name CA)"

    batchloader = BatchLoader(reference=trj_reference,
                              trj_list=nc_files[args.trajectory_start:args.trajectory_length],
                              trajectory_stride=args.trajectory_stride,
                              batch_size=args.batch_size,
                              dt_ns=args.dt_ns,
                              transforms=transforms,
                              pattern=pattern_h4_tail
                              )

    # set trajectory analyzer to calc msd
    analyzer = AnalyzerWrapper(ExtractAtomCoords,
                               atomname_provider=atomname_provider,
                               pattern='all',
                               )

    # set writer to save the coordinates of NH vectors
    writer = BatchCsvWriter(output_directory=args.output_directory,
                               header=["time_ns", "x", "y", "z"]
                               )

    # process batches of trajectory and save the results
    batchanalyzer = BatchAnalyzer(batchloader=batchloader,
                                  analyzer=analyzer,
                                  writer=writer)
    batchanalyzer.analyse()
