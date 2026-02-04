import argparse
import os

from glob import glob

from MDAnalysis import Universe

from process_utils.select import get_sec_str_ca_pattern
from process_utils.analysis import Hist3dAnalysis
from process_utils.transform import AssembleQuaternaryStructure, TransformWrapper, HistonePseudosymmetryAlign
from process_utils.batch_process import BatchLoader, BatcHist3dWriter, BatchAnalyzer


class AnalyzerWrapper:
    def __init__(self, analysis, *args, **kwargs):
        self.analysis = analysis

        self.__dict__.update(kwargs)
        self.attr_names = list(kwargs.keys())

        for arg in args:
            setattr(self, arg, arg)
            self.attr_names.append(arg)

    def __call__(self, ag):
        attrs = [getattr(self, arg) for arg in self.attr_names]
        return self.analysis(ag, **dict(zip(self.attr_names, attrs)))


def rotate_reference(reference: Universe,
                     histone_chains_1,
                     histone_chains_2,
                     dna_chains
                     ) -> Universe:
    reference_rotated = reference.copy()

    # reference_rotated.segments.segids = ["E", "F", "G", "H", "A", "B", "C", "D", "I", "J"]
    reference_rotated.segments.segids = [*histone_chains_2, *histone_chains_1, *dna_chains]
    chainids = []
    for segment in reference_rotated.segments:
        chainids.extend([segment.segid] * segment.atoms.n_atoms)
    reference_rotated.add_TopologyAttr("chainIDs", chainids)
    reference_rotated.atoms = reference_rotated.atoms.sort(key="segids")

    reference_rotated.atoms.write("tmp.pdb")
    reference_rotated = Universe("tmp.pdb", topology_format="PDB")
    os.remove("tmp.pdb")

    return reference_rotated


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--path-to-xray-reference', required=True)
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

    # set xray reference and rotated xray reference
    xray_reference = Universe(args.path_to_xray_reference, topology_format="PDB")
    rotated_xray_reference = rotate_reference(xray_reference,
                                              histone_chains_1=args.protein_chains[:4],
                                              histone_chains_2=args.protein_chains[4:],
                                              dna_chains=args.dna_chains
                                              )

    # set atom selection to assemble quaternary structure and get rotated copy of histone chain:
    sec_str_ca = get_sec_str_ca_pattern(xray_reference, chain_ids=args.protein_chains)
    dna_pattern = f"(chainID {' '.join(args.dna_chains)}) and (name N1 N9)"
    sec_str_ca_and_dna_pattern = f"({dna_pattern}) or ({sec_str_ca})"

    # set  CA atom selection to calculate 3d hist
    first_rid, last_rid = args.residue_of_interest_histone.split("-")
    resids_of_interest = set(list(range(int(first_rid), int(last_rid) + 1)))
    pattern_h4_tail = f"chainID {args.histone_chains.replace(',', ' ')} " \
                      f" and resid {' '.join(list(map(str, resids_of_interest)))} and (name CA)"

    #  load trajectory reference
    trj_reference = Universe(args.path_to_trajectory_reference, topology_format="PDB")

    # set path to trajectory files
    nc_files = glob(os.path.join(args.path_to_trajectory, '*.nc'))
    nc_files.sort()

    # set trajectory transforms
    transforms = [
        # 1. assemble the nucleosome particle that may appear divided at the boundaries
        # due to periodic boundary condition in the MD simulation.
        TransformWrapper(transform=AssembleQuaternaryStructure,
                         reference=trj_reference,
                         chain_ids=args.protein_chains,
                         atom_selector=sec_str_ca_and_dna_pattern),
        # 2. overlay all MD frames by superimposing them onto the reference
        # via the secondary-structure Cα atoms from the histone core.
        TransformWrapper(transform=HistonePseudosymmetryAlign,
                         reference=xray_reference,
                         rotated_reference=rotated_xray_reference,
                         pattern=sec_str_ca)
    ]

    # set loader for trajectory to batch processing due to file open limits
    batchloader = BatchLoader(reference=trj_reference,
                              trj_list=nc_files[args.trajectory_start:args.trajectory_length],
                              trajectory_stride=args.trajectory_stride,
                              batch_size=args.batch_size,
                              dt_ns=args.dt_ns,
                              transforms=transforms,
                              pattern=pattern_h4_tail,
                              )

    # set trajectory analyzer to extract coordinates of NH vectors
    analyzer = AnalyzerWrapper(Hist3dAnalysis,
                               delta=1,
                               gridcenter=xray_reference.atoms.center_of_mass(),
                               xdim=xray_reference.dimensions[:3].max() / 2 + 100,
                               ydim=xray_reference.dimensions[:3].max() / 2 + 100,
                               zdim=xray_reference.dimensions[:3].max() / 2 + 100,
                               )

    # set writer to save the coordinates of NH vectors
    writer = BatcHist3dWriter(output_directory=args.output_directory)

    # process batches of trajectory and save the results
    batchanalyzer = BatchAnalyzer(batchloader=batchloader,
                                  analyzer=analyzer,
                                  writer=writer)
    batchanalyzer.analyse()
