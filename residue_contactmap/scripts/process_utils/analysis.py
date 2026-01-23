from scipy.spatial import cKDTree
import numpy as np
from typing import Callable
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis.base import AnalysisBase


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
        return self.analysis(ag, *attrs)


class ResidueBinaryContactMapper(AnalysisBase):
    def __init__(self,
                 ag: AtomGroup,
                 selector_partner_1: Callable,
                 selector_partner_2: Callable,
                 distance_cutoff: float,
                 residuename_provider: Callable,
                 **kwargs):
        super(ResidueBinaryContactMapper, self).__init__(ag.universe.trajectory,
                                                         **kwargs)
        self._ag = ag

        self.selector_partner_1 = selector_partner_1
        self.selector_partner_2 = selector_partner_2

        self.partner_1 = self.selector_partner_1(self._ag)
        self.partner_2 = self.selector_partner_2(self._ag)

        self.distance_cutoff = distance_cutoff
        self.residuename_provider = residuename_provider

    def _prepare(self):
        for residue_partner_1 in self.partner_1.residues:
            key = self.residuename_provider(residue_partner_1)
            self.results[key] = np.zeros((self.n_frames, self.partner_2.n_residues + 2))

    def _single_frame(self):

        tree_partner_1 = cKDTree(self.partner_1.positions)
        tree_partner_2 = cKDTree(self.partner_2.positions)
        dist_dok_matrix = tree_partner_1.sparse_distance_matrix(tree_partner_2, self.distance_cutoff)

        partener_1_index, partener_2_index = dist_dok_matrix.nonzero()
        residue_pairs = set(zip([atom_1.residue for atom_1 in self.partner_1[partener_1_index]],
                                [atom_2.residue for atom_2 in self.partner_2[partener_2_index]]))

        for key in self.results.keys():
            self.results[key][self._frame_index, :2] = self._ts.frame, self._trajectory.time

        for residue_partner_1, residue_partner_2 in residue_pairs:
            key = self.residuename_provider(residue_partner_1)
            self.results[key][self._frame_index, 2:][residue_partner_2.resindex - self.partner_2[0].resindex] = 1
