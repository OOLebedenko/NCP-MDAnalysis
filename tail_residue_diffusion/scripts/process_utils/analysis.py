import numpy as np
from typing import Callable

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis import AtomGroup


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


class ExtractAtomCoords(AnalysisBase):
    def _prepare(self):
        for atom in self.atoms:
            atom_name = self.atomname_provider(atom)
            self.results.update({atom_name: np.zeros((self.n_frames, 5))})
            self.atom_names.append(atom_name)

    def __init__(self,
                 ag: AtomGroup,
                 atomname_provider: Callable,
                 pattern='all',
                 **kwargs):
        """
        :param atomgroup: group of atoms associated with trajectory
        :param selector: selector of atom pairs: atom_1 (origin of vector) and atom_2 (end of vector)
        :param atomname_provider: callable object (function) to set output names for files with vector coordinates
        """
        super(ExtractAtomCoords, self).__init__(ag.universe.trajectory,
                                                **kwargs)
        self.atomname_provider = atomname_provider
        self.atoms = ag.select_atoms(pattern)
        self.atom_names = []

    def _single_frame(self):

        for atom_name, atom in dict(zip(self.atom_names, self.atoms)).items():
            self.results[atom_name][self._frame_index, :2] = self._ts.frame, self._trajectory.time
            self.results[atom_name][self._frame_index, 2:] = atom.position
