import numpy as np

from MDAnalysis.analysis import align
from MDAnalysis.transformations import fit_rot_trans
from MDAnalysis.transformations.base import TransformationBase


class TransformWrapper:
    def __init__(self, transform, *args, **kwargs):
        self.transfom = transform

        self.__dict__.update(kwargs)
        self.attr_names = list(kwargs.keys())

        for arg in args:
            setattr(self, arg, arg)
            self.attr_names.append(arg)

    def __call__(self, ag):
        attrs = [getattr(self, arg) for arg in self.attr_names]
        return self.transfom(ag, *attrs)


def metric_tensor(v1, v2, v3):
    basis = np.array([v1, v2, v3])
    g_ij = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            g_ij[i, j] = np.dot(basis[i], basis[j])
    return g_ij


class Alignment:

    def __init__(self, mobile, ref, weights=None):
        self.mobile = mobile
        self.ref = ref
        self.weights = weights

    @property
    def matrix3d(self):
        mobile0 = self.mobile.positions - self.mobile.atoms.center_of_mass()
        ref0 = self.ref.positions - self.ref.atoms.center_of_mass()
        R, _ = align.rotation_matrix(mobile0, ref0, weights=self.weights)
        return R

    @property
    def vector3d(self):
        return -np.dot(self.mobile.center_of_mass(), self.matrix3d.T) + self.ref.center_of_mass()


class AssembleQuaternaryStructure(TransformationBase):
    """Assembles multichain structure that may appear divided at the boundaries
    due to periodic boundary condition based on reference configuration"""

    def __init__(self,
                 ag,
                 reference,
                 chain_ids,
                 atom_selector,
                 max_threads=1,
                 parallelizable=True
                 ):

        super().__init__(
            max_threads=max_threads, parallelizable=parallelizable
        )

        self.ag = ag
        self.reference = reference
        self.atom_selector = atom_selector

        self.ref_mol_selection = []
        self.mobile_mol_selection = []
        self.reference_coms = []

        for chain_id in chain_ids:
            ref_mol = reference.select_atoms(f"chainID {chain_id}")
            ref_sel = ref_mol.select_atoms(atom_selector)
            self.ref_mol_selection.append(ref_sel)

            mol = ag.select_atoms(f"chainID {chain_id}")
            mobile_sel = mol.select_atoms(atom_selector)
            self.mobile_mol_selection.append(mobile_sel)

            self.reference_coms.append(ref_sel.center_of_mass())

        self.reference_coms = np.array(self.reference_coms)

    def _transform(self, ts):
        """Apply transformation to timestep"""
        # 1. Align first chain to reference
        alignment = Alignment(self.ref_mol_selection[0],
                              self.mobile_mol_selection[0])

        # 2. Transform reference coordinates using computed alignment
        ref_mean_coords = self.reference_coms @ alignment.matrix3d.T + alignment.vector3d

        # 3. Apply periodic correction to other chains
        v1, v2, v3 = ts.triclinic_dimensions
        g_ij = metric_tensor(v1, v2, v3)  # Compute metric tensor once per frame
        g_inv = np.linalg.inv(g_ij)  # Invert once per frame

        # shift rest of molecules to match reference coordinates as close as possible
        for mol_idx in range(1, len(self.mobile_mol_selection)):

            # Current center of mass
            current_com = self.mobile_mol_selection[mol_idx].center_of_mass()

            # Displacement vector
            delta = ref_mean_coords[mol_idx] - current_com

            # Calculate periodic shifts
            projections = np.array([delta.dot(v1), delta.dot(v2), delta.dot(v3)])
            i, j, k = np.round(g_inv @ projections).astype(int)

            # Apply shift if needed
            if np.any((i, j, k)):
                shift = i * v1 + j * v2 + k * v3
                ts.positions[self.mobile_mol_selection[mol_idx][0].segment.atoms.ix] += shift
        return ts


class fit_rot_trans_by_pattern(fit_rot_trans):

    def __init__(self, ag, reference, pattern='all', *args, **kwargs):
        super().__init__(ag, reference, *args, **kwargs)

        self.mobile_all_atoms_com = self.mobile.atoms.center(self.weights)

        ag_ca = ag.select_atoms(pattern)
        ref_ca = reference.select_atoms(pattern)
        self.ref_com_all = self.ref.center(self.weights)

        self.ref, self.mobile = align.get_matching_atoms(ref_ca.atoms,
                                                         ag_ca.atoms)

        self.weights = align.get_weights(self.ref.atoms,
                                         weights=self.weights)

        self.ref_com = self.ref.center(self.weights)
        self.ref_coordinates = self.ref.atoms.positions - self.ref_com
