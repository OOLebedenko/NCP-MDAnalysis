from typing import Iterable, Callable
from MDAnalysis import Universe, AtomGroup
from MDAnalysis.analysis import dssp


def by_pattern_selector(pattern: str) -> Callable:
    """
    Create a selector function that selects atoms based on a MDAnalysis selection pattern.

    This is a higher-order function that returns a callable selector function.
    Useful for creating partner selectors for ResidueBinaryContactMapper.

    Parameters
    ----------
    pattern : str
        MDAnalysis selection string (e.g., "protein and segid A", "resid 1:100",
        "name CA", "around 5.0 resid 50").

    Returns
    -------
    Callable
        A function that takes an AtomGroup and returns a new AtomGroup containing
        atoms matching the selection pattern.

    Examples
    --------
    >>> # Create a selector for chain A
    >>> chain_a_selector = by_pattern_selector("protein and segid A")
    >>>
    >>> # Use in ResidueBinaryContactMapper
    >>> mapper = ResidueBinaryContactMapper(
    ...     u.atoms,
    ...     selector_partner_1=by_pattern_selector("protein and segid A"),
    ...     selector_partner_2=by_pattern_selector("protein and segid B"),
    ...     distance_cutoff=4.5
    ... )
    >>>
    >>> # Use selector directly on an AtomGroup
    >>> selected_atoms = chain_a_selector(u.atoms)
    >>> print(f"Selected {selected_atoms.n_atoms} atoms")

    Notes
    -----
    The returned function uses MDAnalysis' select_atoms() method internally.
    Refer to MDAnalysis documentation for valid selection patterns.
    """

    def select(ag: AtomGroup) -> AtomGroup:
        """
        Select atoms from an AtomGroup using the predefined pattern.

        Parameters
        ----------
        ag : AtomGroup
            The AtomGroup to select atoms from.

        Returns
        -------
        AtomGroup
            New AtomGroup containing only atoms matching the selection pattern.

        Raises
        ------
        SelectionError
            If the selection pattern is invalid or selects no atoms.
        """
        # Use MDAnalysis' selection engine to select atoms based on the pattern
        return ag.atoms.select_atoms(pattern)

    return select


def get_sec_str_pattern(reference: Universe,
                        chain_ids: Iterable[str]
                        ) -> str:
    """
    :param reference: srtucture to assign secondary-structured elements
    :param chain_ids: list of chain names: ["A"] or ["A", "B"]
    :return: string pattern for selection of atoms from secondary structured regions in specified chains
    """
    sec_str_patterns = []
    for chain in chain_ids:

        dssp_results = dssp.DSSP(reference.select_atoms(f"chainID {chain}")).run().results
        ss_indexes = dssp_results.dssp_ndarray[0][:, 1:].sum(axis=1).astype(bool)
        ss_residues = dssp_results.resids[ss_indexes]

        if len(ss_residues) > 0:
            ss_selection = f"(chainID {chain} and resid {' '.join(ss_residues.astype(str))})"
            sec_str_patterns.append(ss_selection)

    return f"({' or '.join(sec_str_patterns)})" if sec_str_patterns else ""


def get_sec_str_ca_pattern(reference: Universe,
                           chain_ids: Iterable[str]
                           ) -> str:
    """
    :param reference: srtucture to assign secondary-structured elements
    :param chain_ids: list of chain names: ["A"] or ["A", "B"]
    :return: string pattern for selection of Ca atoms from secondary structured regions in specified chains
    """
    selection_sec_str = get_sec_str_pattern(reference=reference,
                                            chain_ids=chain_ids)
    return f"name CA and {selection_sec_str}"
