"""Protein-ligand interaction detection via ``peppr.ContactMeasurement``.

:meth:`InteractionSet.find` wraps peppr's per-interaction-type methods
into a single :class:`InteractionSet` carrying receptor/ligand
``AtomArray`` references plus the detected H-bonds, salt bridges,
π-cation interactions, and π-stacking interactions.

The :class:`InteractionSet` can be passed to :class:`portein.Pymol` via
its ``interactions`` field to overlay the interactions on a PyMOL render.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from biotite.structure import AtomArray


@dataclass(frozen=True)
class HBond:
    """Hydrogen bond between a donor heavy atom and an acceptor heavy atom."""

    donor_atom_idx: int
    """Index in the donor side's AtomArray."""
    acceptor_atom_idx: int
    """Index in the acceptor side's AtomArray."""
    donor_is_receptor: bool
    """If True the receptor donates; if False the ligand donates."""


@dataclass(frozen=True)
class SaltBridge:
    """Salt bridge between an oppositely-charged receptor and ligand atom."""

    receptor_atom_idx: int
    ligand_atom_idx: int


@dataclass(frozen=True)
class PiCation:
    """π-cation interaction between an aromatic ring and a cationic atom."""

    cation_atom_idx: int
    """Index of the cation in the side that owns it."""
    ring_atom_indices: tuple[int, ...]
    """Indices of the aromatic ring atoms in the other side."""
    cation_in_receptor: bool
    """If True the cation is in the receptor (and the ring in the ligand)."""


@dataclass(frozen=True)
class PiStacking:
    """π-stacking interaction between two aromatic rings."""

    receptor_ring_atom_indices: tuple[int, ...]
    ligand_ring_atom_indices: tuple[int, ...]
    kind: int
    """``biotite.structure.PiStacking`` enum value:
    1 = parallel/face-to-face, 2 = T-shaped/edge-to-face."""


@dataclass
class InteractionSet:
    """Bundle of detected protein-ligand interactions plus the source structures.

    Atom indices in the individual interaction objects are local to the
    referenced ``receptor`` and ``ligand`` AtomArrays — not into a global
    combined structure. This lets the rendering layer (or any other consumer)
    resolve them to PyMOL selections via the atoms' ``chain_id``, ``res_id``,
    and ``atom_name`` annotations.
    """

    receptor: AtomArray
    ligand: AtomArray
    hbonds: list[HBond] = field(default_factory=list)
    salt_bridges: list[SaltBridge] = field(default_factory=list)
    pi_cation: list[PiCation] = field(default_factory=list)
    pi_stacking: list[PiStacking] = field(default_factory=list)

    @property
    def combined(self) -> AtomArray:
        """``receptor + ligand`` as one AtomArray, ready for ``ProteinConfig.pdb_file``."""
        return self.receptor + self.ligand

    @classmethod
    def find(
        cls,
        receptor: AtomArray,
        ligand: AtomArray,
        *,
        rotate: bool = False,
        cutoff: float = 8.0,
        ph: float = 7.4,
        use_resonance: bool = True,
        use_tautomers: bool = True,
    ) -> InteractionSet:
        """Detect all supported protein-ligand interactions.

        Wraps :class:`peppr.ContactMeasurement`. Both AtomArrays must
        carry bond information — use
        ``biotite.structure.connect_via_residue_names`` if your input
        came from a PDB file without explicit CONECT records.

        Parameters
        ----------
        receptor, ligand
            Protein and ligand atoms. Must contain heavy atoms only and
            have a ``BondList`` attached.
        rotate
            If True, after detecting interactions also rotate the
            receptor and ligand so the 2D projection maximizes the
            spread of the interaction atoms. The returned set points at
            the rotated structures — ``result.combined`` can be passed
            straight to ``ProteinConfig(pdb_file=..., rotate=False)``.
        cutoff
            Binding-site cutoff (Å) passed to ``peppr.ContactMeasurement``.
        ph
            Environmental pH for charge estimation.
        use_resonance, use_tautomers
            Forwarded to ``peppr.ContactMeasurement``.

        Returns
        -------
        InteractionSet
            All detected interactions, with references to the (possibly
            rotated) ``receptor`` and ``ligand`` arrays.
        """
        import peppr

        cm = peppr.ContactMeasurement(
            receptor=receptor,
            ligand=ligand,
            cutoff=cutoff,
            ph=ph,
            use_resonance=use_resonance,
            use_tautomers=use_tautomers,
        )

        receptor_donates, ligand_donates = cm.find_hbonds()
        hbonds: list[HBond] = []
        for rec_i, lig_i in receptor_donates:
            hbonds.append(HBond(donor_atom_idx=int(rec_i), acceptor_atom_idx=int(lig_i), donor_is_receptor=True))
        for rec_i, lig_i in ligand_donates:
            hbonds.append(HBond(donor_atom_idx=int(lig_i), acceptor_atom_idx=int(rec_i), donor_is_receptor=False))

        salt_bridges: list[SaltBridge] = [
            SaltBridge(receptor_atom_idx=int(rec_i), ligand_atom_idx=int(lig_i))
            for rec_i, lig_i in cm.find_salt_bridges()
        ]

        pi_cation: list[PiCation] = []
        for rec_idx, lig_idx, cation_in_receptor in cm.find_pi_cation_interactions():
            if cation_in_receptor:
                pi_cation.append(
                    PiCation(
                        cation_atom_idx=int(rec_idx[0]),
                        ring_atom_indices=tuple(int(i) for i in lig_idx),
                        cation_in_receptor=True,
                    )
                )
            else:
                pi_cation.append(
                    PiCation(
                        cation_atom_idx=int(lig_idx[0]),
                        ring_atom_indices=tuple(int(i) for i in rec_idx),
                        cation_in_receptor=False,
                    )
                )

        pi_stacking: list[PiStacking] = [
            PiStacking(
                receptor_ring_atom_indices=tuple(int(i) for i in rec_idx),
                ligand_ring_atom_indices=tuple(int(i) for i in lig_idx),
                kind=int(kind),
            )
            for rec_idx, lig_idx, kind in cm.find_stacking_interactions()
        ]

        result = cls(
            receptor=receptor,
            ligand=ligand,
            hbonds=hbonds,
            salt_bridges=salt_bridges,
            pi_cation=pi_cation,
            pi_stacking=pi_stacking,
        )

        if rotate:
            import biotite.structure as bio_struc

            from portein.rotate import get_best_transformation

            coords = result.atom_coords
            if len(coords) >= 2:
                # Combine then rotate then split to guarantee receptor and
                # ligand share exactly the same transformation.
                n_receptor = len(receptor)
                combined = receptor + ligand
                angles, translation = get_best_transformation(coords)
                combined = bio_struc.rotate(bio_struc.translate(combined, translation), angles)
                result.receptor = combined[:n_receptor]
                result.ligand = combined[n_receptor:]

        return result

    @property
    def atom_coords(self) -> np.ndarray:
        """Coordinates of every atom participating in any interaction.

        Pass to :func:`portein.get_best_transformation` to rotate a
        structure so the interactions occupy the largest 2D area of the
        projection, rather than the protein as a whole. (Or simply pass
        ``rotate=True`` to :meth:`InteractionSet.find`, which uses this
        internally.)

        Returns
        -------
        np.ndarray, shape ``(n_atoms, 3)``, dtype float64
            Stacked coordinates of donor and acceptor atoms (H-bonds),
            the charged atoms of salt bridges, the cation and ring atoms
            of π-cation interactions, and both ring sets of π-stacking
            interactions. Atoms appearing in multiple interactions are
            kept (they only weight the optimization slightly).
        """
        coords: list[np.ndarray] = []
        for hb in self.hbonds:
            donor_array = self.receptor if hb.donor_is_receptor else self.ligand
            acceptor_array = self.ligand if hb.donor_is_receptor else self.receptor
            coords.append(donor_array.coord[hb.donor_atom_idx])
            coords.append(acceptor_array.coord[hb.acceptor_atom_idx])
        for sb in self.salt_bridges:
            coords.append(self.receptor.coord[sb.receptor_atom_idx])
            coords.append(self.ligand.coord[sb.ligand_atom_idx])
        for pc in self.pi_cation:
            ring_array = self.ligand if pc.cation_in_receptor else self.receptor
            cation_array = self.receptor if pc.cation_in_receptor else self.ligand
            coords.append(cation_array.coord[pc.cation_atom_idx])
            for idx in pc.ring_atom_indices:
                coords.append(ring_array.coord[idx])
        for ps in self.pi_stacking:
            for idx in ps.receptor_ring_atom_indices:
                coords.append(self.receptor.coord[idx])
            for idx in ps.ligand_ring_atom_indices:
                coords.append(self.ligand.coord[idx])
        if not coords:
            return np.empty((0, 3), dtype=np.float64)
        return np.asarray(coords, dtype=np.float64)
