"""Tests for portein.interactions (peppr-backed PLI detection + PyMOL overlay)."""

from __future__ import annotations

from pathlib import Path

import biotite.structure as struc
import pytest

import portein

DATA_PDB = Path(__file__).resolve().parent.parent / "docs" / "_data" / "7lc2.pdb"


def _have_peppr() -> bool:
    try:
        import peppr  # noqa: F401
    except ImportError:
        return False
    return True


pytestmark = pytest.mark.skipif(not _have_peppr(), reason="peppr not installed")


@pytest.fixture(scope="module")
def kras_receptor_ligand():
    if not DATA_PDB.exists():
        pytest.skip(f"missing test data: {DATA_PDB}")
    pdb = portein.read_structure(str(DATA_PDB))
    chain_a = pdb[pdb.chain_id == "A"]
    chain_a.bonds = struc.connect_via_residue_names(chain_a)
    ligand = chain_a[chain_a.res_name == "GNP"]
    receptor = chain_a[chain_a.res_name != "GNP"]
    return receptor, ligand


def test_find_interactions_recovers_known_interactions(kras_receptor_ligand):
    """7lc2 chain A bound to GNP: expect H-bonds, π-cation, and π-stacking."""
    receptor, ligand = kras_receptor_ligand
    iset = portein.InteractionSet.find(receptor=receptor, ligand=ligand)

    assert iset.receptor is receptor
    assert iset.ligand is ligand
    # At minimum: many H-bonds in the well-bound nucleotide pocket
    assert len(iset.hbonds) > 5
    # GNP's guanine ring is famously π-stacked with the phenylalanine binding-pocket
    # residue and π-cation interactions exist in this complex.
    assert len(iset.pi_stacking) >= 1
    assert len(iset.pi_cation) >= 1


def test_hbond_objects_point_at_valid_atoms(kras_receptor_ligand):
    receptor, ligand = kras_receptor_ligand
    iset = portein.InteractionSet.find(receptor=receptor, ligand=ligand)
    for hb in iset.hbonds:
        donor_array = receptor if hb.donor_is_receptor else ligand
        acceptor_array = ligand if hb.donor_is_receptor else receptor
        # Indices are in bounds
        assert 0 <= hb.donor_atom_idx < len(donor_array)
        assert 0 <= hb.acceptor_atom_idx < len(acceptor_array)
        # Donor/acceptor are heavy atoms (not hydrogen)
        assert donor_array[hb.donor_atom_idx].element != "H"
        assert acceptor_array[hb.acceptor_atom_idx].element != "H"


@pytest.mark.pymol
def test_pymol_renders_with_interactions(kras_receptor_ligand, tmp_path):
    """Full pipeline: detect interactions → ProteinConfig → Pymol → ray."""
    receptor, ligand = kras_receptor_ligand
    iset = portein.InteractionSet.find(receptor=receptor, ligand=ligand)

    combined = receptor + ligand
    protein_config = portein.ProteinConfig(
        pdb_file=combined,
        rotate=True,
        width=300,
        chain_to_color={"A": "lightblue"},
        output_prefix=str(tmp_path / "interactions"),
    )
    pymol = portein.Pymol(
        protein=protein_config,
        layers=[
            [
                portein.PymolConfig(representation="cartoon", selection="all"),
                portein.PymolConfig(representation="sticks", selection="resn GNP", color="green"),
            ],
        ],
        interactions=iset,
    )
    image_file = pymol.run()
    assert Path(image_file).exists() and Path(image_file).stat().st_size > 0
