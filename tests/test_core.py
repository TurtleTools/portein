"""Core Python API: structure I/O, rotation, ProteinConfig."""

import biotite.structure.io as bio_io
import numpy as np

import portein


def test_read_pdb_and_cif_roundtrip(example_pdb_path, tmp_path):
    pdb = portein.read_structure(str(example_pdb_path))
    assert len(pdb) > 0 and "CA" in set(pdb.atom_name)

    cif_path = tmp_path / "tmp.cif"
    bio_io.save_structure(str(cif_path), pdb)
    assert len(portein.read_structure(str(cif_path))) == len(pdb)


def test_rotate_protein_is_rigid(example_pdb_path):
    portein.compile_numba_functions()
    pdb = portein.read_structure(str(example_pdb_path))
    rotated = portein.rotate_protein(pdb)
    assert len(rotated) == len(pdb)
    ca = pdb[pdb.atom_name == "CA"].coord
    ca_rot = rotated[rotated.atom_name == "CA"].coord
    np.testing.assert_allclose(
        np.linalg.norm(ca[0] - ca[-1]),
        np.linalg.norm(ca_rot[0] - ca_rot[-1]),
        rtol=1e-4,
    )


def test_protein_config_from_real_pdb(example_pdb_path):
    cfg = portein.ProteinConfig(pdb_file=str(example_pdb_path), width=500)
    assert cfg.width == 500
    assert cfg.height and cfg.height > 0
    assert isinstance(cfg.chain_to_color, dict) and cfg.chain_to_color
