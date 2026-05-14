"""External renderers: PyMOL and DSSP integration (auto-skipped if missing)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

import portein  # noqa: E402


def test_pymol_missing_gives_actionable_error(monkeypatch):
    from portein.plot import pymol as pymol_mod

    monkeypatch.setattr(pymol_mod, "_PYMOL_IMPORT_ERROR", ImportError("simulated"))
    with pytest.raises(ImportError, match="pymol-open-source"):
        pymol_mod._require_pymol()


@pytest.mark.pymol
def test_pymol_renders_image(example_pdb_path, tmp_path):
    protein = portein.ProteinConfig(
        pdb_file=str(example_pdb_path),
        width=200,
        output_prefix=str(tmp_path / "img"),
    )
    image_file = portein.Pymol(
        protein=protein,
        layers=[portein.PymolConfig(representation="cartoon")],
    ).run()
    assert Path(image_file).exists() and Path(image_file).stat().st_size > 0


@pytest.mark.pymol
def test_pymol_combined_mode_renders_image(example_pdb_path, tmp_path):
    """combine=True takes a different render path — one ray-trace instead of N."""
    protein = portein.ProteinConfig(
        pdb_file=str(example_pdb_path),
        width=200,
        output_prefix=str(tmp_path / "combined"),
    )
    image_file = portein.Pymol(
        protein=protein,
        layers=[
            portein.PymolConfig(representation="cartoon", selection="all", transparency=0.3),
            portein.PymolConfig(representation="sticks", selection="resn GNP"),
        ],
        combine=True,
    ).run()
    assert Path(image_file).exists() and Path(image_file).stat().st_size > 0


@pytest.mark.dssp
def test_secondary_structure_runs(example_pdb_path, tmp_path):
    protein = portein.ProteinConfig(
        pdb_file=str(example_pdb_path),
        width=500,
        output_prefix=str(tmp_path / "ss"),
    )
    ss = portein.SecondaryStructure(protein_config=protein, dpi=72)
    ss.run()
    assert ss.coords is not None and len(ss.coords) > 0
    assert ss.ss_elements
    plt.close("all")
