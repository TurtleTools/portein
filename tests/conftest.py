"""Shared fixtures and skip markers for the portein test suite."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

EXAMPLE_PDB = Path(__file__).resolve().parent.parent / "examples" / "7lc2_rotated.pdb"


@pytest.fixture(scope="session")
def example_pdb_path() -> Path:
    if not EXAMPLE_PDB.exists():
        pytest.skip(f"Example PDB not found: {EXAMPLE_PDB}")
    return EXAMPLE_PDB


def pytest_collection_modifyitems(config, items):
    """Auto-skip `pymol` / `dssp` tests when those tools are missing."""
    try:
        import pymol  # noqa: F401

        have_pymol = True
    except ImportError:
        have_pymol = False
    have_dssp = shutil.which("mkdssp") is not None

    skip_pymol = pytest.mark.skip(reason="pymol-open-source not installed")
    skip_dssp = pytest.mark.skip(reason="mkdssp not installed")
    for item in items:
        if "pymol" in item.keywords and not have_pymol:
            item.add_marker(skip_pymol)
        if "dssp" in item.keywords and not have_dssp:
            item.add_marker(skip_dssp)
