"""Unit tests for the public helpers in portein.image / .color / .compare / .sequence."""

from __future__ import annotations

import biotite.structure as bio_struct
import matplotlib
import matplotlib.colors as m_colors
import numpy as np
from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from portein.color import PLDDT_RANGES, by_plddt  # noqa: E402
from portein.compare import superimpose_by_alignment  # noqa: E402
from portein.image import crop_to_content  # noqa: E402
from portein.sequence import AlignmentPlotter  # noqa: E402


def _make_ca_only_structure(plddt_values: list[float], chain_id: str = "A") -> bio_struct.AtomArray:
    n = len(plddt_values)
    arr = bio_struct.AtomArray(n)
    arr.coord = np.zeros((n, 3), dtype=np.float32)
    arr.atom_name = np.array(["CA"] * n)
    arr.chain_id = np.array([chain_id] * n)
    arr.res_id = np.arange(1, n + 1, dtype=np.int32)
    arr.res_name = np.array(["ALA"] * n)
    arr.element = np.array(["C"] * n)
    arr.set_annotation("b_factor", np.asarray(plddt_values, dtype=np.float32))
    return arr


def _make_padded_image(content: tuple[int, int], pad: int, bg: tuple[int, ...]) -> Image.Image:
    h, w = content
    canvas = np.tile(np.array(bg, dtype=np.uint8), (h + 2 * pad, w + 2 * pad, 1))
    canvas[pad : pad + h, pad : pad + w] = np.array([10, 20, 30, 255], dtype=np.uint8)
    return Image.fromarray(canvas)


def test_crop_to_content_strips_white_padding():
    img = _make_padded_image(content=(5, 7), pad=3, bg=(255, 255, 255, 255))
    cropped = crop_to_content(img)
    assert cropped.size == (7, 5)


def test_crop_to_content_strips_transparent_padding():
    img = _make_padded_image(content=(4, 6), pad=2, bg=(0, 0, 0, 0))
    cropped = crop_to_content(img)
    assert cropped.size == (6, 4)


def test_crop_to_content_returns_original_when_all_background():
    arr = np.full((10, 10, 4), [255, 255, 255, 255], dtype=np.uint8)
    img = Image.fromarray(arr)
    cropped = crop_to_content(img)
    assert cropped.size == img.size


def test_crop_to_content_accepts_numpy_array():
    img = _make_padded_image(content=(3, 3), pad=4, bg=(255, 255, 255, 255))
    cropped = crop_to_content(np.asarray(img))
    assert cropped.size == (3, 3)


def test_by_plddt_bins_residues_by_band():
    structure = _make_ca_only_structure([30.0, 60.0, 80.0, 95.0])
    out = by_plddt(structure)
    assert set(out.keys()) == {"A"}
    a = out["A"]
    assert a["#FE7D45"] == [1]
    assert a["#FFDB13"] == [2]
    assert a["#64CBF3"] == [3]
    assert a["#0053D6"] == [4]


def test_by_plddt_falls_back_to_last_color_when_out_of_range():
    structure = _make_ca_only_structure([150.0])
    out = by_plddt(structure)
    assert out["A"][PLDDT_RANGES[-1][2]] == [1]


def test_by_plddt_splits_chains():
    a = _make_ca_only_structure([95.0, 95.0], chain_id="A")
    b = _make_ca_only_structure([30.0], chain_id="B")
    out = by_plddt(a + b)
    assert out["A"]["#0053D6"] == [1, 2]
    assert out["B"]["#FE7D45"] == [1]


def _ca_with_coords(coords: np.ndarray, chain_id: str = "A") -> bio_struct.AtomArray:
    n = len(coords)
    arr = bio_struct.AtomArray(n)
    arr.coord = coords.astype(np.float32)
    arr.atom_name = np.array(["CA"] * n)
    arr.chain_id = np.array([chain_id] * n)
    arr.res_id = np.arange(1, n + 1, dtype=np.int32)
    arr.res_name = np.array(["ALA"] * n)
    arr.element = np.array(["C"] * n)
    return arr


def test_superimpose_by_alignment_recovers_known_transform():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(10, 3)).astype(np.float32)
    query = _ca_with_coords(coords, chain_id="A")

    angle = np.pi / 4
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ],
        dtype=np.float32,
    )
    translation = np.array([5.0, -2.0, 3.0], dtype=np.float32)
    target_coords = coords @ rotation.T + translation
    target = _ca_with_coords(target_coords, chain_id="C")

    seq = ProteinSequence("A" * 10)
    trace = np.column_stack([np.arange(10), np.arange(10)])
    alignment = Alignment(sequences=[seq, seq], trace=trace, score=0)

    combined, _, paired = superimpose_by_alignment(query, target, alignment)

    assert paired.shape == (10, 2)
    # First 10 atoms in combined are the transformed (and renamed) target.
    transformed_target = combined[:10]
    np.testing.assert_allclose(transformed_target.coord, query.coord, atol=1e-4)
    # Target gets renamed to chain "B" by default, query to chain "A".
    np.testing.assert_array_equal(combined[:10].chain_id, ["B"] * 10)
    np.testing.assert_array_equal(combined[10:].chain_id, ["A"] * 10)


def test_superimpose_by_alignment_can_return_separate():
    rng = np.random.default_rng(2)
    coords = rng.normal(size=(5, 3)).astype(np.float32)
    query = _ca_with_coords(coords, chain_id="X")
    target = _ca_with_coords(coords + 3.0, chain_id="Y")
    seq = ProteinSequence("A" * 5)
    trace = np.column_stack([np.arange(5), np.arange(5)])
    alignment = Alignment(sequences=[seq, seq], trace=trace, score=0)

    out = superimpose_by_alignment(query, target, alignment, combine=False)
    assert len(out) == 4
    target_t, query_r, _, _ = out
    np.testing.assert_array_equal(target_t.chain_id, ["B"] * 5)
    np.testing.assert_array_equal(query_r.chain_id, ["A"] * 5)
    np.testing.assert_allclose(target_t.coord, query.coord, atol=1e-4)


class _FakeBbox:
    def __init__(self, x0=0.0, y0=0.0, width=1.0, height=1.0):
        self.x0 = x0
        self.y0 = y0
        self.width = width
        self.height = height


def _make_alignment(seq_a: str, seq_b: str) -> Alignment:
    sa = ProteinSequence(seq_a.upper())
    sb = ProteinSequence(seq_b.upper())
    trace = np.column_stack([np.arange(len(seq_a)), np.arange(len(seq_b))])
    aln = Alignment(sequences=[sa, sb], trace=trace, score=0)
    # Preserve case info on the original strings so AlignmentPlotter can read it.
    aln.sequences[0] = seq_a
    aln.sequences[1] = seq_b
    return aln


def _patch_color_to_rgb(patch):
    return m_colors.to_rgb(patch.get_facecolor())


def test_alignment_plotter_case_based_coloring():
    fig, ax = plt.subplots()
    aln = _make_alignment("Aa", "Aa")
    plotter = AlignmentPlotter(ax)
    plotter.plot_symbol(_FakeBbox(0, 0), aln, column_i=0, seq_i=0)  # 'A' upper
    plotter.plot_symbol(_FakeBbox(1, 0), aln, column_i=1, seq_i=0)  # 'a' lower

    patches = ax.patches
    assert len(patches) == 2
    assert _patch_color_to_rgb(patches[0]) == m_colors.to_rgb("white")
    assert _patch_color_to_rgb(patches[1]) == m_colors.to_rgb("salmon")
    plt.close(fig)


def test_alignment_plotter_color_fn_by_cell_array():
    fig, ax = plt.subplots()
    aln = _make_alignment("AA", "AA")
    cell_colors = np.array([["red", "blue"], ["green", "yellow"]])
    plotter = AlignmentPlotter(ax, color_fn=lambda s, si, ci, sp, a: cell_colors[si, ci])
    plotter.plot_symbol(_FakeBbox(), aln, column_i=0, seq_i=0)
    plotter.plot_symbol(_FakeBbox(), aln, column_i=1, seq_i=1)

    assert _patch_color_to_rgb(ax.patches[0]) == m_colors.to_rgb("red")
    assert _patch_color_to_rgb(ax.patches[1]) == m_colors.to_rgb("yellow")
    plt.close(fig)


def test_alignment_plotter_color_fn_by_sequence_position():
    fig, ax = plt.subplots()
    sa = ProteinSequence("AAA")
    sb = ProteinSequence("AA")
    trace = np.array([[0, 0], [1, -1], [2, 1]])
    aln = Alignment(sequences=[sa, sb], trace=trace, score=0)
    # color sequence 1 cells by their sequence_position (gap -> "black")
    plotter = AlignmentPlotter(
        ax,
        color_fn=lambda s, si, ci, sp, a: "black" if sp == -1 else ["red", "green", "blue"][sp],
    )
    plotter.plot_symbol(_FakeBbox(), aln, column_i=0, seq_i=1)  # sp=0 -> red
    plotter.plot_symbol(_FakeBbox(), aln, column_i=1, seq_i=1)  # gap -> black
    plotter.plot_symbol(_FakeBbox(), aln, column_i=2, seq_i=1)  # sp=1 -> green

    assert _patch_color_to_rgb(ax.patches[0]) == m_colors.to_rgb("red")
    assert _patch_color_to_rgb(ax.patches[1]) == m_colors.to_rgb("black")
    assert _patch_color_to_rgb(ax.patches[2]) == m_colors.to_rgb("green")
    plt.close(fig)


def test_alignment_plotter_color_fn_by_symbol():
    fig, ax = plt.subplots()
    aln = _make_alignment("AG", "GA")
    scheme = {"A": "tomato", "G": "skyblue"}
    plotter = AlignmentPlotter(ax, color_fn=lambda s, *_: scheme[s.upper()])
    plotter.plot_symbol(_FakeBbox(), aln, column_i=0, seq_i=0)  # 'A'
    plotter.plot_symbol(_FakeBbox(), aln, column_i=1, seq_i=0)  # 'G'

    assert _patch_color_to_rgb(ax.patches[0]) == m_colors.to_rgb("tomato")
    assert _patch_color_to_rgb(ax.patches[1]) == m_colors.to_rgb("skyblue")
    plt.close(fig)


def test_alignment_plotter_default_handles_gaps():
    fig, ax = plt.subplots()
    sa = ProteinSequence("AA")
    sb = ProteinSequence("AA")
    trace = np.array([[0, -1], [1, 0]])
    aln = Alignment(sequences=[sa, sb], trace=trace, score=0)

    plotter = AlignmentPlotter(ax)
    plotter.plot_symbol(_FakeBbox(), aln, column_i=0, seq_i=1)  # gap row

    assert _patch_color_to_rgb(ax.patches[0]) == m_colors.to_rgb("lightgray")
    plt.close(fig)


def _ca_with_resnames(coords: np.ndarray, res_names: list[str], chain_id: str = "A") -> bio_struct.AtomArray:
    n = len(coords)
    arr = bio_struct.AtomArray(n)
    arr.coord = coords.astype(np.float32)
    arr.atom_name = np.array(["CA"] * n)
    arr.chain_id = np.array([chain_id] * n)
    arr.res_id = np.arange(1, n + 1, dtype=np.int32)
    arr.res_name = np.array(res_names)
    arr.element = np.array(["C"] * n)
    return arr


def test_superimpose_by_alignment_handles_structure_missing_residues():
    # Input alignment is between sequences of length 7 ("AGPAGNK") but the
    # structures only carry CA atoms for the first 5 ("AGPAG") — positions
    # 5 and 6 of the input sequence have no structural counterpart.
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(5, 3)).astype(np.float32)
    res_names = ["ALA", "GLY", "PRO", "ALA", "GLY"]
    query = _ca_with_resnames(coords, res_names, chain_id="A")

    angle = np.pi / 6
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ],
        dtype=np.float32,
    )
    target_coords = coords @ rotation.T + np.array([2.0, 1.0, -1.0])
    target = _ca_with_resnames(target_coords, res_names, chain_id="C")

    full_seq = ProteinSequence("AGPAGNK")
    trace = np.column_stack([np.arange(7), np.arange(7)])
    alignment = Alignment(sequences=[full_seq, full_seq], trace=trace, score=0)

    combined, _, paired = superimpose_by_alignment(query, target, alignment)

    # Only the first 5 positions exist in the structures; positions 5/6 drop out.
    assert paired.tolist() == [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]]
    transformed_target = combined[:5]
    np.testing.assert_allclose(transformed_target.coord, query.coord, atol=1e-4)


def test_superimpose_by_alignment_skips_gap_rows():
    coords = np.arange(15, dtype=np.float32).reshape(5, 3)
    query = _ca_with_coords(coords)
    target = _ca_with_coords(coords + 1.0, chain_id="C")
    seq = ProteinSequence("A" * 5)
    trace = np.array([[0, 0], [1, -1], [-1, 2], [3, 3], [4, 4]])
    alignment = Alignment(sequences=[seq, seq], trace=trace, score=0)

    _, _, paired = superimpose_by_alignment(query, target, alignment)
    # Only rows with no -1 should remain: (0,0), (3,3), (4,4)
    assert paired.tolist() == [[0, 0], [3, 3], [4, 4]]
