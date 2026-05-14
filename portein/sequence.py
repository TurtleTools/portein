"""Sequence-alignment plotting helpers built on biotite's SymbolPlotter."""

from __future__ import annotations

from collections.abc import Callable

from biotite.sequence.align import Alignment
from biotite.sequence.graphics import SymbolPlotter
from matplotlib.patches import Rectangle
from matplotlib.typing import ColorType

ColorFn = Callable[[str, int, int, int, Alignment], ColorType]
"""Signature ``(symbol, seq_index, column_index, sequence_position, alignment) -> color``.

* ``symbol`` — the residue letter at this cell, or the plotter's
  ``gap_symbol`` for gap cells.
* ``seq_index`` — which row in the alignment (0..n_sequences-1).
* ``column_index`` — the *aligned* column (0..n_columns-1, includes gaps).
  Same value across sequences in that column.
* ``sequence_position`` — position in the original ungapped sequence
  ``alignment.sequences[seq_index]``, or ``-1`` if this cell is a gap.
* ``alignment`` — the full ``Alignment`` object for advanced access.

Ignore positional args you don't need: ``lambda s, *_: ...``.
"""


def _case_color_fn(symbol: str, seq_i: int, column_i: int, seq_pos: int, alignment: Alignment) -> ColorType:
    if seq_pos == -1:
        return "lightgray"
    return "salmon" if symbol.islower() else "white"


class AlignmentPlotter(SymbolPlotter):
    """Color alignment cells via a user-supplied callable.

    Plug into ``biotite.sequence.graphics.plot_alignment`` via its
    ``symbol_plotter`` argument. The default ``color_fn`` colors cells by
    character case (upper -> white, lower -> salmon) — a convention used by
    some MSA tools to mark unreliable regions — and renders gap cells in
    ``"lightgray"``. Pass any callable matching :data:`ColorFn` to implement
    arbitrary schemes.

    Examples
    --------
    Color each cell by a precomputed ``(n_seqs, n_columns)`` array::

        plotter = AlignmentPlotter(
            ax, color_fn=lambda s, si, ci, sp, a: cell_colors[si, ci]
        )

    Color residues by chemistry, ignoring everything else::

        scheme = {"A": "lightblue", "C": "yellow", ...}
        plotter = AlignmentPlotter(
            ax, color_fn=lambda s, *_: scheme.get(s.upper(), "white")
        )

    Color by per-residue pLDDT (sequence position, not aligned column)::

        plotter = AlignmentPlotter(
            ax, color_fn=lambda s, si, ci, sp, a: "gray" if sp == -1 else plddt_to_color(plddt[si][sp])
        )

    Color by conservation (per aligned column)::

        plotter = AlignmentPlotter(
            ax, color_fn=lambda s, si, ci, sp, a: column_colors[ci]
        )
    """

    def __init__(
        self,
        axes,
        *,
        color_fn: ColorFn | None = None,
        text_color: ColorType = "black",
        gap_symbol: str = "-",
        font_size: float | None = None,
    ) -> None:
        super().__init__(axes)
        self.color_fn: ColorFn = color_fn if color_fn is not None else _case_color_fn
        self.text_color = text_color
        self.gap_symbol = gap_symbol
        self.font_size = font_size

    def plot_symbol(self, bbox, alignment, column_i, seq_i) -> None:
        seq_pos = int(alignment.trace[column_i, seq_i])
        if seq_pos == -1:
            symbol = self.gap_symbol
        else:
            symbol = str(alignment.sequences[seq_i][seq_pos])
        bg = self.color_fn(symbol, seq_i, column_i, seq_pos, alignment)

        self.axes.add_patch(
            Rectangle(
                (bbox.x0, bbox.y0),
                bbox.width,
                bbox.height,
                color=bg,
                zorder=0,
                linewidth=0,
            )
        )
        self.axes.text(
            bbox.x0 + bbox.width / 2,
            bbox.y0 + bbox.height / 2,
            symbol,
            color=self.text_color,
            ha="center",
            va="center",
            fontsize=self.font_size,
            zorder=1,
        )
