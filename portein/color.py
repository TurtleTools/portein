"""Helpers for mapping per-residue values to colors usable by ProteinConfig."""

from __future__ import annotations

from collections import defaultdict

from biotite.structure import AtomArray

PLDDT_RANGES: list[tuple[float, float, str]] = [
    (0.0, 50.0, "#FE7D45"),
    (50.0, 70.0, "#FFDB13"),
    (70.0, 90.0, "#64CBF3"),
    (90.0, 100.0, "#0053D6"),
]
"""Standard AlphaFold pLDDT bins: very low / low / confident / very confident."""


def by_plddt(
    structure: AtomArray,
    ranges: list[tuple[float, float, str]] = PLDDT_RANGES,
    attr: str = "b_factor",
) -> dict[str, dict[str, list[int]]]:
    """Group residues by pLDDT bin.

    Inspects the given per-atom attribute on CA atoms (default ``b_factor``,
    which AlphaFold writes the pLDDT into) and bins each residue into one
    of the supplied color bands. Values outside every band fall into the
    last band's color.

    Parameters
    ----------
    structure
        AtomArray with a per-atom score attribute. Only CA atoms are read.
    ranges
        List of ``(low, high, color)`` tuples. A value ``v`` is placed into
        the band where ``low <= v < high``. Bands are checked in order.
    attr
        Name of the per-atom attribute holding the score.

    Returns
    -------
    Nested dict shaped like ``{chain_id: {color: [residue_ids]}}``, suitable
    to pass straight to :class:`~portein.ProteinConfig.highlight_residues`.
    """
    ca = structure[structure.atom_name == "CA"]
    scores = getattr(ca, attr)

    fallback_color = ranges[-1][2]
    binned: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    for chain_id, res_id, score in zip(ca.chain_id, ca.res_id, scores, strict=True):
        for low, high, color in ranges:
            if low <= score < high:
                binned[chain_id][color].append(int(res_id))
                break
        else:
            binned[chain_id][fallback_color].append(int(res_id))
    return {chain: dict(by_color) for chain, by_color in binned.items()}
