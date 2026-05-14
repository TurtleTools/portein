"""Helpers for comparing and superimposing related structures."""

from __future__ import annotations

import numpy as np
from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment, SubstitutionMatrix, align_optimal
from biotite.structure import AtomArray, superimpose


def _structure_sequence(ca_atoms: AtomArray) -> ProteinSequence:
    """Build a ProteinSequence from the CA residue names of ``ca_atoms``."""
    letters: list[str] = []
    for res_name in ca_atoms.res_name:
        try:
            letters.append(ProteinSequence.convert_letter_3to1(res_name))
        except KeyError:
            letters.append("X")
    return ProteinSequence("".join(letters))


def _input_to_structure_map(
    ca_atoms: AtomArray,
    input_seq: ProteinSequence,
    matrix: SubstitutionMatrix,
) -> np.ndarray:
    """Map each position in ``input_seq`` to the corresponding CA index in
    ``ca_atoms``, or -1 if that input position has no structural counterpart.

    Works by aligning the CA-derived structure sequence against ``input_seq``
    with ``align_optimal`` and inverting the trace.
    """
    struct_seq = _structure_sequence(ca_atoms)
    realignment = align_optimal(struct_seq, input_seq, matrix=matrix)[0]
    mapping = np.full(len(input_seq), -1, dtype=np.int64)
    for struct_pos, input_pos in realignment.trace:
        if struct_pos != -1 and input_pos != -1:
            mapping[input_pos] = struct_pos
    return mapping


def _rename_chains(structure: AtomArray, chain_id: str | None) -> AtomArray:
    if chain_id is None:
        return structure
    renamed = structure.copy()
    renamed.chain_id = np.array([chain_id] * len(renamed))
    return renamed


def superimpose_by_alignment(
    query: AtomArray,
    target: AtomArray,
    alignment: Alignment,
    *,
    query_index: int = 0,
    target_index: int = 1,
    query_chain: str | None = None,
    target_chain: str | None = None,
    query_chain_id: str | None = "A",
    target_chain_id: str | None = "B",
    combine: bool = True,
    matrix: SubstitutionMatrix | None = None,
) -> tuple:
    """Superimpose ``target`` onto ``query`` using paired CA residues from an alignment.

    The supplied ``alignment`` is assumed to be between *input* sequences
    (typically SEQRES). Real structures often have missing residues, so the
    alignment positions don't map directly to CA indices. This helper:

    1. Extracts the CA-derived sequence from each structure (optionally
       restricted to one chain).
    2. Re-aligns each structure sequence to its input sequence with
       ``align_optimal``.
    3. Projects the input alignment trace through both re-alignments to
       obtain matched CA indices.
    4. Runs ``biotite.structure.superimpose`` on the matched CAs and applies
       the transformation to ``target``.

    Which sequence in ``alignment`` corresponds to ``query`` vs ``target`` is
    specified by ``query_index`` and ``target_index``. ``query`` is the fixed
    reference; ``target`` is moved into query's frame.

    Parameters
    ----------
    query
        Fixed reference structure. Not transformed.
    target
        Mobile structure. Transformed into ``query``'s frame.
    alignment
        biotite ``Alignment`` over the input sequences.
    query_index, target_index
        Indices into ``alignment.sequences`` / columns of ``alignment.trace``
        identifying which is which. Defaults match the natural reading
        ``superimpose_by_alignment(query, target, ...)`` with an alignment
        built as ``align_optimal(query_seq, target_seq, ...)``. For
        cigar-built alignments where the reference is sequence 0, pass
        ``query_index=1, target_index=0``.
    query_chain, target_chain
        If the structure spans multiple chains, restrict to this chain id
        when building the structure sequence and for the returned atoms.
        Pass ``None`` (default) to use all atoms.
    query_chain_id, target_chain_id
        Chain ids assigned to query and target atoms after superposition,
        so the two structures can be distinguished in downstream plotting.
        Pass ``None`` to leave a structure's chain ids unchanged.
    combine
        If True (default), return a single combined AtomArray of
        ``target_transformed + query_renamed`` plus the transformation and
        paired indices (3-tuple). If False, return ``target_transformed``
        and ``query_renamed`` separately (4-tuple) so the caller can load
        them as separate objects.
    matrix
        Substitution matrix for the structure-vs-input re-alignment.
        Defaults to BLOSUM62.

    Returns
    -------
    If ``combine`` is True
        ``(combined, transformation, paired_indices)``
    If ``combine`` is False
        ``(target_transformed, query_renamed, transformation, paired_indices)``

    ``paired_indices`` is an ``(n_pairs, 2)`` array of
    ``[query_ca_index, target_ca_index]`` into the chain-filtered CA-only
    sub-arrays of query and target.
    """
    query_for_seq = query if query_chain is None else query[query.chain_id == query_chain]
    target_for_seq = target if target_chain is None else target[target.chain_id == target_chain]
    query_ca = query_for_seq[query_for_seq.atom_name == "CA"]
    target_ca = target_for_seq[target_for_seq.atom_name == "CA"]

    if matrix is None:
        matrix = SubstitutionMatrix.std_protein_matrix()

    query_map = _input_to_structure_map(query_ca, alignment.sequences[query_index], matrix)
    target_map = _input_to_structure_map(target_ca, alignment.sequences[target_index], matrix)

    query_positions = alignment.trace[:, query_index]
    target_positions = alignment.trace[:, target_index]

    paired_list: list[tuple[int, int]] = []
    for qp, tp in zip(query_positions, target_positions, strict=True):
        if qp == -1 or tp == -1:
            continue
        qci = int(query_map[qp])
        tci = int(target_map[tp])
        if qci == -1 or tci == -1:
            continue
        paired_list.append((qci, tci))

    if not paired_list:
        raise ValueError("No CA pairs survived alignment-to-structure projection — cannot superimpose.")

    paired = np.array(paired_list, dtype=np.int64)
    query_indices = paired[:, 0]
    target_indices = paired[:, 1]

    _, transformation = superimpose(query_ca[query_indices], target_ca[target_indices])
    transformed_target = transformation.apply(target)

    query_renamed = _rename_chains(query, query_chain_id)
    target_renamed = _rename_chains(transformed_target, target_chain_id)

    if combine:
        combined = target_renamed + query_renamed
        return combined, transformation, paired
    return target_renamed, query_renamed, transformation, paired
