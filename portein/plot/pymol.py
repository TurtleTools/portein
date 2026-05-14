from dataclasses import dataclass
from pathlib import Path

import yaml
from matplotlib import colors as m_colors
from PIL import Image

from portein import config
from portein.interactions import InteractionSet
from portein.plot import image_utils

try:
    from pymol import cmd, util
except ImportError as _pymol_import_error:
    cmd = None
    util = None
    _PYMOL_IMPORT_ERROR = _pymol_import_error
else:
    _PYMOL_IMPORT_ERROR = None


def _require_pymol():
    if _PYMOL_IMPORT_ERROR is not None:
        raise ImportError(
            "PyMOL is required for portein.Pymol but is not installed. "
            "Install it with:\n"
            "    mamba install -c conda-forge pymol-open-source\n"
            "(or `conda install -c conda-forge pymol-open-source`)"
        ) from _PYMOL_IMPORT_ERROR


@dataclass
class Pymol:
    protein: config.ProteinConfig
    layers: list[config.PymolConfig | list[config.PymolConfig]]
    """Render groups, in order, alpha-composited in PIL.

    Each entry is either:

    * a single :class:`PymolConfig` — its own ``cmd.ray()`` pass; the
      layer's ``transparency`` becomes a flat 2D PIL alpha applied to
      the entire resulting PNG. Use for elements that should sit
      uniformly in front of or behind everything else (e.g., a
      translucent cartoon backdrop).
    * a nested ``list[PymolConfig]`` — all included layers are set up
      in one PyMOL scene and rendered with a single ``cmd.ray()`` pass.
      PyMOL's z-buffer interleaves their geometries per pixel, and each
      layer's ``transparency`` becomes a per-selection PyMOL setting
      (``cartoon_transparency`` / ``transparency`` / ``stick_transparency``).
      Use when overlapping geometries should depth-interleave correctly
      (e.g., ligand sticks weaving through binding-site residue sticks).
    """
    buffer: float = None
    """buffer around the molecule in Angstroms"""
    zoom_to: str = None
    """zoom to selection"""
    center_to: str = None
    """center to selection"""
    interactions: InteractionSet | None = None
    """Optional :class:`portein.interactions.InteractionSet` to overlay on
    the render. Each detected interaction is drawn as a PyMOL distance
    object (dashed line). Drawing happens at scene setup so the lines
    appear in every subsequent ray-traced layer."""

    @staticmethod
    def from_yaml(protein_config: Path, pymol_config: Path, buffer: float = None):
        protein = config.ProteinConfig.from_yaml(protein_config)
        layers: list[config.PymolConfig | list[config.PymolConfig]]
        if pymol_config is None:
            layers = [config.PymolConfig()]
        else:
            loader = yaml.SafeLoader
            loader.add_constructor("!PymolConfig", config.PymolConfig._from_yaml)
            with open(pymol_config) as f:
                layers = yaml.load(f, Loader=loader)
        return Pymol(protein, layers, buffer)

    def run(self, remove_intermediate_files=True):
        _require_pymol()
        cmd.reinitialize()
        cmd.load(self.protein.pdb_file)
        self.set_pymol_settings()
        groups = self._build_groups()
        pngs_with_alpha: list[tuple[str, float]] = []
        for index, group in enumerate(groups):
            png_path = self._render_group(group, index)
            # Only single-layer groups apply transparency as flat PIL alpha
            # on the resulting PNG; multi-layer groups handle transparency
            # per-selection inside the single ray-trace.
            transparency = group[0].transparency if len(group) == 1 else 0.0
            pngs_with_alpha.append((png_path, transparency))
        image_file = self._composite_groups(pngs_with_alpha, remove_intermediate_files)
        cmd.reinitialize()
        return image_file

    def _build_groups(self) -> list[list[config.PymolConfig]]:
        """Normalize ``self.layers`` (mixed PymolConfig + nested lists) into
        a flat list of render groups. Each group is a list of layers that
        share one ``cmd.ray()`` pass.
        """
        groups: list[list[config.PymolConfig]] = []
        for entry in self.layers:
            if isinstance(entry, list):
                groups.append(list(entry))
            else:
                groups.append([entry])
        return groups

    _REP_TRANSPARENCY_SETTING = {
        "cartoon": "cartoon_transparency",
        "surface": "transparency",
        "sticks": "stick_transparency",
        "spheres": "sphere_transparency",
        "ribbon": "ribbon_transparency",
    }

    def _render_group(self, group: list[config.PymolConfig], index: int) -> str:
        """Render one group of layers in a single ``cmd.ray()`` pass.

        Per-layer ``transparency`` in a multi-layer group becomes a PyMOL
        per-selection setting (``cartoon_transparency``, ``transparency``
        for surface, ``stick_transparency``, etc.) so transparencies are
        depth-aware. In a single-layer group, ``transparency`` is left to
        the PIL alpha step (see ``_composite_groups``).
        """
        cmd.hide("everything")
        multi = len(group) > 1
        for layer in group:
            self.change_settings(layer.pymol_settings, selection=layer.selection)
            self.show(layer=layer)
            if multi and layer.selection != "highlight":
                setting = self._REP_TRANSPARENCY_SETTING.get(layer.representation)
                if setting is not None:
                    cmd.set(setting, layer.transparency, layer.selection)
        # cmd.hide("everything") above also hid the distance objects created
        # by _draw_interactions during scene setup. Re-show them so the
        # interaction overlay appears in every group's render — the user
        # is responsible for rendering the participating residues
        # (add a layer with selection="portein_interaction_residues").
        if self.interactions is not None:
            cmd.show("dashes", "portein_*")
        dpi = group[0].dpi if group else 300
        cmd.ray(self.protein.width, self.protein.height)
        png_file = f"{self.protein.output_prefix}_{index}_pymol.png"
        cmd.png(png_file, width=self.protein.width, height=self.protein.height, dpi=dpi)
        return png_file

    def _composite_groups(self, pngs_with_alpha: list[tuple[str, float]], remove_intermediate: bool) -> str:
        """Alpha-composite one PNG per group into the final image."""
        composited = []
        for png_path, transparency in pngs_with_alpha:
            img = Image.open(png_path)
            if transparency > 0:
                img = image_utils.put_alpha(img, transparency)
            composited.append(img)
            if remove_intermediate:
                Path(png_path).unlink()
        image = Image.new("RGBA", composited[0].size)
        for img in composited:
            image = Image.alpha_composite(image, img)
        image_file = f"{self.protein.output_prefix}_pymol.png"
        image.save(image_file)
        return image_file

    def set_pymol_settings(self):
        cmd.bg_color("white")
        cmd.remove("solvent")
        cmd.set("max_threads", 8)
        # When sticks are shown alongside cartoon for the same residues,
        # this setting routes the cartoon ribbon around the side-chain
        # atoms instead of through them — the sticks no longer overlap
        # awkwardly with the backbone trace.
        cmd.set("cartoon_side_chain_helper", 1)
        if self.center_to is not None:
            cmd.center(self.center_to)
        if self.zoom_to is not None:
            cmd.zoom(self.zoom_to, buffer=self.buffer)
        elif self.buffer is not None:
            cmd.zoom(buffer=self.buffer)
        for chain, color in self.protein.chain_to_color.items():
            cmd.color(f"0x{m_colors.to_hex(color).lower()[1:]}", f"chain {chain}")
        if self.protein.chain_transparency:
            for chain, alpha in self.protein.chain_transparency.items():
                cmd.set("cartoon_transparency", alpha, f"chain {chain}")
                cmd.set("transparency", alpha, f"chain {chain}")
                cmd.set("stick_transparency", alpha, f"chain {chain}")
        # Apply per-residue highlight colors at scene setup so they survive
        # subsequent layer renders. Without this, only layers using
        # selection="highlight" (which calls show_highlight) carry the custom
        # colors — a single-pass render that shows multiple chains then
        # can't render highlight residues with their assigned colors and
        # have the structures depth-interleave correctly.
        self._apply_highlight_colors()
        if self.interactions is not None:
            self._draw_interactions()

    @staticmethod
    def _atom_pymol_selection(atom_array, atom_idx: int) -> str:
        """Build a PyMOL selection string targeting a single atom of an AtomArray.

        Uses chain_id + res_id + atom_name, which are stable across PDB load
        and unique within a structure (heavy atoms only). Wrapped in
        parentheses so it composes safely in larger expressions.
        """
        atom = atom_array[atom_idx]
        return f"(chain {atom.chain_id} and resi {int(atom.res_id)} and name {atom.atom_name})"

    @classmethod
    def _atom_set_pymol_selection(cls, atom_array, indices) -> str:
        """Build a PyMOL selection string targeting the union of multiple atoms."""
        return " or ".join(cls._atom_pymol_selection(atom_array, int(i)) for i in indices)

    def _draw_interactions(self) -> None:
        """Overlay :class:`InteractionSet` items as PyMOL distance objects.

        Color and dash conventions follow PLIP: blue for H-bonds, yellow
        for salt bridges, orange for π-cation, green/smudge for parallel/
        T-shaped π-stacking. Pseudoatoms are placed at ring centroids so
        π interactions are drawn between centroids/cations rather than
        between individual ring atoms.
        """
        inter = self.interactions
        if inter is None:
            return

        # Dashes need to survive ray-tracing visibly — default radius is too
        # thin at high resolution.
        cmd.set("dash_radius", 0.12)
        cmd.set("dash_round_ends", 1)

        # Hydrogen bonds — donor to acceptor, blue
        hbond_names: list[str] = []
        for i, hb in enumerate(inter.hbonds):
            donor_array = inter.receptor if hb.donor_is_receptor else inter.ligand
            acceptor_array = inter.ligand if hb.donor_is_receptor else inter.receptor
            donor_sel = self._atom_pymol_selection(donor_array, hb.donor_atom_idx)
            acceptor_sel = self._atom_pymol_selection(acceptor_array, hb.acceptor_atom_idx)
            name = f"portein_hbond_{i}"
            cmd.distance(name, donor_sel, acceptor_sel)
            hbond_names.append(name)
        for name in hbond_names:
            cmd.set("dash_color", "blue", name)
            cmd.hide("labels", name)

        # Salt bridges — yellow with wider gap
        sb_names: list[str] = []
        for i, sb in enumerate(inter.salt_bridges):
            rec_sel = self._atom_pymol_selection(inter.receptor, sb.receptor_atom_idx)
            lig_sel = self._atom_pymol_selection(inter.ligand, sb.ligand_atom_idx)
            name = f"portein_saltbridge_{i}"
            cmd.distance(name, rec_sel, lig_sel)
            sb_names.append(name)
        for name in sb_names:
            cmd.set("dash_color", "yellow", name)
            cmd.set("dash_gap", 0.5, name)
            cmd.hide("labels", name)

        # π-cation — ring centroid pseudoatom to cation, orange.
        # The pseudoatom is placed by PyMOL using the *loaded* structure's
        # coordinates (which may have been rotated by ProteinConfig), not
        # the InteractionSet's AtomArrays — otherwise rotation would
        # desynchronize centroid positions from the rest of the scene.
        for i, pc in enumerate(inter.pi_cation):
            ring_array = inter.ligand if pc.cation_in_receptor else inter.receptor
            cation_array = inter.receptor if pc.cation_in_receptor else inter.ligand
            ring_sel = self._atom_set_pymol_selection(ring_array, pc.ring_atom_indices)
            ring_ps = f"portein_picat_ring_{i}"
            cmd.pseudoatom(ring_ps, selection=ring_sel, mode="rms")
            cation_sel = self._atom_pymol_selection(cation_array, pc.cation_atom_idx)
            name = f"portein_picat_{i}"
            cmd.distance(name, ring_ps, cation_sel)
            cmd.set("dash_color", "orange", name)
            cmd.set("dash_gap", 0.3, name)
            cmd.set("dash_length", 0.6, name)
            cmd.hide("labels", name)
            cmd.hide("everything", ring_ps)

        # π-stacking — centroid to centroid, green (parallel) / smudge (T-shaped)
        for i, ps in enumerate(inter.pi_stacking):
            rec_sel = self._atom_set_pymol_selection(inter.receptor, ps.receptor_ring_atom_indices)
            lig_sel = self._atom_set_pymol_selection(inter.ligand, ps.ligand_ring_atom_indices)
            rec_ps = f"portein_pistack_r_{i}"
            lig_ps = f"portein_pistack_l_{i}"
            cmd.pseudoatom(rec_ps, selection=rec_sel, mode="rms")
            cmd.pseudoatom(lig_ps, selection=lig_sel, mode="rms")
            name = f"portein_pistack_{i}"
            cmd.distance(name, rec_ps, lig_ps)
            # biotite PiStacking enum: 1 = PARALLEL, 2 = T_SHAPED
            color = "green" if ps.kind == 1 else "smudge"
            cmd.set("dash_color", color, name)
            cmd.set("dash_gap", 0.3, name)
            cmd.set("dash_length", 0.6, name)
            cmd.hide("labels", name)
            cmd.hide("everything", rec_ps)
            cmd.hide("everything", lig_ps)

        # Build a named selection of every receptor residue that participates
        # in any interaction so we can show them as sticks alongside the dashes.
        # Without this, the dashes terminate inside invisible cartoon residues
        # and look like they go nowhere.
        residues_by_chain: dict[str, set[int]] = {}

        def _record(atom_array, atom_idx: int) -> None:
            atom = atom_array[atom_idx]
            residues_by_chain.setdefault(str(atom.chain_id), set()).add(int(atom.res_id))

        for hb in inter.hbonds:
            if hb.donor_is_receptor:
                _record(inter.receptor, hb.donor_atom_idx)
            else:
                _record(inter.receptor, hb.acceptor_atom_idx)
        for sb in inter.salt_bridges:
            _record(inter.receptor, sb.receptor_atom_idx)
        for pc in inter.pi_cation:
            if pc.cation_in_receptor:
                _record(inter.receptor, pc.cation_atom_idx)
            else:
                for idx in pc.ring_atom_indices:
                    _record(inter.receptor, idx)
        for ps in inter.pi_stacking:
            for idx in ps.receptor_ring_atom_indices:
                _record(inter.receptor, idx)

        if residues_by_chain:
            sel_parts = [
                f"(chain {ch} and resi {'+'.join(str(r) for r in sorted(rids))})"
                for ch, rids in residues_by_chain.items()
            ]
            cmd.select("portein_interaction_residues", " or ".join(sel_parts))

    def _apply_highlight_colors(self):
        """Create named selections for highlight residues and apply their colors.

        Does not show anything — visibility is handled per-layer.
        """
        for chain, colors in self.protein.highlight_residues.items():
            for color, residues in colors.items():
                effective_color = color if color is not None else self.protein.chain_to_color[chain]
                color_hex = m_colors.to_hex(effective_color).lower()[1:]
                sel_name = f"highlight_{chain}_{color_hex}"
                cmd.select(sel_name, f"chain {chain} and resi {'+'.join(str(x) for x in residues)}")
                cmd.color(f"0x{color_hex}", sel_name)

    # PyMOL settings that support per-atom scoping for representation styling.
    # When found in a layer's pymol_settings, they're applied with
    # cmd.set(name, value, layer.selection) so they don't bleed into later
    # layers' rendering of the same atoms with different selections.
    _PER_ATOM_SETTINGS = frozenset(
        {
            "stick_radius",
            "stick_ball_ratio",
            "cartoon_transparency",
            "transparency",
            "stick_transparency",
            "sphere_transparency",
            "ribbon_radius",
            "ribbon_width",
            "sphere_scale",
            "mesh_radius",
            "label_size",
            "label_color",
        }
    )

    @classmethod
    def change_settings(cls, pymol_settings, selection: str = "all"):
        for key, value in pymol_settings.items():
            if key in ["specular", "depth_cue"]:
                if value:
                    cmd.set(key)
                else:
                    cmd.unset(key)
            elif key in cls._PER_ATOM_SETTINGS:
                cmd.set(key, value, selection)
            else:
                cmd.set(key, value)

    def show(self, layer: config.PymolConfig):
        if layer.selection == "highlight":
            self.show_highlight(layer=layer)
        else:
            cmd.show(layer.representation, layer.selection)
        if layer.spectrum is not None:
            cmd.spectrum(
                expression=layer.spectrum,
                selection=layer.selection,
                palette=layer.color,
                minimum=0,
                maximum=100,
            )
        elif layer.color is not None:
            cmd.color(f"0x{m_colors.to_hex(layer.color).lower()[1:]}", layer.selection)
        if layer.representation == "sticks":
            util.cnc("rep sticks")

    def show_highlight(self, layer: config.PymolConfig):
        """Show all named highlight selections in this layer.

        Selection names and per-residue colors are pre-applied by
        ``_apply_highlight_colors`` during scene setup; this method only
        toggles visibility and applies any per-layer color override.
        """
        for chain, colors in self.protein.highlight_residues.items():
            for color, _ in colors.items():
                effective_color = color if color is not None else self.protein.chain_to_color[chain]
                color_hex = m_colors.to_hex(effective_color).lower()[1:]
                sel_name = f"highlight_{chain}_{color_hex}"
                cmd.show(layer.representation, sel_name)
                if layer.color is not None:
                    cmd.color(f"0x{m_colors.to_hex(layer.color).lower()[1:]}", sel_name)
