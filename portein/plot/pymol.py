from dataclasses import dataclass
from pathlib import Path

import yaml
from matplotlib import colors as m_colors
from PIL import Image

from portein import config
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
    layers: list[config.PymolConfig]
    """order of representations, selections, and transparencies to layer on top of each other"""
    buffer: float = None
    """buffer around the molecule in Angstroms"""
    zoom_to: str = None
    """zoom to selection"""
    center_to: str = None
    """center to selection"""
    combine: bool = False
    """If False (default), render the portein way: each layer is its own
    ``cmd.ray()`` pass, and the PNGs are alpha-composited in PIL. Layer
    transparency is uniform 2D alpha across each layer's image; layers
    cannot depth-interleave with each other.
    If True, render the PyMOL way: set up all layers in a single scene
    and ``cmd.ray()`` once. PyMOL's z-buffer interleaves overlapping
    geometries per-pixel. Per-layer ``transparency`` then becomes a
    PyMOL ``cartoon_transparency`` / ``transparency`` / ``stick_transparency``
    setting on each layer's selection (depth-aware) instead of a flat
    PIL alpha. Choose this when overlapping structures or representations
    should weave through each other rather than stack as flat decals."""

    @staticmethod
    def from_yaml(protein_config: Path, pymol_config: Path, buffer: float = None):
        protein = config.ProteinConfig.from_yaml(protein_config)
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
        if self.combine:
            image_file = self._draw_combined()
        else:
            for index, layer in enumerate(self.layers):
                self.draw(layer, index)
            image_file = self.layer_images(remove_intermediate_files)
        cmd.reinitialize()
        return image_file

    _REP_TRANSPARENCY_SETTING = {
        "cartoon": "cartoon_transparency",
        "surface": "transparency",
        "sticks": "stick_transparency",
        "spheres": "sphere_transparency",
        "ribbon": "ribbon_transparency",
    }

    def _draw_combined(self) -> str:
        """Set up all layers in one PyMOL scene and ray-trace once.

        Per-layer ``transparency`` becomes a PyMOL per-selection setting
        matching the layer's representation (``cartoon_transparency``,
        ``transparency`` for surface, ``stick_transparency``, etc.) so
        transparencies are depth-aware and don't bleed across representations.
        Later layers' settings on overlapping atoms override earlier ones —
        order the layers from most general to most specific.
        """
        cmd.hide("everything")
        for layer in self.layers:
            self.change_settings(layer.pymol_settings)
            self.show(layer=layer)
            if layer.selection != "highlight":
                setting = self._REP_TRANSPARENCY_SETTING.get(layer.representation)
                if setting is not None:
                    cmd.set(setting, layer.transparency, layer.selection)
        dpi = self.layers[0].dpi if self.layers else 300
        cmd.ray(self.protein.width, self.protein.height)
        image_file = f"{self.protein.output_prefix}_pymol.png"
        cmd.png(image_file, width=self.protein.width, height=self.protein.height, dpi=dpi)
        return image_file

    def set_pymol_settings(self):
        cmd.bg_color("white")
        cmd.remove("solvent")
        cmd.set("max_threads", 8)
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

    def draw(self, layer: config.PymolConfig, index: int):
        self.change_settings(layer.pymol_settings)
        cmd.hide("everything")
        self.show(layer=layer)
        cmd.ray(self.protein.width, self.protein.height)
        cmd.png(
            f"{self.protein.output_prefix}_{index}_pymol.png",
            width=self.protein.width,
            height=self.protein.height,
            dpi=layer.dpi,
        )

    @staticmethod
    def change_settings(pymol_settings):
        for key, value in pymol_settings.items():
            if key in ["specular", "depth_cue"]:
                if value:
                    cmd.set(key)
                else:
                    cmd.unset(key)
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

    def layer_images(self, remove_intermediate_files=True):
        pngs = []
        for index, layer in enumerate(self.layers):
            img = Image.open(f"{self.protein.output_prefix}_{index}_pymol.png")
            pngs.append(image_utils.put_alpha(img, layer.transparency))
            if remove_intermediate_files:
                Path(f"{self.protein.output_prefix}_{index}_pymol.png").unlink()
        image = Image.new("RGBA", pngs[0].size)
        for png in pngs:
            image = Image.alpha_composite(image, png)
        image_file = f"{self.protein.output_prefix}_pymol.png"
        image.save(image_file)
        return image_file
