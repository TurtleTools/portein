from pathlib import Path
from matplotlib import colors as m_colors
from dataclasses import dataclass
from PIL import Image
from portein import config
from portein.plot import image_utils
import typing
import yaml

try:
    from pymol import cmd, util
except ImportError:
    print("Warning: pymol not installed, cannot use Pymol class")
    pass


@dataclass
class Pymol:
    protein: config.ProteinConfig
    layers: typing.List[config.PymolConfig]
    """order of representations, selections, and transparencies to layer on top of each other"""
    buffer: float = None
    """buffer around the molecule in Angstroms"""
    zoom_to: str = None
    """zoom to selection"""
    center_to: str = None
    """center to selection"""

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
        cmd.reinitialize()
        cmd.load(self.protein.pdb_file)
        self.set_pymol_settings()
        for index, layer in enumerate(self.layers):
            self.draw(layer, index)
        image_file = self.layer_images(remove_intermediate_files)
        cmd.reinitialize()
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
        for chain, colors in self.protein.highlight_residues.items():
            for color, residues in colors.items():
                if color is None:
                    color = self.protein.chain_to_color[chain]
                if layer.color is not None:
                    color = layer.color
                color = m_colors.to_hex(color).lower()[1:]
                cmd.select(
                    f"highlight_{chain}_{color}",
                    f"chain {chain} and resi {'+'.join(str(x) for x in residues)}",
                )
                cmd.show(layer.representation, f"highlight_{chain}_{color}")
                cmd.color(f"0x{color}", f"highlight_{chain}_{color}")

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
