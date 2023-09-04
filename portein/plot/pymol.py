from pathlib import Path
from matplotlib import colors as m_colors
from dataclasses import dataclass
from PIL import Image
from portein import config
from portein.plot import image_utils
from pymol import cmd, util

@dataclass
class Pymol:
    protein_config: config.ProteinConfig
    pymol_config: config.PymolConfig

    def run(self, remove_intermediate_files=True):
        cmd.reinitialize()
        cmd.load(self.protein_config.pdb_file)
        self.set_pymol_settings()
        self.draw_representations()
        self.layer_images(remove_intermediate_files)
        cmd.reinitialize()

    def set_pymol_settings(self):
        cmd.bg_color("white")
        cmd.remove("solvent")
        cmd.set("max_threads", 8)
        if self.pymol_config.buffer is not None:
            cmd.zoom(buffer=self.pymol_config.buffer)
        for chain, color in self.protein_config.chain_to_color.items():
            cmd.color(f"0x{m_colors.to_hex(color).lower()[1:]}", f"chain {chain}")

    def draw_representations(self):
        for r, representation_config in enumerate(self.pymol_config.layers):
            self.set_representation_config(representation_config)
            cmd.hide("everything")
            self.show_representation(representation_config)
            cmd.ray(self.protein_config.width, self.protein_config.height)
            cmd.png(f"{self.protein_config.output_prefix}_{r}_pymol.png", 
                    width=self.protein_config.width, 
                    height=self.protein_config.height, 
                    dpi=representation_config.dpi)

    def set_representation_config(self, representation_config: config.PymolRepresentationConfig):
        for key, value in representation_config.pymol_settings.items():
            if key in ["specular", "depth_cue"]:
                if value:
                    cmd.set(key)
                else:
                    cmd.unset(key)
            else:
                cmd.set(key, value)

    def show_representation(self, representation_config: config.PymolRepresentationConfig):
        if representation_config.selection == "highlight":
            self.show_highlight_representation(representation_config)
        else:
            cmd.show(representation_config.representation, representation_config.selection)
        if representation_config.color is not None:
            cmd.color(representation_config.color, representation_config.selection)
        if representation_config.representation == "sticks":
            util.cnc("rep sticks")

    def show_highlight_representation(self, representation_config: config.PymolRepresentationConfig):
        for chain, colors in self.protein_config.highlight_residues.items():
            for color, residues in colors.items():
                if color is None:
                    color = self.protein_config.chain_to_color[chain]
                if representation_config.color is not None:
                    color = representation_config.color
                color = m_colors.to_hex(color).lower()[1:]
                cmd.select(f"highlight_{chain}_{color}", f"chain {chain} and resi {'+'.join(str(x) for x in residues)}")
                cmd.show(representation_config.representation, f"highlight_{chain}_{color}")
                cmd.color(f"0x{color}", f"highlight_{chain}_{color}")

    def layer_images(self, remove_intermediate_files=True):
        pngs = []
        for r, representation_config in enumerate(self.pymol_config.layers):
            img = Image.open(f"{self.protein_config.output_prefix}_{r}_pymol.png")
            pngs.append(image_utils.put_alpha(img, representation_config.transparency))
            if remove_intermediate_files:
                Path(f"{self.protein_config.output_prefix}_{r}_pymol.png").unlink() 
        image = Image.new("RGBA", pngs[0].size)
        for png in pngs:
            image = Image.alpha_composite(image, png)
        image.save(f"{self.protein_config.output_prefix}_pymol.png") 