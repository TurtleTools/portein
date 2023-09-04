from dataclasses import dataclass
from matplotlib import colors as m_colors
from pathlib import Path
import subprocess
from PIL import Image
from portein import config
from portein.plot import image_utils


@dataclass
class Illustrate:
    protein_config: config.ProteinConfig
    illustrate_config: config.IllustrateConfig

    def make_command_file(self):
        with open(f"{self.protein_config.output_prefix}_illustrate.inp", "w") as f:
            self.write_omit_hydrogens_and_water(f)
            self.write_highlight_residues(f)
            self.write_chain_colors(f)
            self.write_end_and_transformations(f)

    def write_omit_hydrogens_and_water(self, file):
        file.write(f"""read
{self.protein_config.pdb_file}
HETATM-----HOH-- 0,9999, 0.5,0.5,0.5, 0.0
ATOM  -H-------- 0,9999, 0.5,0.5,0.5, 0.0
ATOM  H--------- 0,9999, 0.5,0.5,0.5, 0.0
""")

    def write_highlight_residues(self, file):
        for chain in self.protein_config.chain_to_color:
            if chain in self.protein_config.highlight_residues:
                self.write_residue_ranges(file, chain)

    def write_residue_ranges(self, file, chain):
        for color, residue_ranges in self.protein_config.chain_to_residue_range_color[chain].items():
            color = m_colors.to_rgb(color)
            for start, end in residue_ranges:
                self.write_residue_range(file, chain, start, end, color, color)

    def write_residue_range(self, file, chain, start, end, calpha_color, sidechain_color):
        file.write(f"ATOM  -C-------{chain} {start},{end+1}, {calpha_color[0]:.1f},{calpha_color[1]:.1f},{calpha_color[2]:.1f}, {self.illustrate_config.carbon_radius}\n")
        file.write(f"ATOM  -S-------{chain} {start},{end+1}, {calpha_color[0]:.1f},{calpha_color[1]:.1f},{calpha_color[2]:.1f}, {self.illustrate_config.carbon_radius}\n")
        file.write(f"ATOM  ---------{chain} {start},{end+1}, {sidechain_color[0]:.1f},{sidechain_color[1]:.1f},{sidechain_color[2]:.1f}, {self.illustrate_config.sidechain_radius}\n")

    def write_chain_colors(self, file):
        for chain in self.protein_config.chain_to_color:
            calpha_color = self.protein_config.chain_to_color[chain]
            sidechain_color = image_utils.alpha_blending(calpha_color, self.illustrate_config.sidechain_transparency)
            self.write_residue_range(file, chain, 0, 9999, calpha_color, sidechain_color)

    def write_end_and_transformations(self, file):
        file.write(f"""END
center
{self.illustrate_config.center}
trans
{','.join(map(lambda x: f"{x}", self.illustrate_config.translation))}
scale
{self.illustrate_config.scale}
xrot
{self.illustrate_config.rotation[0]}
yrot
{self.illustrate_config.rotation[1]}
zrot
{self.illustrate_config.rotation[2]}
wor
{','.join(f"{x}" for x in m_colors.to_rgb(self.illustrate_config.background_color))},{','.join(f"{x}" for x in m_colors.to_rgb(self.illustrate_config.fog_color))},{self.illustrate_config.fog_front_transparency},{self.illustrate_config.fog_back_transparency}
{int(self.illustrate_config.shadow)},{self.illustrate_config.shadow_cone_fraction},{self.illustrate_config.shadow_cone_angle},{self.illustrate_config.shadow_cone_difference},{self.illustrate_config.shadow_cone_max}
{self.illustrate_config.padding[0]},{self.illustrate_config.padding[1]}
illustrate
{self.illustrate_config.contour_outline_min},{self.illustrate_config.contour_outline_max},{self.illustrate_config.contour_ikernel},{self.illustrate_config.contour_difference_min},{self.illustrate_config.contour_difference_max}
{self.illustrate_config.subunit_outline_min},{self.illustrate_config.subunit_outline_max}
{self.illustrate_config.residue_outline_min},{self.illustrate_config.residue_outline_max},{self.illustrate_config.residue_difference}
calculate
{self.protein_config.output_prefix}_illustrate.pnm""")
        
            
    def run(self, remove_intermediate_files: bool = True):
        self.make_command_file()
        with open(f"{self.protein_config.output_prefix}_illustrate.inp", 'r') as file:
            subprocess.run([self.illustrate_config.illustrate_binary], stdin=file, check=True, stdout=subprocess.DEVNULL)
        subprocess.check_call([self.illustrate_config.convert_binary, f"{self.protein_config.output_prefix}_illustrate.pnm", f"{self.protein_config.output_prefix}_illustrate.png"])
        img = Image.open(f"{self.protein_config.output_prefix}_illustrate.png")
        img = img.rotate(90, expand=True)
        img.save(f"{self.protein_config.output_prefix}_illustrate.png")
        if remove_intermediate_files:
            Path(f"{self.protein_config.output_prefix}_illustrate.inp").unlink()
            Path(f"{self.protein_config.output_prefix}_illustrate.pnm").unlink()