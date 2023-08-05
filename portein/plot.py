import matplotlib.pyplot as plt
import numpy as np
import prody as pd
from matplotlib import patches as m_patches
from matplotlib import transforms as m_transforms
from matplotlib import colors as m_colors
from matplotlib import colormaps
import typing as ty
from pathlib import Path
from dataclasses import dataclass, field
from itertools import groupby
from operator import itemgetter
import subprocess
from PIL import Image
from pymol import cmd, util

from portein.config import HelixConfig, TurnConfig, SheetConfig, PorteinConfig
from portein.rotate import get_best_transformation, apply_transformation

SS_DICT = {
    "H": "H",
    "G": "H",
    "I": "H",
    "B": "E",
    "E": "E",
    "T": "T",
    "S": "T",
    "C": "T",
}


def plot_portrait(
        pdb: ty.Union[str, Path, pd.AtomGroup],
        config: PorteinConfig = None,
        height=12,
        width=None,
        ax=None,
):
    """
    Plot 2D portrait of a protein

    Parameters
    ----------
    pdb
        Can be a PDB ID, a PDB file, or a prody AtomGroup object
    config
        PorteinConfig object
    height
        figure height (auto-calculated from width if set to None)
    width
        figure width (auto-calculated from height if set to None)
    ax
        matplotlib ax to use, if None, makes new figure with specified height and width
    Returns
    -------
    matplotlib Figure, matplotlib Axes, 2D points corresponding to each residue
    """
    if config is None:
        config = PorteinConfig.default()
    if type(pdb) == str or type(pdb) == Path:
        structure = pd.parsePDB(pdb)
    else:
        structure = pdb
    dssp_file = pd.execDSSP(pdb)
    structure = pd.parseDSSP(dssp_file, structure)
    structure_alpha = structure.select("calpha")
    coords = structure_alpha.getCoords()
    matrix = get_best_transformation(coords)
    coords = apply_transformation(coords, matrix)[:, :2]
    ss_list = structure_alpha.getSecstrs()
    ss_elements = get_ss_elements(ss_list)
    if ax is None:
        fig, ax = plt.subplots(1, figsize=find_size(coords, height, width))
    ax.axis("off")
    min_x, min_y, max_x, max_y = np.inf, np.inf, -np.inf, -np.inf
    for (ss, start_i, end_i) in ss_elements:
        if ss == "H":
            if config.helix.as_cylinder:
                ss = "HC"
            else:
                ss = "HW"
        start_x, start_y = coords[start_i, 0], coords[start_i, 1]
        end_x, end_y = coords[end_i, 0], coords[end_i, 1]
        min_x, min_y, max_x, max_y = update_limits(
            start_x, start_y, end_x, end_y, min_x, min_y, max_x, max_y
        )
        for patch in make_patch(ss, start_x, start_y, end_x, end_y, ax, config):
            ax.add_patch(patch)
    ax.set_xlim(min_x - 1, max_x + 1)
    ax.set_ylim(min_y - 1, max_y + 1)
    for direction in ["left", "right", "top", "bottom"]:
        ax.spines[direction].set_visible(False)
    return ax, coords


def make_helix_wave(config: HelixConfig, length):
    patches = []
    start_theta, end_theta = 0, 180
    for arc_start in np.arange(0, length, config.wave_arc_length):
        origin = (arc_start + 0.25, 0)
        patches.append(
            m_patches.Arc(
                origin,
                config.wave_arc_length,
                config.wave_arc_height,
                linewidth=config.wave_arc_width,
                # Add a bit to each angle to avoid sharp cuts
                # that show as white lines in plot
                theta1=start_theta - 1,
                theta2=end_theta + 1,
                edgecolor=config.color,
            )
        )
        start_theta += 180
        end_theta += 180
    return patches


def make_helix_cylinder(config: HelixConfig, length):
    patches = []
    # Origin is *center* of ellipse
    origin = (config.cylinder_ellipse_length / 2, 0)
    # First ellipse
    patches.append(
        m_patches.Ellipse(
            origin,
            config.cylinder_ellipse_length,
            config.cylinder_ellipse_height,
            linewidth=config.outline_width,
            edgecolor=config.outline_color,
            facecolor=config.color,
            alpha=config.opacity,
        )
    )

    # Rectangle(s)
    origin = (
        config.cylinder_ellipse_length / 2,
        -config.cylinder_ellipse_height / 2,
    )  # origin is lower left: make it v-cntr
    patches.append(
        m_patches.Rectangle(
            origin,
            length - config.cylinder_ellipse_length,
            config.cylinder_rectangle_height,
            edgecolor=config.outline_color,
            facecolor=config.color,
            linewidth=config.outline_width,
            alpha=config.opacity,
        )
    )

    # Second ellipse
    origin = (length - config.cylinder_ellipse_length / 2, 0)
    patches.append(
        m_patches.Ellipse(
            origin,
            config.cylinder_ellipse_length,
            config.cylinder_ellipse_height,
            linewidth=config.outline_width,
            edgecolor=config.outline_color,
            facecolor=config.color,
            alpha=config.opacity,
        )
    )
    return patches


def make_sheet(config: SheetConfig, length):
    return [
        m_patches.FancyArrow(
            0,
            0,  # x, y of tail
            length,
            0,  # dx, dy=0 -> flat arrow
            length_includes_head=True,
            head_length=length / 4,
            head_width=config.head_height - 0.001,
            width=config.tail_height,
            facecolor=config.color,
            edgecolor=config.outline_color,
            linewidth=config.outline_width,
            alpha=config.opacity,
        )
    ]


def make_turn(config: TurnConfig, length):
    origin = (length / 2, config.height / 2)
    return [
        m_patches.Circle(
            (0, 0),
            config.circle_radius,
            color=config.circle_color,
            alpha=config.opacity,
        ),
        m_patches.Arc(
            origin,
            length,
            config.height,
            linewidth=config.arc_width,
            # Add a bit to each angle to avoid sharp cuts
            # that show as white lines in plot
            theta1=0,
            theta2=180,
            edgecolor=config.arc_color,
            alpha=config.opacity,
        ),
        m_patches.Circle(
            (length, 0),
            config.circle_radius,
            color=config.circle_color,
            alpha=config.opacity,
        ),
    ]


def rotate_patch(patch, start_x, start_y, end_x, end_y, ax):
    rotation_angle = np.arctan2(end_y - start_y, end_x - start_x)
    rotation = m_transforms.Affine2D().rotate_deg(np.rad2deg(rotation_angle))
    translation = m_transforms.Affine2D().translate(start_x, start_y)
    transform = rotation + translation + ax.transData
    patch.set_transform(transform)
    return patch


def make_patch(name, sx, sy, ex, ey, ax, config: PorteinConfig):
    length = np.sqrt((ex - sx) ** 2 + (ey - sy) ** 2)
    if name == "HC":
        patches = make_helix_cylinder(config.helix, length)
    elif name == "HW":
        patches = make_helix_wave(config.helix, length)
    elif name == "E":
        patches = make_sheet(config.sheet, length)
    elif name == "T":
        patches = make_turn(config.turn, length)
    else:
        raise ValueError("name must be one of HC, HW, E, T")
    return [rotate_patch(patch, sx, sy, ex, ey, ax) for patch in patches]


def get_ss_elements(ss_list):
    ss_blocks = []
    prev_ss = None
    prev_i = 0
    for i, ss in enumerate(ss_list):
        ss = SS_DICT.get(ss, "T")
        if prev_ss is None:
            prev_ss = ss
        if ss != prev_ss:
            ss_blocks.append((prev_ss, prev_i, i))
            prev_ss = ss
            prev_i = i
    ss_blocks.append((SS_DICT.get(ss_list[-1], "T"), prev_i, len(ss_list) - 1))
    return ss_blocks


def update_limits(sx, sy, ex, ey, min_x, min_y, max_x, max_y):
    min_x = min(sx, min_x)
    min_x = min(ex, min_x)
    max_x = max(sx, max_x)
    max_x = max(ex, max_x)
    min_y = min(sy, min_y)
    min_y = min(ey, min_y)
    max_y = max(sy, max_y)
    max_y = max(ey, max_y)
    return min_x, min_y, max_x, max_y


def find_size(points, height: ty.Optional[float], width: ty.Optional[float]):
    """
    Find figure size given coordinates and size of one dimension.
    Calculates height from width if height=None and width from height if width=None

    Parameters
    ----------
    points
    height
    width

    Returns
    -------
    (width, height)
    """
    assert (
            width is not None or height is not None
    ), "Either one of height or width must be set"
    if height is not None and width is not None:
        return width, height
    min_x, max_x = np.min(points[:, 0]), np.max(points[:, 0])
    min_y, max_y = np.min(points[:, 1]), np.max(points[:, 1])
    if height is None:
        return width, (max_y - min_y) * width / (max_x - min_x)
    else:
        return (max_x - min_x) * height / (max_y - min_y), height

def put_alpha(img, transparency):
    im2 = img.copy()
    im2.putalpha(int(255 * (1 - transparency)))
    img.paste(im2, img)
    return img


def alpha_blending(color, alpha):
    return tuple( (1. -  alpha) + np.array(color)*alpha )

def get_chain_colors(chain_colormap, pdb_file):
    if isinstance(chain_colormap, str):
        chain_to_color = {}
        chains = sorted(set(pd.parsePDB(str(pdb_file)).getChids()))
        colormap = colormaps.get(chain_colormap, None)
        for i, chain in enumerate(chains):
            if colormap is None:
                chain_to_color[chain] = m_colors.to_rgb(chain_colormap)
            else:
                chain_to_color[chain] = colormap(i)
    else:
        chain_to_color = {c: m_colors.to_rgb(chain_colormap[c]) for c in chain_colormap}
    return chain_to_color
    
def get_residues_colors(highlight_residues):
    chain_to_residue_range_color = {}
    for chain in highlight_residues:
        chain_to_residue_range_color[chain] = {}
        for color, residues in highlight_residues[chain].items():
            residues = sorted(residues)
            residue_ranges = []
            for _, g in groupby(enumerate(residues), lambda ix: ix[0] - ix[1]):
                group = list(map(itemgetter(1), g))
                residue_ranges.append((group[0], group[-1]))
            chain_to_residue_range_color[chain][color] = residue_ranges
    return chain_to_residue_range_color

@dataclass
class Illustrate:
    pdb_file: str
    """PDB file to read"""
    illustrate_binary: str = "illustrate"
    """path to illustrate binary"""
    convert_binary: str = "convert"
    """path to convert binary"""
    output_prefix: str = None
    """prefix for output files. If None uses the PDB file name"""
    chain_colormap: ty.Union[str, dict] = "Set3"
    """colormap to use for coloring chains, either a matplotlib colormap or a dictionary of {chain: color}"""
    highlight_residues: dict = field(default_factory=dict)
    """dictionary of {chain: {color: [residues]}}"""
    center: str = 'auto'
    """center of the image, one of 'auto' or 'center'"""
    translation: ty.Tuple[float, float, float] = (0., 0., 0.)
    """translation x, y, z in Angstroms"""
    scale: float = 10.0
    """scale (pixels/Angstrom), controls size of image"""
    rotation: ty.Tuple[float, float, float] = (0., 0., 0.)
    """rotation x, y, z in degrees"""
    background_color = "white"
    """background color"""
    fog_color = "white"
    """fog color"""
    fog_front_transparency: float = 1.0
    """fractional transparency of fog at front of molecule (0.0-1.0, 1.0=no fog)"""
    fog_back_transparency: float = 1.0
    """fractional transparency of fog at back of molecule (0.0-1.0, 1.0=no fog)"""
    sidechain_transparency: float = 0.8
    """fractional transparency of sidechain atoms (0.0-1.0, 1.0=same color as backbone)"""
    shadow: bool = True
    """whether to include shadows"""
    shadow_cone_fraction: float = 0.0023
    """fractional shadowing around each atom. larger=darker (0.0-1.0, typically 0.0023)"""
    shadow_cone_angle: float = 2.0
    """angle of shadowing around each atom. larger=tighter region (typically 2.0)"""
    shadow_cone_difference: float = 1.0
    """shadowing only applied if z-difference greater than this value (Angstroms, typically 1.0)"""
    shadow_cone_max: float = 0.7
    """maximal shadowing amount. smaller=darker (0.0-1.0, typically 0.7)"""
    padding: tuple = (-2, -2)
    """padding around molecule (width, height)"""
    contour_outline_min: float = 3.0
    contour_outline_max: float = 10.0
    """thresholds for gray to black. typically values from about 3.0-20.0, 
    best values for typical atomic illustrations: 3.0, 10.0"""
    contour_ikernel: int = 4
    """kernel for derivative calculation (1,2,3,4 smoothest=4)"""
    contour_difference_min: float = 0.0
    contour_difference_max: float = 5.0
    """range of z-difference used for derivative (Angstroms)
        0.0,1.0 gives outlines around every atom
        0.0,1000.0 gives only outline around molecule
        0.0,5.0 is typical
    """
    subunit_outline_min: float = 3.0
    subunit_outline_max: float = 10.0
    """thresholds for gray to black (typically ~ 3.0-20.0)"""
    residue_outline_min: float = 3.0
    residue_outline_max: float = 10.0
    """thresholds for gray to black (typically ~ 3.0-20.0)"""
    residue_difference: float = 1.0
    """difference in residue numbers to draw outlines"""
    carbon_radius: float = 1.6
    """radius of carbon atoms (Angstroms)"""
    sidechain_radius: float = 1.5
    """radius of sidechain atoms (Angstroms)"""
    width: int = None
    height: int = None

    @property
    def _output_prefix(self):
        return self.output_prefix if self.output_prefix is not None else Path(self.pdb_file).stem

    def make_command_file(self):
        chain_to_color = get_chain_colors(self.chain_colormap, self.pdb_file)
        chain_to_residue_range_color = get_residues_colors(self.highlight_residues)
        with open(f"{self._output_prefix}.inp", "w") as f:
            # omit hydrogens and water
            f.write(f"""read
{self.pdb_file}
HETATM-----HOH-- 0,9999, 0.5,0.5,0.5, 0.0
ATOM  -H-------- 0,9999, 0.5,0.5,0.5, 0.0
ATOM  H--------- 0,9999, 0.5,0.5,0.5, 0.0
""")        
            
            for chain in chain_to_color:
                calpha_color = chain_to_color[chain]
                sidechain_color = alpha_blending(calpha_color, self.sidechain_transparency)
                if chain in self.highlight_residues:
                    for color, residue_ranges in chain_to_residue_range_color[chain].items():
                        color = m_colors.to_rgb(color)
                        for start, end in residue_ranges:
                            f.write(f"ATOM  -C-------{chain} {start},{end+1}, {color[0]:.1f},{color[1]:.1f},{color[2]:.1f}, {self.carbon_radius}\n")
                            f.write(f"ATOM  -S-------{chain} {start},{end+1}, {color[0]:.1f},{color[1]:.1f},{color[2]:.1f}, {self.carbon_radius}\n")
                            f.write(f"ATOM  ---------{chain} {start},{end+1}, {color[0]:.1f},{color[1]:.1f},{color[2]:.1f}, {self.sidechain_radius}\n")
                
            for chain in chain_to_color:
                calpha_color = chain_to_color[chain]
                sidechain_color = alpha_blending(calpha_color, self.sidechain_transparency)
                f.write(f"ATOM  -C-------{chain} 0,9999, {calpha_color[0]:.1f},{calpha_color[1]:.1f},{calpha_color[2]:.1f}, {self.carbon_radius}\n")
                f.write(f"ATOM  -S-------{chain} 0,9999, {calpha_color[0]:.1f},{calpha_color[1]:.1f},{calpha_color[2]:.1f}, {self.carbon_radius}\n")
                f.write(f"ATOM  ---------{chain} 0,9999, {sidechain_color[0]:.1f},{sidechain_color[1]:.1f},{sidechain_color[2]:.1f}, {self.sidechain_radius}\n")
            
            f.write(f"""END
center
{self.center}
trans
{','.join(map(lambda x: f"{x}", self.translation))}
scale
{self.scale}
xrot
{self.rotation[0]}
yrot
{self.rotation[1]}
zrot
{self.rotation[2]}
wor
{','.join(f"{x}" for x in m_colors.to_rgb(self.background_color))},{','.join(f"{x}" for x in m_colors.to_rgb(self.fog_color))},{self.fog_front_transparency},{self.fog_back_transparency}
{int(self.shadow)},{self.shadow_cone_fraction},{self.shadow_cone_angle},{self.shadow_cone_difference},{self.shadow_cone_max}
{self.padding[0]},{self.padding[1]}
illustrate
{self.contour_outline_min},{self.contour_outline_max},{self.contour_ikernel},{self.contour_difference_min},{self.contour_difference_max}
{self.subunit_outline_min},{self.subunit_outline_max}
{self.residue_outline_min},{self.residue_outline_max},{self.residue_difference}
calculate
{self._output_prefix}.pnm""")
            
    def run(self, remove_intermediate_files: bool = True):
        self.make_command_file()
        with open(f"{self._output_prefix}.inp", 'r') as file:
            subprocess.run([self.illustrate_binary], stdin=file, check=True, stdout=subprocess.DEVNULL)
        subprocess.check_call([self.convert_binary, f"{self._output_prefix}.pnm", f"{self._output_prefix}.png"])
        img = Image.open(f"{self._output_prefix}.png")
        img = img.rotate(90, expand=True)
        self.width = self.width if self.width is not None else img.width
        self.height = self.height if self.height is not None else img.height
        img.save(f"{self._output_prefix}.png")
        if remove_intermediate_files:
            Path(f"{self._output_prefix}.inp").unlink()
            Path(f"{self._output_prefix}.pnm").unlink()

@dataclass
class PymolConfig:
    ray_trace_mode: int = 1
    surface_quality: int = 2
    cartoon_sampling: int = 20
    ambient: float = 0.5
    cartoon_discrete_colors: bool = True
    ray_opaque_background: bool = False
    cartoon_fancy_helices: bool = True
    antialias: int = 2
    ray_trace_gain: float = 0
    ray_trace_disco_factor: int = 1
    ray_texture: int = 0
    ray_trace_fog: bool = False
    hash_max: int = 300
    depth_cue: bool = False
    ray_shadows: bool = False
    light_count: int = 1
    specular: bool = False

@dataclass
class RepresentationConfig:
    representation: str
    pymol_config: PymolConfig
    selection: str
    transparency: float
    color: str = None

@dataclass
class Pymol:
    pdb_file: str
    """PDB file to read"""
    pymol_binary: str = "pymol"
    """path to pymol binary"""
    output_prefix: str = None
    """prefix for output files. If None uses the PDB file name"""
    chain_colormap: ty.Union[str, dict] = "Set3"
    """colormap to use for coloring chains, either a matplotlib colormap or a dictionary of {chain: color}"""
    highlight_residues: dict = field(default_factory=dict)
    """dictionary of {chain: {color: [residue numbers]}}, use None to set the color to the chain color"""
    figure_order: ty.List[RepresentationConfig] = None
    """order of representations, selections, and transparencies to layer on top of each other"""
    orient: bool = False
    """re-orient the protein using pymol"""
    width: int = None
    """width of image"""
    height: int = None
    """height of image"""
    buffer: float = None
    
    def __post_init__(self):
        if self.figure_order is None:
            self.figure_order = [RepresentationConfig(representation="surface",
                                                pymol_config=PymolConfig(), 
                                                selection="all",
                                                transparency=0.5), 
                                RepresentationConfig(representation="cartoon",
                                                    pymol_config=PymolConfig(), 
                                                    selection="all",
                                                    transparency=0.),
                                RepresentationConfig(representation="stick",
                                                pymol_config=PymolConfig(),
                                                selection="hetatm",
                                                transparency=0., 
                                                color="skyblue"),]
            if len(self.highlight_residues) > 0:
                self.figure_order.append(RepresentationConfig(representation="stick",
                                                    pymol_config=PymolConfig(), 
                                                    selection="highlight",
                                                    transparency=0.))
        pdb = pd.parsePDB(self.pdb_file)
        if self.width is None or self.height is None:
            self.width, self.height = find_size(pdb.getCoords(), self.width, self.height)

    @property
    def _output_prefix(self):
        return self.output_prefix if self.output_prefix is not None else Path(self.pdb_file).stem
    
    def draw_protein(self):
        cmd.reinitialize()
        cmd.load(self.pdb_file)
        cmd.bg_color("white")
        cmd.remove("solvent")
        cmd.set("max_threads", 8)
        if self.orient:
            cmd.orient()
        if self.buffer is not None:
            cmd.zoom(buffer=self.buffer)
        chain_to_color = get_chain_colors(self.chain_colormap, self.pdb_file)
        for chain, color in chain_to_color.items():
            cmd.color(f"0x{m_colors.to_hex(color).lower()[1:]}", f"chain {chain}")
    
        for representation_config in self.figure_order:
            config = representation_config.pymol_config.__dict__
            for key, value in config.items():
                if key in ["specular", "depth_cue"]:
                    if value:
                        cmd.set(key)
                    else:
                        cmd.unset(key)
                else:
                    cmd.set(key, value)
            cmd.hide("everything")
            if representation_config.selection == "highlight":
                for chain, colors in self.highlight_residues.items():
                    for color, residues in colors.items():
                        if color is None:
                            color = chain_to_color[chain]
                        if representation_config.color is not None:
                            color = representation_config.color
                        color = m_colors.to_hex(color).lower()[1:]
                        cmd.select(f"highlight_{chain}_{color}", f"chain {chain} and resi {'+'.join(str(x) for x in residues)}")
                        cmd.show(representation_config.representation, f"highlight_{chain}_{color}")
                        cmd.color(f"0x{color}", f"highlight_{chain}_{color}")    
            else:
                cmd.show(representation_config.representation, representation_config.selection)
            if representation_config.representation == "stick":
                util.cnc(f"rep stick")
            if representation_config.color is not None:
                cmd.color(representation_config.color, representation_config.selection)
            cmd.ray(self.width, self.height)
            cmd.png(f"{self._output_prefix}_{representation_config.representation}.png", width=self.width, height=self.height, dpi=300)
        self.layer()
        cmd.reinitialize()

    def layer(self):
        pngs = []
        for representation_config in self.figure_order:
            img = Image.open(f"{self._output_prefix}_{representation_config.representation}.png")
            pngs.append(put_alpha(img, representation_config.transparency))
        image = Image.new("RGBA", pngs[0].size)
        for png in pngs:
            image = Image.alpha_composite(image, png)
        image.save(f"{self._output_prefix}.png")