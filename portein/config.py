from dataclasses import dataclass
from pathlib import Path
import typing

from matplotlib.typing import ColorType
from portein.plot import image_utils
from matplotlib import colormaps
from itertools import groupby
from operator import itemgetter
from matplotlib import colors as m_colors
from portein.rotate import rotate_protein
import yaml
import biotite.structure.io as bio_io
from biotite.structure.io import pdb, pdbx
import biotite.structure as bio_struct
from biotite.database import rcsb
from biotite.structure import AtomArray
from tempfile import gettempdir


def read_structure(path: typing.Union[str, Path]) -> AtomArray:
    if Path(path).suffix == ".pdb":
        return pdb.PDBFile.read(path).get_structure(model=1)
    elif Path(path).suffix == ".cif":
        return pdbx.get_structure(path, model=1)
    else:
        try:
            cif_file_path = rcsb.fetch(path, "cif", gettempdir())
            cif_file = pdbx.CIFFile.read(cif_file_path)
            return pdbx.get_structure(cif_file, model=1)
        except Exception as e:
            raise ValueError(f"Unsupported file: {path}") from e


@dataclass
class HelixConfig:
    """
    Config related to drawing helices

    Parameters
    ----------
    as_cylinder
        If True draws helix as a rectangle capped by two ellipses,
        otherwise as consecutive elliptical arcs
    cylinder_ellipse_length
    cylinder_ellipse_height
    cylinder_rectangle_height
    wave_arc_width
    wave_arc_height
    wave_arc_length
    outline_width
    outline_color
    color
    opacity
    """

    as_cylinder: bool = False
    cylinder_ellipse_length: float = 0.5
    cylinder_ellipse_height: float = 1.0 - 0.001
    cylinder_rectangle_height: float = 1.0
    wave_arc_width: float = 3.0
    wave_arc_height: float = 1.0
    wave_arc_length: float = 0.5
    outline_width: float = 3.0
    outline_color: ColorType = "#5c6887"
    color: ColorType = "lightsteelblue"
    opacity: float = 1.0

    def change_width(self, width=1.0):
        self.cylinder_ellipse_length = width / 2
        self.cylinder_ellipse_height = width - 0.001
        self.cylinder_rectangle_height = width
        self.wave_arc_height = width

    @classmethod
    def from_yaml(cls, yaml_str: typing.Union[str, Path]):
        """
        Create from a yaml string or file
        """
        if yaml_str is None:
            return cls()
        loader = yaml.SafeLoader
        loader.add_constructor("!ColorType", ColorType)
        if Path(yaml_str).is_file():
            yaml_str = Path(yaml_str).read_text()
        return cls(**yaml.load(yaml_str, Loader=loader))


@dataclass
class SheetConfig:
    """
    Config related to drawing beta sheets

    Parameters
    ----------
    thickness_factor
    tail_height
    head_height
    outline_width
    outline_color
    color
    opacity
    """

    thickness_factor: float = 1
    tail_height: float = 1.0
    head_height: float = 2.0
    outline_width: float = 3.0
    outline_color: ColorType = "#5c6887"
    color: ColorType = "#999FD0"
    opacity: float = 1.0

    @classmethod
    def from_yaml(cls, yaml_str: typing.Union[str, Path]):
        """
        Create from a yaml string or file
        """
        if yaml_str is None:
            return cls()
        loader = yaml.SafeLoader
        loader.add_constructor("!ColorType", ColorType)
        if Path(yaml_str).is_file():
            yaml_str = Path(yaml_str).read_text()
        return cls(**yaml.load(yaml_str, Loader=loader))


@dataclass
class TurnConfig:
    """
    Config related to drawing turns

    Parameters
    ----------
    thickness_factor
    height
    circle_radius
    circle_color
    arc_width
    arc_color
    opacity
    """

    thickness_factor: float = 0.5
    height: float = 0.5
    circle_radius: float = 0.2
    circle_color: ColorType = "#d1d6e3"
    arc_width: float = 3.0
    arc_color: ColorType = "#d1d6e3"
    opacity: float = 0.8

    @classmethod
    def from_yaml(cls, yaml_str: typing.Union[str, Path]):
        """
        Create from a yaml string or file
        """
        if yaml_str is None:
            return cls()
        loader = yaml.SafeLoader
        loader.add_constructor("!ColorType", ColorType)
        if Path(yaml_str).is_file():
            yaml_str = Path(yaml_str).read_text()
        return cls(**yaml.load(yaml_str, Loader=loader))


@dataclass
class ProteinConfig:
    pdb_file: str
    """PDB file to read"""
    rotate: bool = False
    """rotate the protein to have the best orientation"""
    output_prefix: str = None
    """prefix for output files. If None uses the PDB file name"""
    chain_colormap: typing.Union[str, typing.Dict[str, ColorType]] = "Set3"
    """colormap to use for coloring chains, either a matplotlib colormap or a dictionary of {chain: color}"""
    highlight_residues: typing.Dict[str, typing.Dict[ColorType, typing.List[int]]] = (
        None
    )
    """dictionary of {chain: {color: [residue numbers]}}, use None to set the color to the chain color"""
    width: int = 1000
    """width of image (in pixels)"""
    height: int = None
    """height of image"""
    chain_to_color: typing.Dict[str, ColorType] = None
    """dictionary of {chain: color}"""
    chain_to_residue_range_color: typing.Dict[
        str, typing.Dict[typing.Tuple[int, int], ColorType]
    ] = None
    """dictionary of {chain: {residue_range: color}}"""

    def __post_init__(self):
        self.finish()

    def finish(self):
        if "." in str(self.pdb_file):
            self.pdb_file = Path(self.pdb_file).resolve()
        if self.output_prefix is None:
            self.output_prefix = Path(self.pdb_file).stem
        if self.rotate:
            self.save_rotated()
        self.width, self.height = image_utils.find_size(
            self.pdb.coord, width=self.width, height=self.height
        )
        if self.highlight_residues is None:
            self.highlight_residues = {}
        if self.chain_to_color is None:
            self.chain_to_color = self.get_chain_colors()
        if self.chain_to_residue_range_color is None:
            self.chain_to_residue_range_color = self.get_residues_colors()

    @classmethod
    def from_yaml(cls, yaml_str: typing.Union[str, Path]):
        """
        Create from a yaml string or file
        """
        loader = yaml.SafeLoader
        loader.add_constructor("!ColorType", ColorType)
        if Path(yaml_str).is_file():
            yaml_str = Path(yaml_str).read_text()
        return cls(**yaml.load(yaml_str, Loader=loader))

    @property
    def pdb(self) -> AtomArray:
        return read_structure(self.pdb_file)

    def save_rotated(self):
        pdb = rotate_protein(self.pdb)
        self.output_prefix = f"{self.output_prefix}_rotated"
        bio_io.save_structure(f"{self.output_prefix}.pdb", pdb)
        self.pdb_file = f"{self.output_prefix}.pdb"

    def get_chain_colors(self):
        if isinstance(self.chain_colormap, str):
            chain_to_color = {}
            chains = sorted(set(bio_struct.get_chains(self.pdb)))
            colormap = colormaps.get(self.chain_colormap, None)
            for i, chain in enumerate(chains):
                if colormap is None:
                    chain_to_color[chain] = m_colors.to_rgb(self.chain_colormap)
                else:
                    chain_to_color[chain] = colormap(i)
        else:
            chain_to_color = {
                c: m_colors.to_rgb(self.chain_colormap[c]) for c in self.chain_colormap
            }
        return chain_to_color

    def get_residues_colors(self):
        chain_to_residue_range_color = {}
        for chain in self.highlight_residues:
            chain_to_residue_range_color[chain] = {}
            for color, residues in self.highlight_residues[chain].items():
                residues = sorted(residues)
                residue_ranges = []
                for _, g in groupby(enumerate(residues), lambda ix: ix[0] - ix[1]):
                    group = list(map(itemgetter(1), g))
                    residue_ranges.append((group[0], group[-1]))
                chain_to_residue_range_color[chain][color] = residue_ranges
        return chain_to_residue_range_color


PYMOL_SETTINGS = dict(
    ambient=0.5,
    antialias=2,
    cartoon_discrete_colors=True,
    cartoon_fancy_helices=True,
    cartoon_sampling=20,
    depth_cue=False,
    hash_max=300,
    light_count=1,
    ray_opaque_background=False,
    ray_shadows=False,
    ray_texture=0,
    ray_trace_disco_factor=1,
    ray_trace_fog=False,
    ray_trace_gain=0,
    ray_trace_mode=1,
    specular=False,
    surface_quality=2,
)


@dataclass
class PymolConfig:
    representation: str = "cartoon"
    pymol_settings: dict = None
    selection: str = "all"
    transparency: float = 0.0
    color: ColorType = None
    spectrum: str = None
    dpi: int = 300

    def __post_init__(self):
        if self.representation not in [
            "cartoon",
            "surface",
            "sticks",
            "spheres",
            "lines",
            "mesh",
            "dots",
            "ribbon",
            "nonbonded",
        ]:
            raise ValueError(f"{self.representation} is not a valid representation")
        if self.pymol_settings is None:
            self.pymol_settings = PYMOL_SETTINGS.copy()

    @classmethod
    def from_yaml(cls, yaml_str: typing.Union[str, Path]):
        """
        Create from a yaml string or file
        """
        loader = yaml.SafeLoader
        loader.add_constructor("!ColorType", ColorType)
        if Path(yaml_str).is_file():
            yaml_str = Path(yaml_str).read_text()
        return cls(**yaml.load(yaml_str, Loader=loader))

    @classmethod
    def _from_yaml(cls, loader, node):
        data = loader.construct_mapping(node)
        return cls(**data)


@dataclass
class IllustrateConfig:
    illustrate_binary: str = "illustrate"
    """path to illustrate binary"""
    convert_binary: str = "magick"
    """path to convert binary"""
    center: str = "auto"
    """center of the image, one of 'auto' or 'center'"""
    translation: typing.Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """translation x, y, z in Angstroms"""
    scale: float = 10.0
    """scale (pixels/Angstrom), controls size of image"""
    rotation: typing.Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """rotation x, y, z in degrees"""
    background_color: ColorType = "white"
    """background color"""
    fog_color: ColorType = "white"
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

    @classmethod
    def from_yaml(cls, yaml_str: typing.Union[str, Path]):
        """
        Create from a yaml string or file
        """
        if yaml_str is None:
            return cls()
        loader = yaml.SafeLoader
        loader.add_constructor("!ColorType", ColorType)
        if Path(yaml_str).is_file():
            yaml_str = Path(yaml_str).read_text()
        return cls(**yaml.load(yaml_str, Loader=loader))
