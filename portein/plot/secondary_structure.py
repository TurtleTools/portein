from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches as m_patches
from matplotlib import transforms as m_transforms
import typing as ty
from dataclasses import dataclass
from portein import config
from portein.plot import image_utils
from portein.config import read_structure
from subprocess import SubprocessError
from tempfile import NamedTemporaryFile
from biotite.application.application import AppState, requires_state
from biotite.application.localapp import LocalApp, cleanup_tempfile, get_version
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdb.convert import set_structure

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


class DsspApp(LocalApp):
    r"""
    [Modified from biotite.application.dssp.DsspApp]

    Annotate the secondary structure of a protein structure using the
    *DSSP* software.

    Internally this creates a :class:`Popen` instance, which handles
    the execution.

    DSSP differentiates between 8 different types of secondary
    structure elements:

       - C: loop, coil or irregular
       - H: :math:`{\alpha}`-helix
       - B: :math:`{\beta}`-bridge
       - E: extended strand, participation in :math:`{\beta}`-ladder
       - G: 3 :sub:`10`-helix
       - I: :math:`{\pi}`-helix
       - T: hydrogen bonded turn
       - S: bend

    Parameters
    ----------
    atom_array : AtomArray
        The atom array to be annotated.
    bin_path : str, optional
        Path of the *DDSP* binary.

    Examples
    --------

    >>> app = DsspApp(atom_array)
    >>> app.start()
    >>> app.join()
    >>> print(app.get_sse())
    ['C' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'T' 'T' 'G' 'G' 'G' 'G' 'T' 'C' 'C' 'C'
     'C' 'C']
    """

    def __init__(self, atom_array, bin_path="mkdssp"):
        super().__init__(bin_path)

        # mkdssp requires also the
        # 'occupancy', 'b_factor' and 'charge' fields
        # -> Add these annotations to a copy of the input structure
        self._array = atom_array.copy()
        categories = self._array.get_annotation_categories()
        if "charge" not in categories:
            self._array.set_annotation(
                "charge", np.zeros(self._array.array_length(), dtype=int)
            )
        if "b_factor" not in categories:
            self._array.set_annotation(
                "b_factor", np.zeros(self._array.array_length(), dtype=float)
            )
        if "occupancy" not in categories:
            self._array.set_annotation(
                "occupancy", np.ones(self._array.array_length(), dtype=float)
            )
        try:
            # The parameters have changed in version 4
            self._new_cli = get_version(bin_path)[0] >= 4
        except SubprocessError:
            # In older versions, the no version is returned with `--version`
            # -> a SubprocessError is raised
            self._new_cli = False
        self._in_file = NamedTemporaryFile("w", suffix=".pdb", delete=False)
        self._out_file = NamedTemporaryFile("r", suffix=".dssp", delete=False)

    def run(self):
        in_file = PDBFile()
        set_structure(in_file, self._array)
        in_file.write(self._in_file)
        self._in_file.flush()
        if self._new_cli:
            self.set_arguments([self._in_file.name, self._out_file.name, "-v"])
        else:
            self.set_arguments(["-i", self._in_file.name, "-o", self._out_file.name])
        super().run()

    def evaluate(self):
        super().evaluate()
        lines = self._out_file.read().split("\n")
        # Index where SSE records start
        sse_start = None
        for i, line in enumerate(lines):
            if line.startswith("  #  RESIDUE AA STRUCTURE"):
                sse_start = i + 1
        if sse_start is None:
            raise ValueError("DSSP file does not contain SSE records")
        # Remove "!" for missing residues
        lines = [
            line for line in lines[sse_start:] if len(line) != 0 and line[13] != "!"
        ]
        self._sse = np.zeros(len(lines), dtype="U1")
        # Parse file for SSE letters
        for i, line in enumerate(lines):
            self._sse[i] = line[16]
        self._sse[self._sse == " "] = "C"

    def clean_up(self):
        super().clean_up()
        cleanup_tempfile(self._in_file)
        cleanup_tempfile(self._out_file)

    @requires_state(AppState.JOINED)
    def get_sse(self):
        """
        Get the resulting secondary structure assignment.

        Returns
        -------
        sse : ndarray, dtype="U1"
            An array containing DSSP secondary structure symbols
            corresponding to the residues in the input atom array.
        """
        return self._sse

    @staticmethod
    def annotate_sse(atom_array, bin_path="mkdssp"):
        """
        Perform a secondary structure assignment to an atom array.

        This is a convenience function, that wraps the :class:`DsspApp`
        execution.

        Parameters
        ----------
        atom_array : AtomArray
            The atom array to be annotated.
        bin_path : str, optional
            Path of the DDSP binary.

        Returns
        -------
        sse : ndarray, dtype="U1"
            An array containing DSSP secondary structure symbols
            corresponding to the residues in the input atom array.
        """
        app = DsspApp(atom_array, bin_path)
        app.start()
        app.join()
        return app.get_sse()


@dataclass
class SecondaryStructure:
    protein_config: config.ProteinConfig
    helix_config: config.HelixConfig = None
    sheet_config: config.SheetConfig = None
    turn_config: config.TurnConfig = None
    dpi: int = 300
    coords: ty.Optional[np.ndarray] = None
    ss_elements: ty.Optional[ty.List[ty.Tuple[str, int, int]]] = None

    def __post_init__(self):
        if self.helix_config is None:
            self.helix_config = config.HelixConfig()
        if self.sheet_config is None:
            self.sheet_config = config.SheetConfig()
        if self.turn_config is None:
            self.turn_config = config.TurnConfig()

    def from_yaml(
        protein_config: Path,
        helix_config: Path,
        sheet_config: Path,
        turn_config: Path,
        dpi: int = 300,
    ):
        protein = config.ProteinConfig.from_yaml(protein_config)
        helix = config.HelixConfig.from_yaml(helix_config)
        sheet = config.SheetConfig.from_yaml(sheet_config)
        turn = config.TurnConfig.from_yaml(turn_config)
        return SecondaryStructure(protein, helix, sheet, turn, dpi)

    def run(self, ax=None, linear=False, y_offset=0, overwrite=False):
        """
        Plot 2D portrait of a protein

        Parameters
        ----------
        ax
            matplotlib ax to use, if None, makes new figure with specified height and width
        linear
            If True, draws all secondary structure elements in a single line instead of using protein coordinates
        y_offset
            If linear is True, y coordinate of the line to draw secondary structure elements on
        Returns
        -------
        matplotlib Axes
        """
        structure = read_structure(self.protein_config.pdb_file)
        sse = DsspApp.annotate_sse(structure[~structure.hetero])
        self.coords = structure[(structure.atom_name == "CA") & ~structure.hetero].coord
        assert len(self.coords) == len(sse), (
            f"Number of coordinates and secondary structure elements must be the same, got {len(self.coords)} and {len(sse)}"
        )
        self.ss_elements = get_ss_elements(sse)
        if ax is None:
            fig, ax = plt.subplots(
                1,
                figsize=(
                    self.protein_config.width / self.dpi,
                    self.protein_config.height / self.dpi,
                ),
                dpi=self.dpi,
            )
        ax.axis("off")
        return self.make_patches(ax, linear, y_offset)

    def make_patches(self, ax: plt.Axes, linear: bool = False, y_offset: float = 0):
        min_x, min_y, max_x, max_y = np.inf, np.inf, -np.inf, -np.inf
        for ss, start_i, end_i in self.ss_elements:
            if ss == "H":
                ss = "HC" if self.helix_config.as_cylinder else "HW"
            if linear:
                start_x, start_y = start_i, y_offset
                end_x, end_y = end_i, y_offset
            else:
                start_x, start_y = self.coords[start_i, 0], self.coords[start_i, 1]
                end_x, end_y = self.coords[end_i, 0], self.coords[end_i, 1]
            min_x, min_y, max_x, max_y = image_utils.update_limits(
                start_x, start_y, end_x, end_y, min_x, min_y, max_x, max_y
            )
            for patch in make_patch(
                ss,
                start_x,
                start_y,
                end_x,
                end_y,
                ax,
                self.helix_config,
                self.sheet_config,
                self.turn_config,
            ):
                ax.add_patch(patch)
        ax.set_xlim(min_x - 1, max_x + 1)
        ax.set_ylim(min_y - 1, max_y + 1)
        for direction in ["left", "right", "top", "bottom"]:
            ax.spines[direction].set_visible(False)
        return ax


def make_helix_wave(config: config.HelixConfig, length):
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


def make_helix_cylinder(config: config.HelixConfig, length):
    # Origin is *center* of ellipse
    origin = (config.cylinder_ellipse_length / 2, 0)
    # First ellipse
    patches = [
        m_patches.Ellipse(
            origin,
            config.cylinder_ellipse_length,
            config.cylinder_ellipse_height,
            linewidth=config.outline_width,
            edgecolor=config.outline_color,
            facecolor=config.color,
            alpha=config.opacity,
        )
    ]

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


def make_sheet(config: config.SheetConfig, length):
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


def make_turn(config: config.TurnConfig, length):
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


def make_patch(name, sx, sy, ex, ey, ax, helix_config, sheet_config, turn_config):
    length = np.sqrt((ex - sx) ** 2 + (ey - sy) ** 2)
    if name == "HC":
        patches = make_helix_cylinder(helix_config, length)
    elif name == "HW":
        patches = make_helix_wave(helix_config, length)
    elif name == "E":
        patches = make_sheet(sheet_config, length)
    elif name == "T":
        patches = make_turn(turn_config, length)
    else:
        raise ValueError("secondary structure element name must be one of HC, HW, E, T")
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
