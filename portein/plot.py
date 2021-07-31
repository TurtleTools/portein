import matplotlib.pyplot as plt
import numpy as np
import prody as pd
from matplotlib import patches as m_patches
from matplotlib import transforms as m_transforms
import typing as ty
from pathlib import Path

from portein.config import HelixConfig, TurnConfig, SheetConfig, PorteinConfig
from portein.rotate import find_best_projection, rotate_to_maximize_bb_height

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
    pdb: ty.Union[str, Path, pd.AtomGroup], config: PorteinConfig, height=12, width=None, ax=None
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
    if type(pdb) == str or type(pdb) == Path:
        structure = pd.parsePDB(pdb)
    else:
        structure = pdb
    dssp_file = pd.execDSSP(pdb)
    structure = pd.parseDSSP(dssp_file, structure)
    structure_alpha = structure.select("calpha")
    coords = structure_alpha.getCoords()
    for i in range(20):
        coords = find_best_projection(coords)
    coords = rotate_to_maximize_bb_height(coords[:, :2])
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
