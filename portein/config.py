from dataclasses import dataclass


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

    as_cylinder: bool
    cylinder_ellipse_length: float
    cylinder_ellipse_height: float
    cylinder_rectangle_height: float
    wave_arc_width: float
    wave_arc_height: float
    wave_arc_length: float
    outline_width: float
    outline_color: str
    color: str
    opacity: float

    @classmethod
    def default(cls, width=1.0, outline_width=3.0, as_cylinder=False):
        return cls(
            as_cylinder,
            cylinder_ellipse_length=width / 2,
            cylinder_ellipse_height=width - 0.001,
            cylinder_rectangle_height=width,
            wave_arc_width=3.0,
            wave_arc_height=width,
            wave_arc_length=0.5,
            outline_width=outline_width,
            outline_color="#5c6887",
            color="lightsteelblue",
            opacity=1.0,
        )


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

    thickness_factor: float
    tail_height: float
    head_height: float
    outline_width: float
    outline_color: str
    color: str
    opacity: float

    @classmethod
    def default(cls, width=1.0, outline_width=3.0):
        return cls(
            thickness_factor=1,
            tail_height=width,
            head_height=width * 2,
            outline_width=outline_width,
            outline_color="#5c6887",
            color="#999FD0",
            opacity=1.0,
        )


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

    thickness_factor: float
    height: float
    circle_radius: float
    circle_color: str
    arc_width: float
    arc_color: str
    opacity: float

    @classmethod
    def default(cls, width=1.0, outline_width=3.0):
        return cls(
            thickness_factor=1 / 2,
            height=width / 2,
            circle_radius=0.2,
            circle_color="#d1d6e3",
            arc_width=outline_width,
            arc_color="#d1d6e3",
            opacity=0.8,
        )


@dataclass
class PorteinConfig:
    """
    Config for controlling appearance of secondary structural elements

    See docs for `HelixConfig`, `SheetConfig` and `TurnConfig`
    """

    helix: HelixConfig
    sheet: SheetConfig
    turn: TurnConfig

    @classmethod
    def default(cls, width=1.0, outline_width=3.0):
        """
        Generate the default config

        Parameters
        ----------
        width
            Controls the width of each element
        outline_width
            Controls the line width of the outline for helices and sheets,
            and the arc line width for turns

        Returns
        -------
        PorteinConfig object
        """
        return cls(
            helix=HelixConfig.default(width, outline_width),
            sheet=SheetConfig.default(width, outline_width),
            turn=TurnConfig.default(width, outline_width),
        )
