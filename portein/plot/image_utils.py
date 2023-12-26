import numpy as np
from PIL import Image
import typing as ty


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


def put_alpha(img: Image, transparency: float):
    im2 = img.copy()
    im2.putalpha(int(255 * (1 - transparency)))
    img.paste(im2, img)
    return img


def alpha_blending(color, alpha: float):
    return tuple((1.0 - alpha) + np.array(color) * alpha)


def find_size(
    points, width: ty.Optional[float] = None, height: ty.Optional[float] = None
):
    """
    Find figure size given coordinates and size of one dimension.
    Calculates height from width if height=None and width from height if width=None

    Parameters
    ----------
    points
    width
    height

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
