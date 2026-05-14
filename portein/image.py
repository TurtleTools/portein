"""Public image helpers — operations users compose around the renderers."""

from __future__ import annotations

import numpy as np
from PIL import Image


def crop_to_content(
    img: Image.Image | np.ndarray,
    bg: tuple[int, ...] = (255, 255, 255, 255),
    *,
    treat_transparent_as_bg: bool = True,
) -> Image.Image:
    """Crop uniform-color borders from a rendered image.

    Removes outer rows and columns whose pixels are all background. Useful
    after PyMOL/Illustrate renders to tighten the bounding box around the
    molecule.

    Parameters
    ----------
    img
        Input image as a PIL Image or numpy array (H, W, 3) / (H, W, 4).
    bg
        Background pixel value to strip. Length must match the image's
        channel count. Defaults to opaque white in RGBA, which also matches
        plain white in RGB after truncation.
    treat_transparent_as_bg
        If True and the image has an alpha channel, pixels with alpha=0
        are also treated as background regardless of RGB value. PyMOL
        renders against a transparent background commonly produce such
        pixels around the molecule.

    Returns
    -------
    Cropped PIL Image. If every pixel is background, the original image
    is returned unchanged.
    """
    arr = np.asarray(img.convert("RGBA") if isinstance(img, Image.Image) else img)

    bg_arr = np.asarray(bg[: arr.shape[-1]])
    background_mask = np.all(arr[..., : len(bg_arr)] == bg_arr, axis=-1)
    if treat_transparent_as_bg and arr.shape[-1] == 4:
        background_mask |= arr[..., 3] == 0

    keep_rows = np.any(~background_mask, axis=1)
    keep_cols = np.any(~background_mask, axis=0)
    if not (np.any(keep_rows) and np.any(keep_cols)):
        return img if isinstance(img, Image.Image) else Image.fromarray(arr)

    cropped = arr[keep_rows][:, keep_cols]
    return Image.fromarray(cropped)
