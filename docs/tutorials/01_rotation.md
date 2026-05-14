---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Optimal 2D orientation

Portein uses some linear algebra (for [Optimal rotation of 3D model for 2D projection](https://stackoverflow.com/a/2970340) and [Rotating an object to maximize bounding box height](https://stackoverflow.com/a/47844156)) to find the best 2D projection for the input protein's 3D coordinates.

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np

import portein
pdb = portein.read_structure("../_data/7lc2.pdb")
old_coords = pdb[pdb.atom_name == "CA"].coord

pdb_oriented = portein.rotate_protein(pdb)
new_coords = pdb_oriented[pdb_oriented.atom_name == "CA"].coord
```

`portein.find_size` returns a figure size that preserves the structure's
aspect ratio given either a width or a height.

```{code-cell} ipython3
old_width, old_height = portein.find_size(old_coords, height=5)
new_width, new_height = portein.find_size(new_coords, height=5)

fig, ax = plt.subplots(
    1,
    2,
    figsize=(old_width + new_width, new_height),
    gridspec_kw={"width_ratios": [old_width, new_width]},
)
for axis, coords, title in [
    (ax[0], old_coords, "Before rotation"),
    (ax[1], new_coords, "After rotation"),
]:
    axis.plot(coords[:, 0], coords[:, 1], "-", c="black")
    axis.scatter(
        coords[:, 0],
        coords[:, 1],
        c=np.arange(coords.shape[0]),
        s=50,
        cmap="Blues",
        edgecolors="gray",
    )
    axis.set_title(title, fontsize=16)
    axis.axis("off")
plt.tight_layout()
```

## From the command line

Save an oriented version of any structure with:

```sh
portein rotate 7lc2
```

It writes `7lc2_rotated.cif` next to the input.
