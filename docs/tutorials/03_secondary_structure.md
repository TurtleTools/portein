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

# Secondary-structure topology diagrams

Portein runs DSSP to assign each residue to a secondary-structure element
(SSE), then draws each element with a parameterized glyph:

* helices as either waves or cylinders (toggle via
  {class}`portein.HelixConfig`'s `as_cylinder`)
* β-strands as arrows
* turns as arcs joined by circles

Adapted from [this gist](https://gist.github.com/JoaoRodrigues/f9906b343d3acb38e39f2b982b02ecb0).

```{code-cell} ipython3
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt

import portein

output_dir = Path(tempfile.mkdtemp())
```

## Default topology

```{code-cell} ipython3
protein_config = portein.ProteinConfig(
    pdb_file="../_data/7lc2.pdb",
    rotate=True,
    width=1000,
    output_prefix=str(output_dir / "secondary_structure"),
)
pss = portein.SecondaryStructure(
    protein_config=protein_config,
    helix_config=portein.HelixConfig(),
    sheet_config=portein.SheetConfig(),
    turn_config=portein.TurnConfig(),
    dpi=100,
)
pss.run();
```

## Highlight specific residues

`SecondaryStructure.run()` returns a Matplotlib Axes — overlay any
matplotlib primitive on top.

```{code-cell} ipython3
ax = pss.run()
ax.set_title("Portrait of PDB ID: 7lc2", fontsize=20)
highlight_residues = [30, 35, 25, 10, 11, 12, 13, 14, 15]
ax.scatter(
    pss.coords[highlight_residues, 0],
    pss.coords[highlight_residues, 1],
    color="red",
    s=100,
    edgecolor="black",
    linewidth=2,
);
```

## Linear diagram

Pass `linear=True` for a single-row strip — useful as a header band on top
of a sequence alignment or a per-residue plot.

```{code-cell} ipython3
fig, ax = plt.subplots(1, figsize=(50, 1))
pss.run(ax=ax, linear=True);
```

## From the command line

```sh
portein secondary 7lc2
```

`-h`, `-s`, `-t` accept YAML config files for helix, sheet, and turn
parameters. See `configs/` in the repo for examples.
