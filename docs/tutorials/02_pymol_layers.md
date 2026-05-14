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

# Layered PyMOL renders

Portein composes [PyMOL](https://github.com/schrodinger/pymol-open-source)
ray-traced images by stacking layers — each layer is its own PyMOL
representation rendered independently, then alpha-composited together
with user-controlled transparencies.

Any setting valid for `pymol.cmd.set(...)` can be passed via
{class}`portein.PymolConfig`'s `pymol_settings` dict.

```{code-cell} ipython3
import tempfile
from pathlib import Path

import yaml
from biotite import structure as struct
from biotite.structure import io as bio
from IPython.display import Image

import portein
```

The package ships some default PyMOL settings under `configs/`:

```{code-cell} ipython3
with open("../../configs/pymol_settings.yaml") as f:
    pymol_settings = yaml.safe_load(f)
pymol_settings
```

```{code-cell} ipython3
output_dir = Path(tempfile.mkdtemp())
```

## Simple cartoon

A single-layer cartoon, colored by chain. `rotate=True` orients the structure
before rendering.

```{code-cell} ipython3
protein_config = portein.ProteinConfig(
    pdb_file="../_data/7lc2.pdb",
    rotate=True,
    width=1000,
    chain_colormap="Set3",
    output_prefix=str(output_dir / "7lc2_simple"),
)
pymol = portein.Pymol(
    protein=protein_config,
    layers=[portein.PymolConfig(representation="cartoon", pymol_settings=pymol_settings)],
)
image_file = pymol.run()
Image(image_file)
```

## Four-layer composite

Layer the surface (50% opacity), the cartoon, sticks for highlighted
residues, and sticks colored green for the bound ligand. The
`selection="highlight"` keyword picks up the residues passed via
`highlight_residues` on the `ProteinConfig`.

```{code-cell} ipython3
protein_config = portein.ProteinConfig(
    pdb_file="../_data/7lc2.pdb",
    rotate=True,
    chain_colormap="Set3",
    highlight_residues={
        "A": {"black": [30, 35], "red": list(range(10, 20))},
        "B": {"black": [25], "red": list(range(10, 16))},
    },
    width=1000,
    output_prefix=str(output_dir / "7lc2"),
)
layers = [
    portein.PymolConfig(representation="surface", pymol_settings=pymol_settings, transparency=0.5),
    portein.PymolConfig(representation="cartoon", pymol_settings=pymol_settings),
    portein.PymolConfig(representation="sticks", pymol_settings=pymol_settings, selection="highlight"),
    portein.PymolConfig(
        representation="sticks",
        pymol_settings=pymol_settings,
        selection="resn GNP",
        color="green",
    ),
]
pymol = portein.Pymol(protein=protein_config, layers=layers, buffer=10, combine=True)
image_file = pymol.run()
Image(image_file)
```

## Zoom on a ligand pocket

Pick out the residues near a bound ligand, rotate to maximize the pocket's
projection, save only the proximal chains, and render with a translucent
surface.

```{code-cell} ipython3
import numpy as np

pdb = portein.read_structure("../_data/7lc2.pdb")
ligand = pdb[(pdb.chain_id == "A") & (pdb.res_name == "GNP")]
mask = struct.CellList(pdb, 6).get_atoms(ligand.coord, 6, as_mask=True).any(axis=0)
ligand_pocket = pdb[mask]
proximal_chains = struct.get_chains(ligand_pocket)

# Orient to best-show the pocket
rotation, translation = portein.get_best_transformation(ligand_pocket.coord.astype(np.float64))
pdb_oriented = struct.rotate(struct.translate(pdb, translation), rotation)

# Save only the proximal chains
pocket_path = output_dir / "7lc2_ligand.pdb"
bio.save_structure(str(pocket_path), pdb_oriented[np.isin(pdb_oriented.chain_id, proximal_chains)])

protein_config = portein.ProteinConfig(
    pdb_file=str(pocket_path),
    rotate=False,
    chain_colormap="white",
    width=1000,
    output_prefix=str(output_dir / "7lc2_pocket"),
)
layers = [
    portein.PymolConfig(representation="surface", pymol_settings=pymol_settings, transparency=0.3),
    portein.PymolConfig(representation="cartoon", pymol_settings=pymol_settings),
    portein.PymolConfig(
        representation="sticks",
        pymol_settings=pymol_settings,
        selection="(chain A and resn GNP)",
        color="green",
    ),
]
pymol = portein.Pymol(protein=protein_config, layers=layers, combine=True)
image_file = pymol.run()
Image(image_file)
```

## Layer compositing: 2D decals vs depth-interleaved (`combine`)

By default, `portein.Pymol` ray-traces each layer separately and alpha-composites the PNGs in PIL. That's the **2D decals** behavior — a later layer always sits "on top" of an earlier layer at every pixel where it isn't transparent, regardless of whether the underlying 3D geometry would be in front or behind. It's the right tool when you want a translucent surface applied uniformly across the whole rendered image, or when your layers' geometries don't overlap in 3D.

When two layers' geometries *do* overlap and you want PyMOL's z-buffer to interleave them per pixel — the way PyMOL's native `super` command composites two superposed structures — set `combine=True`. All layers go into one PyMOL scene and a single `cmd.ray()` produces the final image. Per-layer `transparency` then becomes a PyMOL per-selection `cartoon_transparency` / `transparency` / `stick_transparency` setting (depth-aware) instead of a flat PIL alpha.

```{code-cell} ipython3
from IPython.display import Image as IpyImage, display

protein_config = portein.ProteinConfig(
    pdb_file="../_data/7lc2.pdb",
    rotate=True,
    width=600,
    chain_colormap="Set3",
    output_prefix=str(output_dir / "combine_demo"),
)
layers = [
    portein.PymolConfig(representation="surface", pymol_settings=pymol_settings, transparency=0.5),
    portein.PymolConfig(representation="cartoon", pymol_settings=pymol_settings),
]
```

```{code-cell} ipython3
# Default (combine=False): the cartoon PNG sits entirely on top of the
# translucent surface PNG.
protein_config.output_prefix = str(output_dir / "decals")
decals_image = portein.Pymol(protein=protein_config, layers=layers).run()
display(IpyImage(decals_image))
```

```{code-cell} ipython3
# combine=True: one ray-trace. The translucent surface and the cartoon
# share the z-buffer — the cartoon appears "behind frosted glass" wherever
# the surface is in front, and unmuted where the cartoon protrudes.
protein_config.output_prefix = str(output_dir / "combined")
combined_image = portein.Pymol(protein=protein_config, layers=layers, combine=True).run()
display(IpyImage(combined_image))
```

When to pick which:

| Want | Use |
|---|---|
| Uniform 2D translucency over each whole layer | `combine=False` (default) |
| Depth-correct interleaving of overlapping geometries | `combine=True` |
| Two superposed structures rendered like PyMOL's `super` | `combine=True` |
| Surface + sticks where the ligand should occlude correctly | `combine=True` |
| Highlight layer drawn boldly on top regardless of depth | `combine=False` |

## From the command line

Save a YAML with the protein config:

```yaml
pdb_file: 7lc2
rotate: true
width: 1000
chain_colormap: Set3
output_prefix: examples/7lc2_simple
```

Then run:

```sh
portein pymol protein.yaml
portein pymol protein.yaml pymol_layers.yaml --buffer 10
```
