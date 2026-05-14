# Portein

## Portraits of Proteins

Portein plots 3D proteins according to their best 2D projection (best = greatest visible area), and composes PyMOL ray-traced renders, secondary-structure topology diagrams, and sequence visualizations into publication-quality figures.

[![Documentation](https://img.shields.io/badge/docs-jupyter--book-blue)](https://turtletools.github.io/portein/)

➡️  **Read the docs** at <https://turtletools.github.io/portein/> for usage examples.

## Installation

Fresh env:

```sh
mamba env create -f https://raw.githubusercontent.com/TurtleTools/portein/main/environment-user.yml
mamba activate portein
portein --help
```

Existing env:

```sh
mamba env update -n <your-env> -f https://raw.githubusercontent.com/TurtleTools/portein/main/environment-user.yml
```

Or, if your env already has `pymol-open-source` and `dssp`, just the pip part:

```sh
pip install git+https://github.com/TurtleTools/portein.git
```

### Developer install (editable + tests)

```sh
git clone https://github.com/TurtleTools/portein.git
cd portein
mamba env create -f environment.yml
mamba activate portein
pytest
```

### Minimal install (Python API only)

If you only need the rotation logic and don't want PyMOL or DSSP:

```sh
pip install git+https://github.com/TurtleTools/portein.git
```

Then add the renderers if/when you need them:

```sh
mamba install -c conda-forge pymol-open-source dssp
```
