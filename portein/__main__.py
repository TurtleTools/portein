import portein
import matplotlib.pyplot as plt
import typer
from pathlib import Path
import warnings
from typing_extensions import Annotated
from typing import Optional

warnings.filterwarnings("ignore")
portein.compile_numba_functions()


app = typer.Typer(name="portein")


@app.command()
def pymol(
    protein: Annotated[
        Path,
        typer.Argument(
            ..., help="Path to protein config file", exists=True, show_default=False
        ),
    ],
    config: Annotated[
        Optional[Path],
        typer.Argument(..., help="Path to pymol config file", show_default=False),
    ] = None,
    buffer: Annotated[
        Optional[float],
        typer.Option(
            help="Buffer around the molecule in Angstroms", show_default=False
        ),
    ] = None,
):
    """
    Render protein using Pymol
    """
    pymol_class = portein.Pymol.from_yaml(protein, config, buffer)
    image_file = pymol_class.run()
    print(f"Saved image to {image_file}")


@app.command()
def illustrate(
    protein: Annotated[
        Path,
        typer.Argument(
            ..., help="Path to protein config file", exists=True, show_default=False
        ),
    ],
    config: Annotated[
        Optional[Path],
        typer.Argument(..., help="Path to illustrate config file", show_default=False),
    ] = None,
):
    """
    Render protein using illustrate
    """
    illustrate_class = portein.Illustrate.from_yaml(protein, config)
    image_file = illustrate_class.run()
    print(f"Saved image to {image_file}")


@app.command()
def rotate(
    protein: Annotated[
        str,
        typer.Argument(
            ..., help="PDB ID or path to PDB file", exists=True, show_default=False
        ),
    ],
    output_prefix: Annotated[
        str,
        typer.Option(
            "--output-prefix, -o",
            help="Prefix of output file (default: same as input file)",
            show_default=False,
        ),
    ] = None,
):
    """
    Orient protein and save it to a new PDB file
    """

    protein = portein.ProteinConfig(protein, rotate=True, output_prefix=output_prefix)
    print(f"Saved protein to {protein.pdb_file}")


@app.command()
def secondary(
    protein: Annotated[
        str,
        typer.Argument(
            ..., help="PDB ID or path to PDB file", exists=True, show_default=False
        ),
    ],
    norotate: Annotated[
        bool,
        typer.Option(
            "--norotate",
            "-n",
            help="Do not re-orient the protein",
        ),
    ] = False,
    output_prefix: Annotated[
        str,
        typer.Option(
            "--output-prefix, -o",
            help="Prefix of output file (default: same as input file)",
            show_default=False,
        ),
    ] = None,
    width: Annotated[
        int,
        typer.Option(
            help="Width of image in pixels",
        ),
    ] = 3000,
    height: Annotated[
        int,
        typer.Option(
            help="Height of image in pixels",
        ),
    ] = None,
    helix: Annotated[
        Optional[Path],
        typer.Option(
            "--helix", "-h", help="Path to helix config file", show_default=False
        ),
    ] = None,
    sheet: Annotated[
        Optional[Path],
        typer.Option(
            "--sheet", "-s", help="Path to sheet config file", show_default=False
        ),
    ] = None,
    turn: Annotated[
        Optional[Path],
        typer.Option(
            "--turn", "-t", help="Path to turn config file", show_default=False
        ),
    ] = None,
    dpi: Annotated[int, typer.Option(help="DPI of image")] = 300,
    keep_files: Annotated[
        bool,
        typer.Option(
            help="Keep intermediate DSSP files",
        ),
    ] = False,
):
    """
    Plot secondary structure topology diagram
    """
    protein = portein.ProteinConfig(protein, rotate=not norotate, output_prefix=output_prefix, width=width, height=height)
    ss = portein.SecondaryStructure(protein, 
                                    portein.HelixConfig.from_yaml(helix),
                                    portein.SheetConfig.from_yaml(sheet),
                                    portein.TurnConfig.from_yaml(turn),
                                    dpi=dpi)
    ss.run()
    image_file = f"{ss.protein_config.output_prefix}_secondary_structure.png"
    plt.savefig(image_file, dpi=dpi, bbox_inches="tight")
    if not keep_files:
        ss.cleanup()
    print(f"Saved image to {image_file}")


if __name__ == "__main__":
    app()
