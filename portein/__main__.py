from portein import plot_portrait, PorteinConfig
import matplotlib.pyplot as plt
import click
import warnings

warnings.filterwarnings("ignore")


@click.command("portein", context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("pdb", metavar="<pdb>")
@click.argument("output", metavar="<output>")
def main(pdb, output):
    """
    Save a 2D portrait of a 3D protein structure

    \b
    <pdb> is a PDB ID or PDB file path
    <output> is the output image file path with extension
    """
    plot_portrait(pdb, PorteinConfig.default(), 12)
    plt.savefig(output)


if __name__ == "__main__":
    main()
