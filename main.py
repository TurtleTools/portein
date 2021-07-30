import sys
from portein import plot_portrait, PorteinConfig


def main():
    _, pdb_id, filename = sys.argv
    fig, _, _ = plot_portrait(pdb_id, PorteinConfig.default(1), 12)
    fig.savefig(filename)


if __name__ == "__main__":
    main()
