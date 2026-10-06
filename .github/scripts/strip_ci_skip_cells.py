"""Write a copy of a notebook without its '# ci-skip: ...' cells.

Used to run the tutorials headlessly in CI, without the cells that need a running
Cytoscape. The tracked notebook is unchanged, so it still shows them.

    python .github/scripts/strip_ci_skip_cells.py tutorials/tutorial-pss.ipynb tutorials/.ci-tutorial-pss.ipynb
"""
import sys

import nbformat

MARKER = "# ci-skip:"


def main(src_path, dst_path):
    nb = nbformat.read(src_path, as_version=4)
    nb.cells = [
        cell for cell in nb.cells
        if not (cell.cell_type == "code" and cell.source.lstrip().startswith(MARKER))
    ]
    nbformat.write(nb, dst_path)


if __name__ == "__main__":
    main(*sys.argv[1:3])
