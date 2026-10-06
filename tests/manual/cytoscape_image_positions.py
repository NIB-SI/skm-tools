"""Manual check of cytoscape_utils.show_node_images positions, in a running Cytoscape.

Creates a network "image positions" with two nodes, "png" and "svg". Each shows five
images (slots 1-5), labelled with the position they should be at: above, below, left,
right and center of the node. Check in Cytoscape that each label is where it says.

    python tests/manual/cytoscape_image_positions.py [folder for the images]
"""

import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

from skm_tools import cytoscape_utils as cu

POSITIONS = ["above", "below", "left", "right", "center"]
COLOURS = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"]

folder = Path(sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="image-positions-"))
folder.mkdir(parents=True, exist_ok=True)

for position, colour in zip(POSITIONS, COLOURS):
    fig = plt.figure(figsize=(1.2, 0.4))
    fig.text(0.5, 0.5, position, ha="center", va="center", fontsize=14, color="white",
             bbox=dict(facecolor=colour, edgecolor="black"))
    fig.savefig(folder / f"{position}.png", transparent=True)
    plt.close(fig)
    (folder / f"{position}.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="120" height="40">'
        f'<rect width="120" height="40" fill="{colour}" stroke="black"/>'
        f'<text x="60" y="27" font-size="18" text-anchor="middle" fill="white">{position}</text>'
        "</svg>"
    )

# py4cytoscape can't create a network without edges
g = nx.Graph([("png", "svg")])
suid = cu.load_network(g, "image positions", collection="image positions")
style = cu.copy_style("default", "image positions", networks=[suid])

# far apart, so the images of the two nodes don't overlap
cu.layout_from_coords(suid, pd.DataFrame({"x": [0, 500], "y": [0, 0]}, index=["png", "svg"]))

for slot, position in enumerate(POSITIONS, start=1):
    images = {"png": folder / f"{position}.png", "svg": folder / f"{position}.svg"}
    column = f"image_{position}"
    cu.load_node_images(images, column, network=suid, unique_dir=folder / "unique")
    cu.show_node_images(style, column, slot=slot, position=position, size=60)

print(f"Network 'image positions' (SUID {suid}); images in {folder}")
