'''Images and charts for network nodes: match image files to nodes, make unique copies, and
build chart definitions.

These helpers don't need Cytoscape; :mod:`skm_tools.cytoscape_utils` shows the results on the
nodes (:func:`skm_tools.cytoscape_utils.load_node_images`,
:func:`skm_tools.cytoscape_utils.show_node_images`).
'''

import logging
import re
import shutil
import uuid
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)


def node_file_key(name):
    '''Default key for matching file names to node names (see :func:`match_files_to_nodes`).

    Replaces the characters that can't (or shouldn't) be in file names with ``_``, e.g.
    ``"WRKY33[fc00166]"`` -> ``"WRKY33_fc00166_"`` and ``"11-/12-OH-JA"`` ->
    ``"11-_12-OH-JA"``; gene ids stay as they are.

    Parameters
    ----------
    name : str
        Node name, or file name without extension.

    Returns
    -------
    str

    Examples
    --------
    >>> node_file_key("WRKY33[fc00166]")
    'WRKY33_fc00166_'
    '''
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(name))


def match_files_to_nodes(folder, nodes, key=node_file_key, aliases=None, extensions=("png", "svg")):
    '''Match image files (e.g. one plot per gene or metabolite) to nodes by name.

    Parameters
    ----------
    folder : str or pathlib.Path
        Folder with the files.
    nodes : iterable
        Node names, e.g. ``g.nodes()`` or ``py4cytoscape.get_all_nodes()``.
    key : callable
        Applied to node names and file names (without extension); they match if the keys
        are equal. Default :func:`node_file_key`.
    aliases : dict, optional
        File name (without extension) -> node name, for files not named after their node
        (e.g. ``{"Pro": "Proline accumulation"}``).
    extensions : iterable of str
        File extensions to use (default png and svg).

    Returns
    -------
    dict
        Node -> file path. Nodes with the same key get the same file. Files that match no
        node are logged (INFO).
    '''
    aliases = aliases or {}
    nodes_by_key = defaultdict(list)
    for n in nodes:
        nodes_by_key[key(n)].append(n)

    matched, unmatched = {}, []
    for ext in extensions:
        for path in sorted(Path(folder).glob(f"*.{ext}")):
            k = key(aliases.get(path.stem, path.stem))
            if k in nodes_by_key:
                matched.update({n: path for n in nodes_by_key[k]})
            else:
                unmatched.append(path.name)

    if unmatched:
        logger.info("%d files in %s match no node: %s", len(unmatched), folder, unmatched)
    return matched


def unique_image_copies(images, folder, to_png=False):
    '''Copy images to new, unique file names, optionally converting SVG to PNG.

    Cytoscape caches images, and doesn't always notice that a file is a different one (e.g.
    the same file name in another folder) or has changed: it may show an old image instead.
    Every call makes new copies, with names that no other file has had (``<name>_<random
    id>``), so Cytoscape has to load them again. Old copies are not deleted.

    Parameters
    ----------
    images : dict
        Node -> image path.
    folder : str or pathlib.Path
        Folder for the copies (created if needed).
    to_png : bool
        Convert SVG images to PNG. Needs ``cairosvg`` (``pip install cairosvg``).

    Returns
    -------
    dict
        Node -> path of the copy. Nodes with the same image share one copy.
    '''
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    copy_of = {}
    for path in {Path(p) for p in images.values()}:
        convert = to_png and path.suffix.lower() == ".svg"
        copy = folder / f"{path.stem}_{uuid.uuid4().hex[:12]}{'.png' if convert else path.suffix}"
        if convert:
            import cairosvg
            cairosvg.svg2png(url=str(path), write_to=str(copy))
        else:
            shutil.copyfile(path, copy)
        copy_of[path] = copy
    return {node: copy_of[Path(p)] for node, p in images.items()}


def chart_column(df, columns, colors, chart="barchart", value_range=None, labels=None,
                 na_value=0, **options):
    '''Chart definitions for each row of `df`, to show as node charts in Cytoscape.

    Builds `enhancedGraphics <https://apps.cytoscape.org/apps/enhancedgraphics>`_ chart
    strings (e.g. ``barchart: colorlist="..." valuelist="..."``); load them as a node table
    column and show them like images, with :func:`skm_tools.cytoscape_utils.show_node_images`.
    Needs the enhancedGraphics app in Cytoscape.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per node.
    columns : list of str
        Columns of `df` to chart, in order.
    colors : str or list of str
        Hex colour per column, or one for all.
    chart : str
        enhancedGraphics chart type, e.g. ``"barchart"`` (default), ``"linechart"``,
        ``"heatstripchart"``, ``"piechart"``.
    value_range : tuple of float, optional
        (min, max) of the value axis (default: per chart, from its values).
    labels : list of str, optional
        Label per column (default: the column names). enhancedGraphics splits lists on
        commas and can't escape them, so labels can't contain ``,`` or ``"``.
    na_value : float
        Value for missing data (default 0).
    **options
        Further enhancedGraphics options, e.g. ``separation=2``, ``ybase=1``.

    Returns
    -------
    pandas.Series
        Chart string per row of `df`.

    Raises
    ------
    ValueError
        If a label contains ``,`` or ``"``.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({"t1": [1.5], "t2": [-0.5]}, index=["AT2G38470"])
    >>> chart_column(df, ["t1", "t2"], "#E41A1C", value_range=(-2, 2)).iloc[0]
    'barchart: colorlist="#E41A1C,#E41A1C" valuelist="1.5,-0.5" labellist="t1,t2" range="-2,2"'
    '''
    if isinstance(colors, str):
        colors = [colors] * len(columns)
    labels = [str(label) for label in (labels or columns)]
    bad = [label for label in labels if "," in label or '"' in label]
    if bad:
        raise ValueError(f"Chart labels can't contain ',' or '\"' (enhancedGraphics can't escape "
                         f"them): {bad}. Pass other `labels`.")

    fixed = f'colorlist="{",".join(colors)}"'
    rest = f'labellist="{",".join(labels)}"'
    if value_range is not None:
        rest += f' range="{value_range[0]:g},{value_range[1]:g}"'
    rest += "".join(f" {k}={v}" for k, v in options.items())

    values = df[columns].astype(float).fillna(na_value)
    return values.apply(
        lambda row: f'{chart}: {fixed} valuelist="{",".join(f"{v:g}" for v in row)}" {rest}',
        axis=1,
    )
