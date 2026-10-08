'''skm_tools: toolbox for analysis and visualisation of the SKM networks (PSS and CKN).

The Cytoscape modules (:mod:`skm_tools.cytoscape_utils`, :mod:`skm_tools.cytoscape_pdf_utils`)
need the optional extras and are not imported here.

The functions log what they did (e.g. how many nodes a filter removed) with the standard
:mod:`logging` module, under the ``skm_tools`` logger. To see the messages::

    import logging
    logging.basicConfig(level="INFO")
'''

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("skm-tools")
except PackageNotFoundError:  # not installed, e.g. run from the source folder
    __version__ = "unknown"

import skm_tools.pss
import skm_tools.ckn
import skm_tools.annotations
import skm_tools.paths
import skm_tools.neighbors
import skm_tools.experimental_data
import skm_tools.cuts
import skm_tools.translate
import skm_tools.serialize
import skm_tools.persistence
import skm_tools.dinar
