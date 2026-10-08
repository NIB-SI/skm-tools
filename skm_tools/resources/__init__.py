'''Bundled resources: the SKM Cytoscape styles (SKM and SKM-reactions).'''

import os

# visualisation
STYLE_XML = "skm-styles.xml"

# apply_builtin_style name -> Cytoscape style name
BUILTIN_STYLES = {
    "skm": "SKM",                       # CKN, PSS interaction network, PSS gene networks
    "skm-reactions": "SKM-reactions",   # PSS reaction graph
}


def get_style_xml_path():
    '''Path to the bundled Cytoscape style file (skm-styles.xml).

    Returns
    -------
    str
    '''
    module_path = os.path.abspath(__file__)
    resource_path = os.path.join(os.path.dirname(module_path), STYLE_XML)
    return resource_path
