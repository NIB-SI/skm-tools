'''Bundled resources: the SKM Cytoscape styles (PSS-default and CKN-default).'''

import os

# visualisation
STYLE_XML = "skm-default-styles.xml"
PSS_DEFAULT_STYLE = "PSS-default"
CKN_DEFAULT_STYLE = "CKN-default"
BUILTIN_STYLES = ['pss', 'ckn']

def get_style_xml_path():
    '''Path to the bundled Cytoscape style file (skm-default-styles.xml).

    Returns
    -------
    str
    '''
    module_path = os.path.abspath(__file__)
    resource_path = os.path.join(os.path.dirname(module_path), STYLE_XML)
    return resource_path
