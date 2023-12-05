'''skm_tools resources '''

import os

# visualisation
STYLE_XML = "skm-default-styles.xml"
PSS_DEFAULT_STYLE = "PSS-default"
CKN_DEFAULT_STYLE = "CKN-default"
BUILTIN_STYLES = ['pss', 'ckn']

def get_style_xml_path():
    '''file path on disk '''
    module_path = os.path.abspath(__file__)
    return os.path.join(os.path.dirname(module_path), STYLE_XML)
