'''URLs for downloading networks from skm.nib.si, and the default local file names.

The default file names are the names the browser gives the downloads, without the parts
that change between exports (the export date of PSS, the version of the translations).
'''

# CKN v2.0.1 (on the dev server; not yet on production as of 2026-10-07)
CKN_VERSION = 'v2.0.1-2026.10'
CKN_EDGE_URL = f'https://skm.nib.si/downloads/ckn/{CKN_VERSION}/edges'
CKN_NODE_URL = f'https://skm.nib.si/downloads/ckn/{CKN_VERSION}/nodes'
CKN_EDGE_FILE = f'AtCKN-{CKN_VERSION}.tsv.gz'
CKN_NODE_FILE = f'AtCKN-{CKN_VERSION}_node-annot.tsv.gz'

# PSS network exports (pss-export), public, live version: /pss/<access>/<release>/<format>/<file>,
# gene networks /pss/<access>/<release>/gene-network/<species>/<file>
# (not on skm.nib.si until the web app with the new exports goes live)
_PSS = 'https://skm.nib.si/downloads/pss/public/live'
PSS_REACTION_GRAPH_EDGE_URL = f'{_PSS}/reaction-graph/edges'
PSS_REACTION_GRAPH_NODE_URL = f'{_PSS}/reaction-graph/nodes'
PSS_INTERACTION_NETWORK_EDGE_URL = f'{_PSS}/interaction-network/edges'
PSS_INTERACTION_NETWORK_NODE_URL = f'{_PSS}/interaction-network/nodes'
PSS_GENE_NETWORK_EDGE_URL = _PSS + '/gene-network/{}/edges'
PSS_GENE_NETWORK_NODE_URL = _PSS + '/gene-network/{}/nodes'
PSS_REACTION_GRAPH_EDGE_FILE = 'pss-public-reaction-graph-edges.tsv'
PSS_REACTION_GRAPH_NODE_FILE = 'pss-public-reaction-graph-nodes.tsv'
PSS_INTERACTION_NETWORK_EDGE_FILE = 'pss-public-interaction-network-edges.tsv'
PSS_INTERACTION_NETWORK_NODE_FILE = 'pss-public-interaction-network-nodes.tsv'
PSS_GENE_NETWORK_EDGE_FILE = 'pss-public-gene-network-{}-edges.tsv'
PSS_GENE_NETWORK_NODE_FILE = 'pss-public-gene-network-{}-nodes.tsv'

# Translation files, e.g. https://skm.nib.si/downloads/translations/stu
GENE_TRANSLATION_URL = 'https://skm.nib.si/downloads/translations/{}'
GENE_TRANSLATION_FILE = 'translation_ath_to_{}.tsv.gz'
