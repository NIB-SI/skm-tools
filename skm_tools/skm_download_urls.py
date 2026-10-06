'''URLs for downloading networks from skm.nib.si'''

# CKN v2.0.1 (on the dev server; not yet on production as of 2026-10-06)
CKN_EDGE_URL = 'https://skm.nib.si/downloads/ckn/v2.0.1-2026.10/edges'
CKN_NODE_URL = 'https://skm.nib.si/downloads/ckn/v2.0.1-2026.10/nodes'

# PSS network exports (pss-export), public version
# (provisional: not on skm.nib.si until the web app with the new exports goes live)
PSS_REACTION_GRAPH_EDGE_URL = 'https://skm.nib.si/downloads/pss/public/reaction-graph-edges'
PSS_REACTION_GRAPH_NODE_URL = 'https://skm.nib.si/downloads/pss/public/reaction-graph-nodes'
PSS_INTERACTION_NETWORK_EDGE_URL = 'https://skm.nib.si/downloads/pss/public/interaction-network-edges'
PSS_INTERACTION_NETWORK_NODE_URL = 'https://skm.nib.si/downloads/pss/public/interaction-network-nodes'
# e.g. https://skm.nib.si/downloads/pss/public/gene-network-ath-edges
PSS_GENE_NETWORK_EDGE_URL = 'https://skm.nib.si/downloads/pss/public/gene-network-{}-edges'
PSS_GENE_NETWORK_NODE_URL = 'https://skm.nib.si/downloads/pss/public/gene-network-{}-nodes'

# Translation files
# e.g. https://skm.nib.si/downloads/translations/parm
GENE_TRANSLATION_URL = 'https://skm.nib.si/downloads/translations/{}'