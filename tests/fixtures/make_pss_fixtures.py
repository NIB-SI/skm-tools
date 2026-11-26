'''Make the PSS test fixtures: a few reactions of the public PSS exports (as the SKM web app serves them).

    python tests/fixtures/make_pss_fixtures.py <folder with the exports> <export date>

e.g. the folder with pss-2026-11-26-public-reaction-graph-edges.tsv etc. (pss-export, public, no species filter;
gene network ath). Keeps the edges of REACTIONS and the nodes they connect, rows and columns as in the exports.
'''
import sys
from pathlib import Path

import pandas as pd

REACTIONS = ['rx00001', 'rx00002', 'rx00005', 'rx00006', 'rx00007', 'rx00008', 'rx00012', 'rx00019', 'rx00026',
             'rx00096', 'rx00126', 'rx00160',
             'rx01007']     # participant identifiers (reaction graph)
FILES = {'reaction_graph': 'reaction-graph', 'interaction_network': 'interaction-network',
         'gene_network_ath': 'gene-network-ath'}


def read(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def main(folder, export_date):
    out = Path(__file__).parent
    for fixture, name in FILES.items():
        edges = read(Path(folder) / f'pss-{export_date}-public-{name}-edges.tsv')
        nodes = read(Path(folder) / f'pss-{export_date}-public-{name}-nodes.tsv')
        edges = edges[edges['reaction_id'].isin(REACTIONS)]
        nodes = nodes[nodes['id'].isin(set(edges['source']) | set(edges['target']))]
        edges.to_csv(out / f'pss_{fixture}_edges.tsv', sep='\t', index=False)
        nodes.to_csv(out / f'pss_{fixture}_nodes.tsv', sep='\t', index=False)
        print(fixture, len(edges), 'edges', len(nodes), 'nodes')


if __name__ == '__main__':
    main(*sys.argv[1:])
