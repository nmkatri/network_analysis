import string_db
import pandas as pd
import os
import pathlib
from typing import Union
import sqlite3
import matplotlib.pyplot as plt
import networkx as nx
import time


# Calculate node centralities (degree, eigenvector, betweenness, closeness), average degree for the network, and edge
# betweenness centrality. Return dictionary containing the results and, by default, export them in TXT format too.
def calculate_network_properties(network, export=True):
    # A dictionary of the properties to calculate and the underlying functions to be called.
    properties_mappings = {
        'average_network_degree': lambda net: net.number_of_edges() / net.number_of_nodes(),
        'degree_centrality': nx.degree_centrality,
        'eigenvector_centrality': nx.eigenvector_centrality,
        'betweenness_centrality': nx.betweenness_centrality,
        'closeness_centrality': nx.closeness_centrality,
        'edge_betweenness_centrality': nx.edge_betweenness_centrality
    }

    # Calculate all properties and save the results in a dictionary for better readability.
    results = {}
    for prop, func in properties_mappings.items():
        results[prop] = func(network)

    # Export the results in TXT format (if user chooses to).
    if export:
        for prop, data in results.items():
            # Append time of export to avoid unwanted overwriting.
            time_str = time.strftime('%d%b%Y-%H%M%S')
            with open(f'output-{prop}-{time_str}.txt', 'w') as f:
                if isinstance(data, dict):
                    # (a) We are exporting a dictionary...
                    # Before starting, sort dictionary in descending order!
                    data = dict(sorted(data.items(),
                                       key=lambda x: x[1], reverse=True))
                    for key, value in data.items():
                        # (i) If key is of type 'str', it should be a node...
                        if isinstance(key, str):
                            output_str = f'{key} {value}\n'
                        # (ii) ...or it's a tuple of nodes (making an edge).
                        else:
                            output_str = f'{key[0]} {key[1]} {value}\n'
                        f.write(output_str)
                else:
                    # (b) ...otherwise it's just a float.
                    # TODO: That's a single float in a file... meh!
                    f.write(f'{data}\n')

    return results


# Import every spreadsheets found in the working path as a list of DataFrames.
def import_genes():
    script_path = pathlib.Path(__file__).parent.absolute()

    # PyCharm BUG: pathlib.Path objects not considered 'PathLike'
    fn: Union[pathlib.Path, os.PathLike]
    # Retrieve all Excel files, but skip any starting with 'output' or 'temp'.
    sheet_paths = [fn for fn in script_path.glob('*.xlsx')
                   if not os.path.basename(fn).startswith(('output', 'temp'))]

    _sheet_list = []
    for path in sheet_paths:
        _sheet = pd.read_excel(path)
        _sheet_list.append(_sheet)
    _num_of_sheets = len(_sheet_list)

    return _sheet_list, _num_of_sheets


if __name__ == '__main__':
    # Build STRING database for species-ID: 10090 (mouse).
    string_db.create_string_db(['10090'])

    # Import mouse genes from spreadsheet and create list.
    sheet_list = import_genes()
    sheet = sheet_list[0][0]
    genes = list(sheet.Mouse_gene)
    genes = list(map(str.strip, genes))

    # Use STRING-DB to build and expand the PPI network:
    # (a) First match gene names with their ENSEMBL-IDs.
    conn = sqlite3.connect('string.db')
    sql_get_ids = """
        SELECT ensembl_id
        FROM p_aliases
        WHERE alias IN ({seq})
        """.format(seq=','.join(['?'] * len(genes)))

    cursor = conn.execute(sql_get_ids, genes)
    gene_ids = [x[0] for x in cursor.fetchall()]

    # (b) Then get PPI interactions between given genes.
    sql_get_ppi_1 = """
        SELECT ensembl_id_1, ensembl_id_2
        FROM pp_links
        WHERE ensembl_id_1 IN ({seq})
          AND ensembl_id_2 IN ({seq})
          AND combined_score > 400
        """.format(seq=','.join(['?'] * len(genes)))

    cursor = conn.execute(sql_get_ppi_1, 2*gene_ids)
    ppi_initial = cursor.fetchall()

    # (c) Expand PPI interactions by one neighbourhood.
    sql_get_ppi_2 = """
        SELECT ensembl_id_1, ensembl_id_2
        FROM pp_links
        WHERE ensembl_id_1 IN ({seq})
          AND combined_score > 400
        """.format(seq=','.join(['?'] * len(genes)))

    cursor = conn.execute(sql_get_ppi_2, gene_ids)
    ppi_expanded = cursor.fetchall()

    # Create graphs for initial and expanded networks, then calculate their properties.
    nw0 = nx.Graph()
    nw0.add_edges_from(ppi_initial)
    props_nw0 = calculate_network_properties(nw0)

    nw1 = nx.Graph()
    nw1.add_edges_from(ppi_expanded)
    props_nw1 = calculate_network_properties(nw1)
