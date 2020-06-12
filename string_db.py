"""
SQL database handing based on code found on:
https://github.com/scanisius/discover-notebooks/blob/master/lib/nbsupport/stringdb.py
"""

import urllib.request
from pathlib import Path
import gzip
from contextlib import closing
import sqlite3
import re

script_path = Path(__file__).parent.absolute()

# String-DB base URL and database names to be imported. It's important that 'protein.links.detailed' comes first and is
# consequently being parsed first as well, otherwise exception might occur when creating the SQL database.
base_fns = ['protein.links.detailed.v11.0.txt.gz',
            'protein.aliases.v11.0.txt.gz']
base_url = 'https://stringdb-static.org/download/'


# Return True if string represents a number.
def is_str_number(string):
    try:
        int(string)
        return True
    except ValueError:
        return False


# Print file download progress.
def download_progress(chunk_num, chunk_size, total_size):
    downloaded = chunk_num * chunk_size
    progress = int((downloaded/total_size)*100)
    print(f'Download progress: {str(progress)}%', end='\r')


# Use STRING-DB file header to identify the file that's being read. Return True if file contains detailed PPI links, or
# False when file's identified to contain protein aliases. When neither holds true, an exception's raised.
def contains_ppi_links(header):
    # List of all the column names expected to be read in each case.
    pp_links_header = ['protein1', 'protein2', 'combined_score']
    p_aliases_header = ['string_protein_id', 'alias', 'source']

    # Make a list of the column names contained in the header given.
    header = header.replace('## ', '').split()

    if set(pp_links_header) <= set(header):
        return True
    elif set(p_aliases_header) <= set(header):
        return False
    else:
        raise Exception('Parse error: unexpected file header.')


# Find and import STRING-DB database files found under ../data/. User has the option to provide a list of NCBI taxonomy
# IDs to only import data related to specific species. If such data don't exist, function attempts to download them.
def import_txtgz(db_cur, ncbi_taxonomy_ids=None):
    if not ncbi_taxonomy_ids:
        ncbi_taxonomy_ids = ['']

    # List of all the filenames that need to be imported or downloaded.
    txtgz_fns = []
    for each in ncbi_taxonomy_ids:
        txtgz_fns.extend(
            map(lambda x: '.'.join(filter(None, [each, x])), base_fns))

    url_opener = urllib.request.build_opener()
    url_opener.addheaders = [('User-agent', 'Mozilla/5.0')]
    urllib.request.install_opener(url_opener)
    for fn in txtgz_fns:
        txtgz_path = Path.joinpath(script_path, 'data', fn)

        # If file does not exist, attempt to download it...
        if not txtgz_path.exists():
            ext = ''
            if is_str_number(fn.split('.')[0]):
                ext = '.'.join(fn.split('.')[1:-2]) + '/'
            url = base_url + ext + fn
            # TODO: Check that download wasn't interrupted or raise exception.
            urllib.request.urlretrieve(url, txtgz_path,
                                       reporthook=download_progress)
            print(f'{fn} downloaded.')

        # ...and then import all into the SQL database.
        with gzip.open(txtgz_path, 'rt') as f:
            # Extract header, then use it to identify file contents.
            hdr = f.readline()
            if contains_ppi_links(hdr):
                # (a) This file contains PPI information.
                proteins = set()
                for line in f:
                    records = line.split(' ')
                    proteins.add(records[0])

                    _, protein_1 = records[0].split('.')
                    _, protein_2 = records[1].split('.')
                    combined_score = int(records[9])

                    db_cur.execute(
                        "INSERT INTO pp_links VALUES (?, ?, ?)",
                        (protein_1, protein_2, combined_score)
                    )

                for protein in proteins:
                    species, ensembl_id = protein.split('.')

                    db_cur.execute(
                        "INSERT INTO proteins VALUES (?, ?)",
                        (ensembl_id, species)
                    )

            else:
                # (b) This file contains protein aliases.
                for line in f:
                    records = re.split(r'\t+', line)
                    _, ensembl_id = records[0].split('.')
                    alias, sources = records[1:3]

                    db_cur.execute(
                        "INSERT INTO p_aliases VALUES (?, ?, ?)",
                        (ensembl_id, alias, sources)
                    )


# The 'proteins' table holds Ensembl-IDs and associated species information. The 'pp_links' table holds PPI data based
# on STRING's combined score. The 'p_aliases' table associates aliases from multiple sources with their Ensembl-ID. No
# need to issue a COMMIT command, sqlite3 operates in auto-commit mode.
def sqlite_create_tables(db_cur):

    db_cur.executescript("""
        CREATE TABLE IF NOT EXISTS proteins
        (
            ensembl_id TEXT NOT NULL UNIQUE,
            species    TEXT NOT NULL,
            PRIMARY KEY (ensembl_id)
        );

        CREATE TABLE IF NOT EXISTS pp_links
        (
            ensembl_id_1   TEXT    NOT NULL,
            ensembl_id_2   TEXT    NOT NULL,
            combined_score INTEGER NOT NULL,
            PRIMARY KEY (ensembl_id_1, ensembl_id_2),
            FOREIGN KEY (ensembl_id_1)
                REFERENCES proteins (ensembl_id)
        );

        CREATE TABLE IF NOT EXISTS p_aliases
        (
            ensembl_id TEXT NOT NULL,
            alias      TEXT NOT NULL,
            sources    TEXT NOT NULL,
            PRIMARY KEY (ensembl_id, alias, sources),
            FOREIGN KEY (ensembl_id)
                REFERENCES proteins (ensembl_id)
        );""")


# Append StringDB data into SQL database, if one doesn't already exist.
# TODO: Append new data when database already present.
def create_string_db(ncbi_taxonomy_ids=None):
    db_path = Path.joinpath(script_path, 'string.db')
    # If 'string.db' already exists, terminate.
    if db_path.exists():
        return
    else:
        print('Creating STRING database...')

    with closing(sqlite3.connect('string.db')) as conn:
        with conn:  # auto-commit
            with closing(conn.cursor()) as cur:
                # First create the tables needed...
                sqlite_create_tables(cur)
                # ...then append to those tables.
                import_txtgz(cur, ncbi_taxonomy_ids)

    print('STRING database built successfully!')
