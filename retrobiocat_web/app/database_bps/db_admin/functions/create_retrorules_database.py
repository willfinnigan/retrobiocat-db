import sqlite3
from rdkit import Chem
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator
import os
import numpy as np
from pathlib import Path
import gzip
import shutil
from retrobiocat_web.app.database_bps.db_admin.functions.download_starting_files import download_file

data_folder = str(Path(__file__).parents[3]) + '/retro/data/retrorules'

def make_rr_dir():
    Path(data_folder).mkdir(parents=True, exist_ok=True)

def download_retrorules_files():
    print('DOWNLOADING RETRORULES FILES')
    rr_path = f'{data_folder}/rr.gz'
    rr_unzip = f'{data_folder}/retrorules_rr02_flat_all.tsv'
    download_file("https://retrorules.org/dl/preparsed/rr02/rp3/nohs", rr_path)

    with gzip.open(rr_path, 'rb') as f_in:
        with open(rr_unzip, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def delete_downloaded_files():
    rr_file = f'{data_folder}/rr.gz'
    if os.path.exists(rr_file):
        os.remove(rr_file)

    rr_unzip = f'{data_folder}/retrorules_rr02_flat_all.tsv'
    if os.path.exists(rr_unzip):
        os.remove(rr_unzip)

def delete_old_db():
    db_file = f'{data_folder}/retrorules.db'
    if os.path.exists(db_file):
        os.remove(db_file)

def create_db():
    print('Creating retrorules database')
    conn = sqlite3.connect(f'{data_folder}/retrorules.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE RETRORULES
              ([Rule_ID] text PRIMARY KEY, 
              [Legacy_ID] text, 
              [Reaction_ID] text,
              [Diameter] integer,
              [Rule_order] integer,
              [Rule_SMARTS] text,
              [Substrate_ID] text,
              [Substrate_SMILES] text,
              [Product_IDs] text,
              [Product_SMILES] text,
              [Rule_SMILES] text,
              [Rule_SMARTS_lite] text,
              [Score] real,
              [Score_normalized] real,
              [Reaction_EC_number] text,
              [Reaction_direction] text,
              [Rule_relative_direction] integer,
              [Rule_usage] text)''')
    conn.commit()

def data_to_db():
    print('Loading data into retrorules database')
    df = pd.read_csv(f'{data_folder}/retrorules_rr02_flat_all.tsv', sep='\t')
    conn = sqlite3.connect(f'{data_folder}/retrorules.db')
    df.to_sql('RETRORULES', conn, if_exists='replace')

def create_indexes():
    print('Creating indexes in retrorules database')
    conn = sqlite3.connect(f'{data_folder}/retrorules.db')
    c = conn.cursor()
    col_name = 'Substrate_ID'
    table = 'RETRORULES'
    c.execute(f"CREATE INDEX index_{col_name} ON {table}({col_name})")
    conn.commit()

    conn = sqlite3.connect(f'{data_folder}/retrorules.db')
    c = conn.cursor()
    col_name = 'Rule_usage'
    table = 'RETRORULES'
    c.execute(f"CREATE INDEX index_{col_name} ON {table}({col_name})")
    conn.commit()

    conn = sqlite3.connect(f'{data_folder}/retrorules.db')
    c = conn.cursor()
    col_name = 'Diameter'
    table = 'RETRORULES'
    c.execute(f"CREATE INDEX index_{col_name} ON {table}({col_name})")
    conn.commit()

def create_smi_file():
    print('Creating retrorules smi file')
    fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator()
    df = pd.read_csv(f'{data_folder}/retrorules_rr02_flat_all.tsv', sep='\t')
    smi_only = df[['Substrate_ID', 'Substrate_SMILES']]
    smi_only = smi_only.drop_duplicates()

    list_mols = []
    for i, smi in enumerate(list(smi_only.Substrate_SMILES)):
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            print(f'{smi} - {i} of {len(list(smi_only.Substrate_SMILES))} - error loading')
            mol = np.nan
        list_mols.append(mol)

    smi_only['mols'] = list_mols
    smi_only.dropna(inplace=True)
    smi_only['fp'] = smi_only['mols'].apply(fp_gen.GetFingerprint)
    smi_only.drop(columns=['mols'], inplace=True)

    smi_only.to_hdf(f'{data_folder}/rr_smi.h5', key='df', mode='w')

def run():
    make_rr_dir()
    download_retrorules_files()
    delete_old_db()
    create_db()
    data_to_db()
    create_indexes()
    create_smi_file()
    delete_downloaded_files()
    print('Retrorules db loaded')

if __name__ == '__main__':
    run()
