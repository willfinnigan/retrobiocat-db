import os

import requests
from pathlib import Path
import zipfile

data_folder = str(Path(__file__).parents[4]) + '/data'
analysis_folder = str(Path(__file__).parents[4]) + '/analysis/analysis_data'

def download_file(url, filename):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(filename, "wb") as fileobj:
            for chunk in response.iter_content(chunk_size=1024):
                fileobj.write(chunk)

def create_other_folders():
    ssn = f"{analysis_folder}/ssn"
    blasts = f"{analysis_folder}/all_by_all_blast"
    Path(ssn).mkdir(parents=True, exist_ok=True)
    Path(blasts).mkdir(parents=True, exist_ok=True)

def download_building_blocks():
    print('DOWNLOADING BUILDING BLOCKS')
    Path(f'{data_folder}/buyability').mkdir(parents=True, exist_ok=True)
    bb_path = data_folder + '/buyability/source_mols.db'
    print(bb_path)
    download_file("https://figshare.com/ndownloader/files/37376029", bb_path)

def download_aizynthfinder_files():
    print('DOWNLOADING AISYNTHFINDER FILES')
    Path(f'{data_folder}/aizynthfinder').mkdir(parents=True, exist_ok=True)
    model_path = data_folder + '/aizynthfinder/uspto_model.hdf5'
    filter_path = data_folder + '/aizynthfinder/filter_policy_all.hdf5'
    templates_path = data_folder + '/aizynthfinder/uspto_templates.hdf5'
    download_file("https://figshare.com/ndownloader/files/23086454", model_path)
    download_file("https://figshare.com/ndownloader/files/23086457", templates_path)
    download_file("https://figshare.com/ndownloader/files/25584743", filter_path)



def download_ringbreaker_files():
    print('DOWNLOADING RINGBREAKER FILES')
    folder = f'{data_folder}/ringbreaker'
    Path(folder).mkdir(parents=True, exist_ok=True)

    download_file("https://figshare.com/ndownloader/files/22789826", f"{folder}/model.zip")
    with zipfile.ZipFile(f"{folder}/model.zip", 'r') as zip_ref:
        zip_ref.extractall(folder)

    download_file("https://figshare.com/ndownloader/files/22789832", f"{folder}/data.zip")
    with zipfile.ZipFile(f"{folder}/data.zip", 'r') as zip_ref:
        zip_ref.extractall(folder)

    os.remove(f"{folder}/model.zip")
    os.remove(f"{folder}/data.zip")

if __name__ == '__main__':
    #download_aizynthfinder_files()
    #download_building_blocks()
    download_ringbreaker_files()
