from pathlib import Path
import pandas as pd
from rdkit import Chem
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig


data_folder = str(Path(__file__).parents[4]) + '/data/retrorules'

class RetroRules_Processor():
    cofactors_df = None

    def __init__(self, config=None, log_level='WARNING'):
        self.logger = add_logger('RetroBioCat_Expander', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

    @classmethod
    def load_cofactor_df(cls):
        if cls.cofactors_df is None:
            cls.cofactors_df = pd.read_csv(f"{data_folder}/retrorules_cofactors.tsv",
                                            sep='\t',
                                            names=['inchi', 'name', 'reactions'])

    def process(self, products_dict):
        if self.config.rr_remove_cofactors:
            products_dict, cofactors = self.remove_cofactors(products_dict)

        return products_dict

    def remove_cofactors(self, products_dict):
        self.load_cofactor_df()
        cofactors = {}
        for id in products_dict:
            for i, products in enumerate(products_dict[id]):
                for smi in products:
                    mol = Chem.MolFromSmiles(smi)
                    inchi = Chem.inchi.MolToInchi(mol)
                    if inchi in list(self.cofactors_df['inchi']):
                        products_dict[id][i].remove(smi)
                        cofactors[id] = self.cofactors_df['name']
        return products_dict, cofactors