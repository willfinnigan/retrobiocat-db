from pathlib import Path

from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.sqlite_source_mol.connect_sqlitedb import \
    SQLite_Database
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.sqlite_source_mol.query_sqlitedb import \
    DB_Query_SQLite


data_folder = str(Path(__file__).parents[3]) + '/data/buyability'

class StartingMaterialEvaluator():
    vendor_urls = {'mcule': 'https://mcule.com/[[ID]]',
                   'sigma': 'https://www.sigmaaldrich.com/GB/en/search/[[ID]]?focus=products&page=1&perpage=30&sort=relevance&term=[[ID]]&type=product',
                   'lifechem': 'https://shop.lifechemicals.com/compound/[[ID]]',
                   'apollo': 'https://store.apolloscientific.co.uk/search?search=[[ID]]',
                   'alfa': 'https://www.alfa.com/en/catalog/[[ID]]',
                   'zinc': 'https://zinc.docking.org/substances/[[ID]]',
                   'flurochem': 'http://www.fluorochem.co.uk/Products/Product?code=[[ID]]',
                   'ecmdb': 'https://ecmdb.ca/compounds/[[ID]]'}

    def __init__(self, config=RetrosynthesisConfig()):

        self.config = config
        db_path = data_folder + '/source_mols.db'
        self.database = SQLite_Database(db_path)
        self.query = DB_Query_SQLite(self.database)
        self.cache_column_names = {}
        self.cache_vendor_names = {}

    def eval(self, smi, mode='building_blocks', vendors=None):
        if mode not in ['building_blocks', 'metabolites']:
            Exception(f'Mode specified ({mode}) is invalid')
            return 0, {}

        if self.is_mol_chiral_and_this_is_not_allowed(smi):
            return 0, {}

        if vendors is not None:
            vendors = [v for v in vendors if v in self.vendor_names(mode)]  # only pass vendors which are actually in db

        result = self.query.smiles_lookup(smi, mode, vendors=vendors)
        if result is not None:
            info = self._process_info(result, mode)
            return 1, info
        return 0, {}

    def is_mol_chiral_and_this_is_not_allowed(self, smi):
        if self.config.source_mols_can_be_chiral:
            return False
        elif '@' in smi:
            return True
        else:
            return False

    def column_names(self, mode='building_blocks'):
        if mode not in self.cache_column_names:
            self.cache_column_names[mode] = self.query.get_column_names(mode)
        return self.cache_column_names[mode]

    def vendor_names(self, mode='building_blocks'):
        if mode not in self.cache_vendor_names:
            columns = self.column_names(mode)
            vendors = []
            for col in columns:
                if '_id' in col:
                    vendors.append(col.replace('_id', ''))
            self.cache_vendor_names[mode] = vendors

        return self.cache_vendor_names[mode]


    def _process_info(self, result, mode):
        columns = self.column_names(mode)
        vendors = self.vendor_names(mode)
        info = {k: v for k, v in zip(columns, result)}
        info.pop('id')
        info = {k: v for k, v in info.items() if v is not None}

        vendor_info = {}
        for col, value in info.items():
            for vendor in vendors:
                if vendor in col:
                    if vendor not in vendor_info:
                        vendor_info[vendor] = {}
                    vendor_info[vendor][col.replace(f"{vendor}_", '')] = value
                    if ('_id' in col) and (vendor in self.vendor_urls):
                        url = self.vendor_urls[vendor].replace('[[ID]]', value)
                        vendor_info[vendor]['url'] = url
        return vendor_info

if __name__ == '__main__':
    sme = StartingMaterialEvaluator()
    available, info = sme.eval('CCC(=O)C(=O)O', mode='metabolites')
    columns = sme.column_names()
    print(available)
