from retrobiocat_web import logging
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig


class SmallMolRemover():

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if config is None:
            self.config = RetrosynthesisConfig()

        self.logger = logging.add_logger('SmallMolRemover', level=log_level)

    def remove(self, products):
        removed = []
        kept = []
        for mol in products:
            if mol.smi in self.config.small_precursors:
                removed.append(mol)
            else:
                kept.append(mol)

        if len(removed) != 0:
            self.logger.debug(f"Removed {len(removed)} small mols from {products} - removed = {removed}, kept = {kept}")

        return kept, removed


if __name__ == '__main__':
    remover = SmallMolRemover(log_level='DEBUG')
    products = ['CCCC', 'N']
    kept, removed = remover.remove(products)