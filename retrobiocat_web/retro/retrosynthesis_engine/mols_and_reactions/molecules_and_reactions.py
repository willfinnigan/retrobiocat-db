import uuid

from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.complexity_evaluator import \
    Complexity_Evaluator
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.starting_material_evaluator import \
    StartingMaterialEvaluator

class MoleculeFactory():

    def __init__(self, target_mol=None, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.logger = add_logger('MoleculeFactory', level=log_level)

        self.target_mol = target_mol

        self.complexity_evaluator = Complexity_Evaluator()
        self.starting_material_evaluator = StartingMaterialEvaluator(config=self.config)

    def make_mol(self, smi, force_rdkit=False):
        self.logger.debug(f"Making molecule for smiles: {smi}")
        if force_rdkit:
            smi = rdkit_smile(smi)

        if smi is None:
            self.logger.error("smi is None when trying to make molecule")
            Exception('smi is None')
            return None

        complexity = self.complexity_evaluator.get_complexity(smi)
        in_stock, in_stock_info = self.starting_material_evaluator.eval(smi,
                                                                        mode=self.config.source_mol_mode,
                                                                        vendors=self.config.source_mol_vendors)

        relative_complexity = 0
        if self.target_mol is not None:
            relative_complexity = complexity - self.target_mol.complexity

        self.logger.debug(f"Created mol - {smi}, complexity={complexity}, relative_complexity={relative_complexity}, in_stock={in_stock}")

        mol = Molecule(smi, complexity,
                       relative_complexity=relative_complexity,
                       in_stock=in_stock,
                       in_stock_info=in_stock_info)

        if mol is None:
            self.logger.error(f"Created molecule is None.  Smi should be {smi}")
            Exception('Mol is None')

        return mol

    def create_target_mol(self, smi, force_rdkit=False):
        self.target_mol = self.make_mol(smi, force_rdkit=force_rdkit)

        if self.config.target_always_not_buyable:
            self.target_mol.in_stock = 0
            self.target_mol.in_stock_info = {}

        return self.target_mol

class Molecule():

    def __init__(self, smi, complexity, relative_complexity=None, in_stock=None, in_stock_info=None):
        self.smi = smi
        self.stats = {'complexity': complexity,
                      'relative_complexity': relative_complexity,
                      'in_stock': in_stock,
                      'in_stock_info': in_stock_info}

        self.complexity = complexity
        self.relative_complexity = relative_complexity
        self.in_stock = in_stock
        self.in_stock_info = in_stock_info


    def __repr__(self):
        return f"{self.smi}"

class Reaction(object):

    def __init__(self, name, reactant, substrates, metadata=None, unique_id=None, rxn_type=''):
        self.name = name
        self.unique_id = unique_id
        if self.unique_id is None:
            self.unique_id = uuid.uuid4()
        self.reactant = reactant
        self.substrates = substrates
        self.complexity_change = self._get_complexity_change()
        self.metadata = metadata
        if self.metadata is None:
            self.metadata = {}
        self.rxn_type = rxn_type
        self.graph_node_name = None

    def __repr__(self):
        return f"{self.name} {self.unique_id}"

    def _get_complexity_change(self):
        substrate_complexities = [substrate.complexity for substrate in self.substrates]
        max_substrate_complexity = max(substrate_complexities)
        complexity_change = max_substrate_complexity - self.reactant.complexity
        return complexity_change

    def reaction_smi(self):
        reaction_smi = ""
        for substrate in self.substrates:
            reaction_smi += substrate.smi
            reaction_smi += '.'

        if reaction_smi[-1] == '.':
            reaction_smi = reaction_smi[0:-1]

        reaction_smi += '>>'
        reaction_smi += self.reactant.smi

        return reaction_smi