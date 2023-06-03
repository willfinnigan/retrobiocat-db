from retrobiocat_web import logging
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.molecules_and_reactions import Molecule, Reaction, MoleculeFactory


class ReactionCreator():
    """
    Take the output of rule applicator and process for addition to network
    """

    def __init__(self, network, config=None, log_level='WARNING'):
        self.network = network

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()
        self.logger = logging.add_logger('ReactionCreator', log_level)

        self.mol_factory = MoleculeFactory(config=config, log_level=log_level)

        if network.target_smiles is not None:
            self.mol_factory.target_mol = self._find_or_create_mol(network.target_smiles)
        else:
            self.logger.debug("Could not make target mol as network target smiles is None")

    def create_new_reactions(self, target_smi, products_dict, metadata=None, rxn_type=''):
        # product dict will = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1'], ['C1=N[C@H](c2ccccc2)CCC1']]}

        if metadata is None:
            metadata = {}

        self.logger.debug(f"Find or create reactant mol - {target_smi}")
        target_mol = self._find_or_create_mol(target_smi)
        new_reactions = []
        for rxn_name, products_list in products_dict.items():
            rxn_metadata = metadata.get(rxn_name, None)
            mols = self._get_mols(products_list)
            reactions = self._create_reactions(target_mol, mols, rxn_name, rxn_type, rxn_metadata)
            new_reactions.extend(reactions)

        self.logger.info(f"Created {len(new_reactions)} reactions")
        for reaction in new_reactions:
            self.logger.debug(f"{reaction.name} - {reaction.reactant}-->{reaction.substrates}")
        return new_reactions

    def _create_reactions(self, target_mol, mols, rxn_name, rxn_type, metadata):
        reactions = []
        for set_of_mols in mols:
            if len(set_of_mols) != 0:
                new_reaction = Reaction(rxn_name, target_mol, set_of_mols, rxn_type=rxn_type, metadata=metadata)
                reactions.append(new_reaction)
        return reactions

    def _get_mols(self, products_list):
        self.logger.debug(f"Getting mols for product list {products_list}")
        all_mols = []
        for list_smis in products_list:
            if None in list_smis:
                self.logger.error(f"None in products list")
                self.logger.error(f"{products_list}")
            mols = []
            for smi in list_smis:
                mols.append(self._find_or_create_mol(smi))
            all_mols.append(mols)
        return all_mols

    def _find_or_create_mol(self, smi):
        if smi is None:
            self.logger.error("Smi is None when trying to find or create molecule")
            Exception("Smi is None when trying to find or create molecule")
            return None

        if smi in self.network.substrate_nodes:
            attributes = self.network.graph.nodes[smi].get('attributes', {})
            complexity = attributes.get('complexity', None)
            relative_complexity = attributes.get('relative_complexity', 0)
            in_stock = attributes.get('is_starting_material', 0)
            in_stock_info = attributes.get('starting_material_info', {})
            if complexity is not None:
                return Molecule(smi, complexity,
                                relative_complexity=relative_complexity,
                                in_stock=in_stock,
                                in_stock_info=in_stock_info)

        return self.mol_factory.make_mol(smi)

