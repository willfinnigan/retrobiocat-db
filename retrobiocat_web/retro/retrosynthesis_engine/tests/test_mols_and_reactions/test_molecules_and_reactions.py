from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.molecules_and_reactions import Molecule, \
    Reaction, MoleculeFactory


class Test_Molecules():

    def test_can_create_a_molecule(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        mol = mol_factory.make_mol('CCC=O')
        assert mol.smi == 'CCC=O'

    def test_created_mol_has_a_complexity(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        mol = mol_factory.make_mol('CCC=O')
        assert mol.complexity > 0

    def test_created_mol_has_in_stock(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        mol = mol_factory.make_mol('CCC=O')
        assert mol.in_stock == 1
        assert 'askcos' in mol.in_stock_info

    def test_created_mol_has_not_in_stock_if_not_for_specified_vendor(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        mol_factory.config.source_mol_vendors = ['askcos']
        mol = mol_factory.make_mol("O=C1CCC[C@@H](c2ccccc2)N1")
        assert mol.in_stock == 0

    def test_created_mol_in_stock_if_all_vendors_specified(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        mol_factory.config.source_mol_vendors = None
        mol = mol_factory.make_mol("O=C1CCC[C@@H](c2ccccc2)N1")
        assert mol.in_stock == 1


class Test_Reactions():

    def test_can_make_a_reaction(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        reactant = mol_factory.make_mol('CCC=O')
        substrates = [mol_factory.make_mol('CCCO')]
        reaction = Reaction('test_rxn', reactant, substrates, metadata=None, unique_id=None, rxn_type='')
        assert reaction.name == 'test_rxn'

    def test_reaction_has_a_complexity_change(self):
        mol_factory = MoleculeFactory(log_level='DEBUG')
        reactant = mol_factory.make_mol('N[C@@H](CC(N1CCN2C(C1)=NN=C2C(F)(F)F)=O)CC3=C(F)C=C(F)C(F)=C3') # complex
        substrates = [mol_factory.make_mol('CCCO')]
        reaction = Reaction('test_rxn', reactant, substrates, metadata=None, unique_id=None, rxn_type='')
        print(f"\n Complexity change = {reaction.complexity_change}")
        assert reaction.complexity_change < 0  # complexity is complex product from simple materials, so should be -ve




