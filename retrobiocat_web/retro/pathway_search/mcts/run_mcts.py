from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.pathway_search.mcts.config import MCTS_Config
from retrobiocat_web.retro.pathway_search.mcts.mcts import MCTS


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    """
    target_smiles = 'CC(=CCCC(C)(C=C)O)C'
    network = Network(target_smiles=target_smiles)

    retro_config = RetrosynthesisConfig()
    mcts_config = MCTS_Config()

    mcts_config.max_search_time = 60
    mcts_config.maxLength = 5
    mcts_config.biosynthesis_expansion = True
    mcts_config.chemistry_expansion = False
    mcts_config.biocatalysis_expansion = False

    mcts = MCTS(network, config=mcts_config, retro_config=retro_config, log_level='DEBUG')

    mcts.run()
    """

    target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
    network = Network(target_smiles=target_smiles)

    retro_config = RetrosynthesisConfig()
    mcts_config = MCTS_Config()

    mcts_config.max_search_time = 120
    mcts_config.maxLength = 2
    mcts_config.biosynthesis_expansion = False
    mcts_config.chemistry_expansion = True
    mcts_config.biocatalysis_expansion = False

    mcts = MCTS(network, config=mcts_config, retro_config=retro_config, log_level='DEBUG')

    mcts.run()


