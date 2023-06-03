import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from scipy.spatial.distance import squareform
from pathlib import Path

from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_distances.route_distance_model import \
    RouteDistanceModel
from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_distances.utils import \
    preprocess_reaction_tree, collate_trees


class Route_Distancer():
    model = None
    model_path = str(Path(__file__).parents[3]) + '/retrosynthesis_engine/data/aizynthfinder/chembl_10k_route_distance_model.ckpt'


    def __init__(self):
        self.smi_fp_dict = {}
        self.load_model()

    @classmethod
    def load_model(cls):
        if cls.model is None:
            cls.model = RouteDistanceModel.load_from_checkpoint(cls.model_path)
            cls.model.eval()

    def get_distances(self, pathways):
        if len(pathways) <= 2:
            return None

        trees = self.get_trees(pathways)
        processed_trees = [preprocess_reaction_tree(tree) for tree in trees]
        tree_data = collate_trees(processed_trees)
        pred_torch = self.model(tree_data)
        pred_np = pred_torch.detach().numpy()
        return squareform(pred_np)

    def get_trees(self, pathways):
        trees = [pathway.aizynthfinder_tree() for pathway in pathways]
        [self._add_fingerprints(tree) for tree in trees]
        return trees

    def _add_fingerprints(self, tree, radius=2, nbits=2048):
        smi = tree["smiles"]
        if smi not in self.smi_fp_dict:
            mol = Chem.MolFromSmiles(tree["smiles"])
            rd_fp = AllChem.GetMorganFingerprintAsBitVect(
                mol, radius=radius, nBits=nbits, useFeatures=False, useChirality=True
            )
            np_fp = np.empty(radius, np.int8)
            DataStructs.ConvertToNumpyArray(rd_fp, np_fp)
            self.smi_fp_dict[smi] = np_fp

        tree["fingerprint"] = self.smi_fp_dict[smi]

        for child in tree.get("children", []):
            self._add_fingerprints(child, radius, nbits)


if __name__ == '__main__':
    model_path = str(Path(__file__).parents[3]) + '/retrosynthesis_engine/data/aizynthfinder' \
                                                  '/chembl_10k_route_distance_model.ckpt '
    print(model_path)
    print("/Users/willfinnigan/PycharmProjects/retrobiocat/retrobiocat_web/retro/retrosynthesis_engine/data/aizynthfinder/chembl_10k_route_distance_model.ckpt")

    distancer = Route_Distancer()