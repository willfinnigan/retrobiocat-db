
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem, DataStructs

from requests.utils import quote
import time
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram

from retrobiocat_web.mongo.model_queries import sequence_queries
from retrobiocat_web.mongo.model_queries.get_activity_data import get_activity_data
from retrobiocat_web.similarity.similarity_tools import get_single_fp


def activity_df(activity_data):
    return pd.DataFrame(list(activity_data))


def get_unique_enzymes(activity_data, enzyme_col):
    enzymes = []
    for data in activity_data:
        if enzyme_col not in data:
            if 'nan' not in enzymes:
                enzymes.append('nan')
        elif data[enzyme_col] not in enzymes:
            enzymes.append(data[enzyme_col])

    return enzymes


def get_unique_smiles(activity_data, smi_col):
    unique_smis = []
    for data in activity_data:
        if smi_col in data:
            smi = data[smi_col]
            if smi not in unique_smis:
                unique_smis.append(smi)

    if "" in unique_smis:
        unique_smis.remove("")

    return unique_smis


def get_substrates_fps(unique_smiles):
    fps = []
    for smi in unique_smiles:
        fp = get_single_fp(smi)
        fps.append(fp)

    return fps


def make_fp_dict(substrates, fps):
    fp_dict = {}
    for i, substrate in enumerate(substrates):
        fp_dict[substrate] = fps[i]
    return fp_dict


def get_enzyme_embeddings(enzymes):
    enzyme_embeddings_map = {}
    seqs = sequence_queries.seq_objs_from_seq_names(enzymes)
    # seqs = Sequence.objects(db.Q(enzyme_name__in=enzymes) & db.Q(unirep__ne=None))
    embeddings = sequence_queries.get_embeddings(seqs)

    for seq, embedding in zip(seqs, embeddings):
        name = seq.enzyme_name
        enzyme_embeddings_map[name] = embedding

    return enzyme_embeddings_map


def get_mol_svgs(substrates, molSize=(250, 100)):
    def create_svg(smi, molSize=molSize):
        mol = Chem.MolFromSmiles(smi)
        d2d = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
        # d2d.SetDrawOptions()
        opts = d2d.drawOptions()
        opts.clearBackground = False
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        svg_text = d2d.GetDrawingText()
        return svg_text

    def get_svg_url(svg):
        return "data:image/svg+xml;charset=utf-8," + quote(svg)

    svg_map = {}
    svg_url_map = {}

    for smi in substrates:
        svg = create_svg(smi)
        svg_map[smi] = svg
        svg_url_map[smi] = get_svg_url(svg)

    return svg_map, svg_url_map


def order_smis(smis, threshold=0.7):
    size = len(smis)
    hmap = np.empty(shape=(size, size))
    fps = get_substrates_fps(smis)
    for i, fp_i in enumerate(fps):
        similarities = DataStructs.BulkTanimotoSimilarity(fp_i, list(fps))
        for j, sim in enumerate(similarities):
            hmap[i, j] = sim

    linked = linkage(hmap, 'single')
    results = dendrogram(linked, no_plot=True, distance_sort='ascending', color_threshold=threshold)

    colours = results['color_list']
    sub_numbered = list(map(int, results['ivl']))

    substrates_ordered = []
    for i in sub_numbered:
        substrates_ordered.append(smis[i])

    return substrates_ordered, colours


def calc_mol_stats(df, smi, mol_type):
    smi_df = df[(df[mol_type] == smi)]
    active_df = smi_df[(smi_df['binary'] == 'True')]
    inactive_df = smi_df[(smi_df['binary'] == 'False')]

    stats = {'active': len(active_df.index),
             'inactive': len(inactive_df.index),
             'active_enzymes': len(active_df.enzyme_name.unique()),
             'active_papers': len(active_df.paper.unique()),
             'inactive_enzymes': len(inactive_df.enzyme_name.unique()),
             'inactive_papers': len(inactive_df.paper.unique())
             }

    return stats


class DataQuery(object):

    def __init__(self, enzyme_type=None, reaction=None,
                 paper_id=None, enzyme_names=None,
                 only_reviewed=True, smi_col='product_1_smiles',
                 smiles=None, remove_negative=None,
                 log_level=0):

        self.activity_data = None
        self.activity_df = None
        self.mol_stats = None
        self.enzyme_type = enzyme_type
        if self.enzyme_type == 'All':
            self.enzyme_type = None
        self.reaction = reaction
        if self.reaction == 'All':
            self.reaction = None
        self.paper_id = paper_id
        self.enzyme_names = enzyme_names
        self.only_reviewed = only_reviewed
        self.smi_col = smi_col
        self.smiles = smiles
        self.convert_safe_hash_back()
        self.log_level = log_level
        self.remove_negative = remove_negative

        self.get_data()

    def convert_safe_hash_back(self):
        if self.smiles != None:
            if type(self.smiles) == str:
                self.smiles = self.smiles.replace("[hash]", "#")

    def unique_smiles(self):
        return get_unique_smiles(self.activity_data, self.smi_col)

    def unique_enzymes(self):
        return get_unique_enzymes(self.activity_data, 'enzyme_name')

    def reaction_names(self):
        df = self.get_activity_df()

        if 'reaction' not in df.columns:
            return None

        return list(df['reaction'].unique())

    def enzyme_counts(self, only_top=100):
        df = self.get_activity_df()
        enzyme_count_total = df['enzyme_name'].value_counts()
        enzyme_to_keep = list(enzyme_count_total.index)
        if len(enzyme_to_keep) > only_top:
            enzyme_to_keep = enzyme_to_keep[0:only_top]

        filtered_df = df[df['enzyme_name'].isin(enzyme_to_keep)]
        mol_count_active = filtered_df[filtered_df['binary'] == 'True']['enzyme_name'].value_counts().to_dict()
        mol_count_high = filtered_df[filtered_df['categorical'] == 'High']['enzyme_name'].value_counts().to_dict()
        mol_count_medium = filtered_df[filtered_df['categorical'] == 'Medium']['enzyme_name'].value_counts().to_dict()
        mol_count_low = filtered_df[filtered_df['categorical'] == 'Low']['enzyme_name'].value_counts().to_dict()
        mol_count_inactive = filtered_df[filtered_df['binary'] == 'False']['enzyme_name'].value_counts().to_dict()

        value_counts = {}
        for enz in enzyme_to_keep:
            high = mol_count_high.get(enz, 0)
            medium = mol_count_medium.get(enz, 0)
            low = mol_count_low.get(enz, 0)
            active = mol_count_active.get(enz, 0) - (high + medium + low)
            inactive = mol_count_inactive.get(enz, 0)
            enzyme_type = df.loc[df['enzyme_name'] == enz, 'enzyme_type'].iloc[0]
            # enzyme_type = df[df['enzyme_name'] == enz]['enzyme_type'].values[0]
            # enzyme_type = df['enzyme_type'][df['enzyme_name'].loc(enz)]
            value_counts[enz] = {'high': high,
                                 'medium': medium,
                                 'low': low,
                                 'active': active,
                                 'inactive': inactive,
                                 'total': high + medium + low + active + inactive,
                                 'enzyme_type': enzyme_type}
        return value_counts

    def mol_counts(self, only_top=np.inf):
        df = self.get_activity_df()
        mol_count_total = df[self.smi_col].value_counts()

        mols_to_keep = list(mol_count_total.index)
        if len(mols_to_keep) > only_top:
            mols_to_keep = mols_to_keep[0:only_top]

        filtered_df = df[df[self.smi_col].isin(mols_to_keep)]
        mol_count_active = filtered_df[filtered_df['binary'] == 'True'][self.smi_col].value_counts().to_dict()
        mol_count_high = filtered_df[filtered_df['categorical'] == 'High'][self.smi_col].value_counts().to_dict()
        mol_count_medium = filtered_df[filtered_df['categorical'] == 'Medium'][self.smi_col].value_counts().to_dict()
        mol_count_low = filtered_df[filtered_df['categorical'] == 'Low'][self.smi_col].value_counts().to_dict()
        mol_count_inactive = filtered_df[filtered_df['binary'] == 'False'][self.smi_col].value_counts().to_dict()

        value_counts = {}
        for smi in mols_to_keep:
            high = mol_count_high.get(smi, 0)
            medium = mol_count_medium.get(smi, 0)
            low = mol_count_low.get(smi, 0)
            active = mol_count_active.get(smi, 0) - (high + medium + low)
            inactive = mol_count_inactive.get(smi, 0)
            value_counts[smi] = {'high': high,
                                 'medium': medium,
                                 'low': low,
                                 'active': active,
                                 'inactive': inactive,
                                 'total': high + medium + low + active + inactive}
        return value_counts

    def get_mol_stats(self):
        if self.mol_stats is None:
            df = self.get_activity_df()
            self.mol_stats = {}
            for smi in self.unique_smiles():
                self.mol_stats[smi] = calc_mol_stats(df, smi, self.smi_col)
        return self.mol_stats

    def get_data(self):
        t0 = time.time()
        if self.remove_negative == 'all':
            include_negative = False
        else:
            include_negative = True

        self.activity_data = get_activity_data(enzyme_type=self.enzyme_type,
                                               reaction=self.reaction,
                                               paper_id=self.paper_id,
                                               enzyme_names=self.enzyme_names,
                                               smiles=self.smiles,
                                               smi_col=self.smi_col,
                                               only_reviewed=self.only_reviewed,
                                               include_negative=include_negative,
                                               as_pymongo=True)

        # removes smiles for which there is only negative data
        if self.remove_negative == 'only_negative':
            self.filter_out_only_negative_data()

        t1 = time.time()
        self.log(f"Queried database in {round(t1 - t0, 2)} seconds")
        self.log(f"{len(self.activity_data)} entries")

    def filter_out_only_negative_data(self):

        unique_smiles = get_unique_smiles(self.activity_data, self.smi_col)
        activity_df = self.activity_df()
        smis_with_no_activity = []
        for smi in unique_smiles:
            filtered_df = activity_df[(activity_df[self.smi_col] == smi)]
            if len(filtered_df[filtered_df['binary'] == 'True']) == 0:
                smis_with_no_activity.append(smi)

        count_filter = 0
        for i, data in enumerate(self.activity_data):
            if data[self.smi_col] in smis_with_no_activity:
                self.activity_data.pop(i)
                count_filter += 1

        self.log(f"Filtered {count_filter} only negative data entries")

    def get_activity_df(self):
        if self.activity_df is None:
            self.activity_df = pd.DataFrame(list(self.activity_data))
            if 'categorical' not in self.activity_df.columns:
                self.activity_df['categorical'] = np.nan
            if 'binary' not in self.activity_df.columns:
                self.activity_df['binary'] = np.nan
        return self.activity_df

    def reactions(self):
        reactions = []
        for data in self.activity_data:
            r = data.get('reaction', None)
            if r is not None and r not in reactions:
                reactions.append(r)
        return reactions

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"DataQuery ({level}): {msg}")


if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection

    make_default_connection()

    dq = DataQuery(enzyme_type='IRED', log_level=1, remove_negative='all', smi_col='substrate_1_smiles')

    stats = dq.get_mol_stats()
    print(stats)