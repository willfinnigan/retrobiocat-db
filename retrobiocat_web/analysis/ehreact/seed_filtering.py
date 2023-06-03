import copy

from rdkit import Chem
from retrobiocat_web.analysis.ehreact.helpers.rdkit import read_in_smiles_unique, preprocess_seeds
from retrobiocat_web.analysis.ehreact.tests.smiles_lists import imi_red_smiles, imi_red_seeds, imi_red_prod_smiles, \
    imi_red_prod_seeds, orc_smiles, orc_seeds, lipase_smiles, lipase_seeds


def filter_smiles_with_smarts(smiles, list_pattern_smis):
    if len(list_pattern_smis) == 0:
        return smiles

    all_matches = []
    for pattern_smi in list_pattern_smis:
        smi_mol_dict = {s: Chem.MolFromSmiles(s) for s in smiles}
        patt = Chem.MolFromSmiles(pattern_smi)
        matches = {s: m for s, m in smi_mol_dict.items() if m.HasSubstructMatch(patt)}.keys()
        all_matches += matches

    all_matches = list(set(all_matches))

    non_matches = []
    if len(all_matches) != smiles:
        for smi in smiles:
            if smi not in all_matches:
                non_matches.append(smi)

    return all_matches

def get_non_matches(smiles, seeds):
    if len(smiles) == 0:
        return [], [], {}

    smiles, smiles_dict, skipped = read_in_smiles_unique(smiles)
    seeds, rule_dict, num_smiles_seed = preprocess_seeds(seeds, smiles, smiles_dict)

    non_matched_smiles = []
    all_matches = []
    if num_smiles_seed != len(smiles):
        for seed, seed_dict in rule_dict.items():
            all_matches += seed_dict['smiles']

        for smi in smiles:
            if smi not in all_matches:
                non_matched_smiles.append(smi)

    ori_smiles = {smi: smiles_dict[smi]['original_smi'] for smi in smiles_dict}

    return non_matched_smiles, all_matches, ori_smiles

def generate_non_match_seed(smiles, seeds):
    non_matches, matches, ori_smiles = get_non_matches(smiles, seeds)
    if len(non_matches) == 0:
        return []

    if len(non_matches) == 1:
        return [non_matches[0]]

    smiles, smiles_dict, skipped = read_in_smiles_unique(non_matches)
    non_match_seeds, rule_dict, num_smiles_seed = preprocess_seeds([], smiles, smiles_dict)
    return non_match_seeds



def are_any_seed_substructures_of_each_other(seeds):
    if len(seeds) <= 1:
        return seeds

    redundent_seeds = []
    for i, seed in enumerate(seeds):
        other_seeds = copy.deepcopy(seeds)
        other_seeds.remove(seed)
        matches = filter_smiles_with_smarts(other_seeds, [seed])
        redundent_seeds += matches

    for seed in redundent_seeds:
        seeds.remove(seed)

    return seeds

if __name__ == '__main__':
    #test_seeds = ['C(=O)O[H]', '[H]COC=O']
    #new_seeds = are_any_seed_substructures_of_each_other(test_seeds)
    #print(new_seeds)

    smis = lipase_smiles
    seeds = lipase_seeds

    all_matches = filter_smiles_with_smarts(smis, seeds)
    print(f"Initial filter - {len(smis)} to start with, {len(all_matches)} after filtering")

    non_matches, matches = get_non_matches(smis, seeds)
    print(f"{len(non_matches)} non_matches using ehreact")
    print(non_matches)

    non_match_seeds = generate_non_match_seed(smis, seeds)
    print(non_match_seeds)

