from rdkit.Chem import AllChem, Draw, rdChemReactions
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from requests.utils import quote

from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs

import time


"""def moltosvg(mol,molSize=(300,300)):
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    opts = drawer.drawOptions()
    opts.padding = 0.15
    drawer.SetFontSize(1.0)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg"""

def moltosvg(mol,molSize=(300,300),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    opts = drawer.drawOptions()
    opts.padding = 0.15
    drawer.SetFontSize(1.0)
    opts.addStereoAnnotation = True
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg

def smile_to_svg_url(smiles, size=(300,300), kekulize=True):
    mol = Chem.MolFromSmiles(smiles)
    svg = moltosvg(mol, molSize=size, kekulize=kekulize)
    url = "data:image/svg+xml;charset=utf-8," + quote(svg)
    return url

def morgan_fingerprint_vis(product_smiles, sim_smiles):
    mol = Chem.MolFromSmiles(sim_smiles)
    refmol = Chem.MolFromSmiles(product_smiles)

    fp1 = SimilarityMaps.GetMorganFingerprint(mol)
    fp2 = SimilarityMaps.GetMorganFingerprint(refmol)

    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, SimilarityMaps.GetMorganFingerprint)
    return DataStructs.FingerprintSimilarity(fp1, fp2)

def smiles_rxn_to_svg(smiles_rxn, rxnSize=(900,300)):
    try:
        rxn = rdChemReactions.ReactionFromSmarts(smiles_rxn, useSmiles=True)
        d = Draw.MolDraw2DSVG(rxnSize[0], rxnSize[1])
        opts = d.drawOptions()
        opts.addStereoAnnotation = True
        d.DrawReaction(rxn)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        return svg
    except:
        return None

