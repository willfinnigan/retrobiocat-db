from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from requests.utils import quote
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.models.mappers import CategoricalColorMapper
from bokeh.palettes import Category20

def create_svg(mol, molSize=(250, 100)):
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

def get_mol_svgs(smis, molSize=(250, 100)):
    svgs = []
    svg_urls = {}

    for smi in smis:
        mol = Chem.MolFromSmiles(smi)
        svg = create_svg(mol, molSize=molSize)
        svgs.append(svg)
        svg_urls[smi] = get_svg_url(svg)

    return svgs, svg_urls
