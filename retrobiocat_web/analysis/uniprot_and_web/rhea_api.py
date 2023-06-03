import requests
import bioservices
import pandas as pd
import io
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
from retrobiocat_web.analysis.drawing import images
from rdkit import Chem
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, UniRef50
from retrobiocat_web.mongo.models.other_models import Chebi, Rhea

def get_rhea(rhea_id):
    rhea = Rhea.objects(code=rhea_id).first()

    if rhea is None:
        url = "https://www.rhea-db.org/rhea?"

        parameter = {
            "query": rhea_id,
            "columns": "equation,chebi-id,ec",
            "format": 'tsv',
            "limit": 10,
        }

        response = requests.get(url, params=parameter)
        df = pd.read_csv(io.StringIO(response.text), sep='\t')
        chebis = df.iloc[0]['ChEBI identifier'].split(';')
        equation = df.iloc[0]['Equation']
        ec = df.iloc[0]['EC number']
        if str(ec) == 'nan':
            ec = []
        else:
            ec = ec.split(';')

        ec = ec_to_ec3(ec)

        rhea = Rhea(code=rhea_id,
                    equation=equation,
                    chebis=chebis,
                    ecs=ec)
        rhea.save()

    return rhea.equation, rhea.chebis, rhea.ecs

def ec_to_ec3(list_ec):
    new_ecs = []
    for ec in list_ec:
        if ec.count('.') == 3:
            ec = ec[0:-3]
        if ec not in new_ecs:
            new_ecs.append(ec)
    return new_ecs

def chebi_to_smi(chebi_code):
    che =  Chebi.objects(code=chebi_code).first()
    if che is None:
        try:
            res = bioservices.ChEBI().getCompleteEntity(chebi_code)
            smi = res.smiles
            che = Chebi(code=chebi_code,
                        smiles=smi)
            che.save()

        except:
            print('ERROR getting ChEBI SMILES')
            return ''

    return che.smiles

def parse_equation(equation, chebis):
    lhs, rhs = equation.split(' = ')
    lhs = lhs.split(' + ')
    rhs = rhs.split(' + ')
    equation_split = lhs + rhs
    smi_map = {}
    if len(equation_split) == len(chebis):
        chebi_map = {k: v for k, v in zip(equation_split, chebis)}
        for k, v in chebi_map.items():
            smi_map[k] = chebi_to_smi(v)
    else:
        smi_map = {k: '' for k in equation_split}
    return lhs, rhs, smi_map


def make_rxn_smi(lhs, rhs, smi_map, rxnSize=(500, 200)):
    rxn_smi = ""
    for r in lhs:
        if rxn_smi != "":
            rxn_smi += "."

        rxn_smi += smi_map[r]

    rxn_smi += '>>'

    for r in rhs:
        if rxn_smi != "":
            rxn_smi += "."

        rxn_smi += smi_map[r]

    svg = smiles_rxn_to_svg(rxn_smi, rxnSize=rxnSize)

    return rxn_smi, svg

def make_smi_dicts(lhs, rhs, smi_map, molSize=(200,200)):
    lhs_dict = {}
    for r in lhs:
        smi = smi_map[r]
        try:
            svg = images.moltosvg_url(Chem.MolFromSmiles(smi), molSize=molSize)
        except:
            svg = ''
        lhs_dict[r] = [r, smi, svg]

    rhs_dict = {}
    for r in rhs:
        smi = smi_map[r]
        try:
            svg = images.moltosvg_url(Chem.MolFromSmiles(smi), molSize=molSize)
        except:
            svg = ''
        rhs_dict[r] = [r, smi, svg]

    return lhs_dict, rhs_dict

def get_rhea_imgs(rhea_id, molSize=(200,200)):
    equation, chebis, ec = get_rhea(rhea_id)
    lhs, rhs, smi_map = parse_equation(equation, chebis)
    lhs_dict, rhs_dict = make_smi_dicts(lhs, rhs, smi_map, molSize=molSize)
    return lhs_dict, rhs_dict

def get_rheas_for_enzyme_type(enzyme_type):
    et_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    rheas = UniRef50.objects(enzyme_type=et_obj).distinct('rhea')
    return rheas

def get_rhea_dicts(rheas, molSize=(125,90)):
    rhea_equations = {}
    rhea_imgs = {}
    for r in rheas:
        equation, chebis, ec = get_rhea(r)
        lhs, rhs, smi_map = parse_equation(equation, chebis)
        lhs_dict, rhs_dict = make_smi_dicts(lhs, rhs, smi_map, molSize=molSize)
        rhea_equations[r] = [equation, ec]
        rhea_imgs[r] = [lhs_dict, rhs_dict]
    return rhea_equations, rhea_imgs

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    #rhea_id = 'RHEA:15133'
    rhea_id = 'RHEA-COMP:9561'

    equation, chebis, ec = get_rhea(rhea_id)
    print(equation)
    print(chebis)
    lhs, rhs, smi_map = parse_equation(equation, chebis)
    print(smi_map)
    print(lhs)
    lhs_dict, rhs_dict = make_smi_dicts(lhs, rhs, smi_map)



