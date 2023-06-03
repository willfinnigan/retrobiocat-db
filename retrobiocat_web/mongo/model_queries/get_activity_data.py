from retrobiocat_web.mongo.functions import process_get_activity_data
from retrobiocat_web.mongo.models.biocatdb_models import Activity, Paper
import mongoengine as db

MONGO_COLS = ['reaction',
           'enzyme_type',
           'enzyme_name',
           'short_citation',
           'html_doi',
           'cascade_num',
           'substrate_1_smiles',
           'substrate_2_smiles',
           'product_1_smiles',
           'temperature',
           'ph',
           'solvent',
           'other_conditions',
           'notes',
           'reaction_vol',
           'formulation',
           'biocat_conc',
           'kcat',
           'km',
           'mw',
           'substrate_1_conc',
           'substrate_2_conc',
           'specific_activity',
           'conversion',
           'conversion_time',
           'categorical',
           'binary',
           'selectivity',
           'paper',
           'id']

def get_activity_data(enzyme_type=None, reaction=None,
                      paper_id=None, enzyme_names=None,
                      formulation=None, solvent=None,
                      ph_left=None, ph_right=None,
                      temp_left=None, temp_right=None,
                      smi_col='product_1_smiles', smiles=None,
                      only_reviewed=True,
                      include_negative=True,
                      include_auto_generated=False,
                      as_pymongo=False):

    if include_auto_generated == False:
        aq = db.Q(auto_generated=False)
    else:
        aq = db.Q()

    if include_negative == True:
        inq = db.Q()
    else:
        inq = db.Q(binary=True)

    etq = db.Q()
    if enzyme_type is not None and enzyme_type.lower() != 'all':
        etq = db.Q(enzyme_type=enzyme_type)

    rq = db.Q()
    if reaction is not None and reaction.lower() != 'all':
        rq = db.Q(reaction=reaction)

    sol_q = db.Q()
    if solvent is not None:
        sol_q = db.Q(solvent=solvent)

    form_q = db.Q()
    if formulation is not None:
        form_q = db.Q(formulation=formulation)

    ph_l_q = db.Q()
    if ph_left is not None:
        ph_l_q = db.Q(ph__gte=float(ph_left))

    ph_r_q = db.Q()
    if ph_right is not None:
        ph_r_q = db.Q(ph__lte=float(ph_right))

    temp_l_q = db.Q()
    if temp_left is not None:
        temp_l_q = db.Q(temperature__gte=float(temp_left))

    temp_r_q = db.Q()
    if temp_right is not None:
        temp_r_q = db.Q(temperature__lte=float(temp_right))

    if paper_id == None:
        pq = db.Q()
    else:
        paper_query = Paper.objects(id=paper_id)
        if len(paper_query) == 0:
            pq = db.Q(paper=None)
        else:
            paper = paper_query[0]
            pq = db.Q(paper=paper)

    if only_reviewed == True:
        orq = db.Q(reviewed=True)
    else:
        orq = db.Q()

    if enzyme_names == None:
        enq = db.Q()
    else:
        if type(enzyme_names) == str:
            enzyme_names.replace(", ", ",")
            enzyme_names = enzyme_names.split(',')
        enq = db.Q(enzyme_name__in=enzyme_names)

    if smiles == None:
        smiq = db.Q()
    else:
        if type(smiles) == str:
            smiles.replace(", ", ",")
            smiles = smiles.split(',')

        if smi_col == 'product_1_smiles':
            smiq = db.Q(product_1_smiles__in=smiles)
        elif smi_col == 'substrate_1_smiles':
            smiq = db.Q(substrate_1_smiles__in=smiles)
        elif smi_col == 'substrate_2_smiles':
            smiq = db.Q(substrate_2_smiles__in=smiles)
        else:
            smiq = db.Q()

    if as_pymongo is False:
        return Activity.objects(etq & rq & orq & pq & enq & smiq & inq & aq & sol_q & form_q & ph_l_q & ph_r_q & temp_l_q & temp_r_q).only(*MONGO_COLS)
    else:
        activity_data = list(Activity.objects(etq & rq & orq & pq & enq & smiq & inq & aq & sol_q & form_q & ph_l_q & ph_r_q & temp_l_q & temp_r_q).only(
            *MONGO_COLS).as_pymongo())
        activity_data = process_get_activity_data.process_activity_data(activity_data)
        return activity_data

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    smiles = ['CCCC=O']
    data = get_activity_data(smiles=smiles, as_pymongo=True)

    print(data)