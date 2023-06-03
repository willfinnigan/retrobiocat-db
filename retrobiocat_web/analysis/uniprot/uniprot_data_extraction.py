import collections

def get_pfams(xml_dict):
    """ Return a list of the pfam id's in a uniprot or uniref representative sequence """

    p_data = xml_dict['uniprot']['entry']['dbReference']
    pfams = dict()
    for list_item in p_data:
        if type(list_item) != str:
            if list_item.get("@type", '') == 'Pfam':
                pfam_id = ''
                pfam_name = ''
                if '@id' in list_item:
                    pfam_id = list_item['@id']
                if 'property' in list_item:
                    for property_item in list_item['property']:
                        if property_item.get('@type', '') == 'entry name':
                            if '@value' in property_item:
                                pfam_name = property_item['@value']

                if pfam_id != '':
                    pfams[pfam_id] = pfam_name

    return pfams

def get_pdbs(xml_dict):
    """ Return a list of the pdbs in a uniprot or uniref representative sequence"""

    p_data = xml_dict['uniprot']['entry']['dbReference']
    pdbs = []
    for list_item in p_data:
        if type(list_item) != str:
            if list_item.get("@type", '') == 'PDB':
                if '@id' in list_item:
                    pdbs.append(list_item['@id'])

    return pdbs

def get_rhea(xml_dict, ignore_rhea_comp=True):
    """ Return a list of the rhea reactions in a uniprot or uniref representative sequence"""

    if 'comment' not in xml_dict['uniprot']['entry']:
        return []

    r_data = xml_dict['uniprot']['entry']['comment']
    rhea = []
    for list_item in r_data:
        if type(list_item) != str:
            if list_item.get("@type", '') == 'catalytic activity':
                if 'reaction' in list_item:
                    if 'dbReference' in list_item['reaction']:
                        for db_ref in list_item['reaction']['dbReference']:
                            if type(db_ref) == collections.OrderedDict:
                                if db_ref.get('@type', '') == 'Rhea':
                                    rhea_id = db_ref.get('@id', '')
                                    if (ignore_rhea_comp == False or 'RHEA-COMP' not in rhea_id) and rhea_id != '':
                                        rhea.append(rhea_id)

    return rhea

def print_xml(xml_dict):
    for key, value in xml_dict['uniprot']['entry'].items():
        print(f"-- {key} --")
        print(value)
        print()

def is_swissprot(xml_dict):
    if xml_dict != {}:
        sp = xml_dict['uniprot']['entry']['@dataset']
        if sp == "Swiss-Prot":
            return True
    return False