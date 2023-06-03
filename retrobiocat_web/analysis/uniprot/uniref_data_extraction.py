import collections


def get_cluster_name(xml_dict):
    if 'UniRef' in xml_dict:
        return xml_dict['UniRef']['entry']['@id']
    return ''

def check_id_match(original_name, xml_dict):
    retrieved_name = get_cluster_name(xml_dict)
    if retrieved_name != original_name:
        return False
    return True

def strip_uniref_name(name):
    name = name.replace('UniRef50_', '')
    name = name.replace('UniRef90_', '')
    name = name.replace('UniRef100_', '')
    return name

def get_uniprot_id(cluster_name):
    if cluster_name is None:
        return None

    uniprot_id = strip_uniref_name(cluster_name)
    if uniprot_id[0:2] == 'UP':
        uniprot_id = ""
    return uniprot_id

def get_uniref_members(xml_dict):
    list_to_process = []
    uniref_90 = set()
    uniref_100 = set()
    uniprot_dict = dict()
    uniprot_kb = set()

    if xml_dict != {}:
        list_to_process.append(xml_dict['UniRef']['entry']['representativeMember']['dbReference']['property'])
        if 'member' in xml_dict['UniRef']['entry']:
            for member in xml_dict['UniRef']['entry']['member']:
                if type(member) == collections.OrderedDict:
                    list_to_process.append(member['dbReference']['property'])

    # print(list_to_process[0])

    for member_list in list_to_process:
        uniprot_acc = ""
        org = ""
        protein_name = ""
        for rep in member_list:
            if rep.get('@type', '') == 'UniRef100 ID':
                uniref_100.add(rep['@value'])
            elif rep.get('@type', '') == 'UniRef90 ID':
                uniref_90.add(rep['@value'])
            elif rep.get('@type', '') == 'UniProtKB accession':
                uniprot_acc = rep['@value']
                uniprot_kb.add(rep['@value'])
            elif rep.get('@type', '') == 'protein name':
                protein_name = rep['@value']
            elif rep.get('@type', '') == 'source organism':
                org = rep['@value']

        if uniprot_acc != '':
            uniprot_dict[uniprot_acc] = [protein_name, org]

    return uniref_90, uniref_100, uniprot_dict, uniprot_kb