import cirpy

def get_smi_from_name(name):
    smi = cirpy.resolve(name, 'smiles')
    return smi


if __name__ == '__main__':
    name = 'Benzoic acid'
    name_wrong = "this is not a molecule"
    print(get_smi_from_name(name))
    print(get_smi_from_name(name_wrong))

