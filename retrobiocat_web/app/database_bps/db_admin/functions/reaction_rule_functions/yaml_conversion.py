import yaml
import collections
import json


def setup_yaml():
  represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
  yaml.add_representer(collections.OrderedDict, represent_dict_order)

setup_yaml()


def steps_to_names(rxn_steps):
    list_of_list_of_names = []
    for group_steps in rxn_steps:
        list_names = []
        for step in group_steps:
            list_names.append(step.name)
        list_of_list_of_names.append(list_names)
    print(list_of_list_of_names)
    return list_of_list_of_names

def reaction_to_yaml_dict(yaml_dict, reaction):
    name = reaction.name
    yaml_dict[name] = {}

    yaml_dict[name]['smarts'] = list(reaction.smarts)
    yaml_dict[name]['steps'] = steps_to_names(reaction.steps)
    yaml_dict[name]['enzymes'] = dict(reaction.cofactors)
    yaml_dict[name]['positive_tests'] = list(reaction.positive_tests)
    yaml_dict[name]['negative_tests'] = list(reaction.negative_tests)
    yaml_dict[name]['type'] = reaction.type
    yaml_dict[name]['product_seeds'] = list(reaction.product_seeds)
    yaml_dict[name]['substrate_1_seeds'] = list(reaction.substrate_1_seeds)
    yaml_dict[name]['substrate_2_seeds'] = list(reaction.substrate_2_seeds)
    if reaction.experimental == True:
        yaml_dict[name]['experimental'] = True

    if reaction.two_step == True:
        yaml_dict[name]['two_step'] = True

    if reaction.requires_absence_of_water == True:
        yaml_dict[name]['requires_absence_of_water'] = True

    yaml_json = json.dumps(yaml_dict)
    yaml_dict = json.loads(yaml_json)

    return yaml_dict

def yaml_dict_to_yaml(yaml_dict):
    yaml_json = json.dumps(yaml_dict)
    yaml_dict = json.loads(yaml_json)
    return yaml.dump(yaml_dict)


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    from retrobiocat_web.mongo.models.reaction_models import Reaction

    redam = Reaction.objects(name='Reductive amination').select_related()[0]
    yaml_dict = reaction_to_yaml_dict({}, redam)
    yaml = yaml_dict_to_yaml(yaml_dict)
    print(yaml)