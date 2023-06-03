import os
import uuid
from datetime import datetime
from typing import List

from retrobiocat_web.app.retrobiocat3.functions.load_save_datastore import DATASTORE_PATH, datastore_json_filepath
from retrobiocat_web.mongo.model_queries.retro3_queries import get_project, get_search, get_pathway_hash
from retrobiocat_web.mongo.models.retro3_models import Project, SearchResult
from retrobiocat_web.retro3.pathway.pathway import Pathway
from pathlib import Path

def create_new_search(project_obj: Project, source_mol_config_dict, mcts_config_dict):
    new_search = SearchResult(mcts_config_dict=mcts_config_dict,
                              source_mol_config_dict=source_mol_config_dict)
    project_obj.search_results.append(new_search)
    new_search.save()
    project_obj.save()
    return new_search

def add_search_result(search_obj: SearchResult,
                      solved_pathways: List[Pathway],
                      enzyme_pathways_pathways: List[Pathway],
                      all_pathways: List[Pathway],
                      network_size: int,
                      max_pathways=10000):


    if len(all_pathways) > max_pathways:
        all_pathways = all_pathways[0:max_pathways]

    search_obj.solved_pathways = [{'reaction_ids': p.save()} for p in solved_pathways]
    search_obj.num_solved_pathways = len(solved_pathways)
    search_obj.all_pathways = [{'reaction_ids': p.save()} for p in all_pathways]
    search_obj.num_all_pathways = len(all_pathways)
    search_obj.solved_pathways_with_enzymes = [{'reaction_ids': p.save()} for p in enzyme_pathways_pathways]
    search_obj.num_solved_pathways_solved_pathways_with_enzymes = len(enzyme_pathways_pathways)
    search_obj.network_size = network_size
    search_obj.datetime_completed = datetime.now()
    search_obj.save()

def update_search_status(search_obj: SearchResult, status: str):
    search_obj.status = status
    search_obj.save()

def update_project_reaction_number(project_obj: Project, num_reactions: int):
    project_obj.num_reactions = num_reactions
    project_obj.last_search_completed = datetime.now()
    project_obj.save()

def create_new_project_objects(target_smi, expansion_config_dict, source_mol_config_dict, user, project_name):

    new_project = Project(target_smi=target_smi,
                          project_name=project_name,
                          expansion_config_dict=expansion_config_dict,
                          source_mol_config_dict=source_mol_config_dict,
                          owner=user)
    new_project.save()

    return new_project

def delete_retro3_project(project_id):
    try:
        filepath = datastore_json_filepath(project_id)
        file_to_rem = Path(filepath)
        file_to_rem.unlink(missing_ok=True)
        project_obj = get_project(project_id)
        project_obj.delete()
        return []
    except Exception as e:
        return [str(e)]

def delete_retro3_search(project_id, search_id):
    try:
        project_obj = get_project(project_id)
        search_obj = get_search(search_id)
        if search_obj in project_obj.search_results:
            project_obj.search_results.remove(search_obj)
        project_obj.save()
        search_obj.delete()
        return []
    except Exception as e:
        return [str(e)]

def restart_retro3_search(search_id):
    try:
        search_obj = get_search(search_id)
        search_obj.solved_pathways_with_enzymes = []
        search_obj.datetime = datetime.utcnow()
        search_obj.datetime_completed = None
        search_obj.status = 'Queued'
        search_obj.network_size = None
        search_obj.all_pathways = None
        search_obj.num_all_pathways = None
        search_obj.solved_pathways = None
        search_obj.num_solved_pathways = None
        search_obj.solved_pathways_with_enzymes = None
        search_obj.num_solved_pathways_solved_pathways_with_enzymes = None
        search_obj.save()
        return []
    except Exception as e:
        return [str(e)]



def save_pathway_to_project(project_id, reaction_ids):
    project_obj = get_project(project_id)
    reaction_ids = sorted(reaction_ids)

    if reaction_ids in project_obj.saved_pathways:
        return ['Error, reaction already saved']

    project_obj.saved_pathways.append(reaction_ids)
    project_obj.save()
    print('pathway saved')
    return []

def unsave_pathway_to_project(project_id, reaction_ids):
    project_obj = get_project(project_id)
    reaction_ids = sorted(reaction_ids)

    if reaction_ids in project_obj.saved_pathways:
        project_obj.saved_pathways.remove(reaction_ids)
        project_obj.save()
        return []

    return ['Error, pathway not saved to begin with']


def add_reaction_to_block_list(project, reaction_id, reaction_svg):
    if reaction_id not in project.blocked_reactions:
        project.blocked_reactions[reaction_id] = reaction_svg
        project.save()

def remove_reaction_from_block_list(project, reaction_id):
    if reaction_id in project.blocked_reactions:
        project.blocked_reactions.pop(reaction_id)
        project.save()

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    Project.drop_collection()
    SearchResult.drop_collection()
