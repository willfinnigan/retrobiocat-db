import time
import zipfile

from retrobiocat_web.app.database_bps.db_admin import bp
from flask import make_response, send_file
from flask_security import roles_required

from retrobiocat_web.app.database_bps.db_admin.functions.reaction_rule_functions import yaml_conversion
from retrobiocat_web.app.database_bps.tables.get_table_data.get_papers_table_data import process_papers_to_table
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import Sequence
from retrobiocat_web.mongo.functions.mongo_dump import execute_mongo_dump
from retrobiocat_web.mongo.model_queries import paper_queries, specificity_data_query
import yaml
import json
import pandas as pd
import shutil

@bp.route('/download_rxn_yaml', methods=['GET'])
@roles_required('admin')
def download_rxn_yaml():
    all_reactions = Reaction.objects()
    yaml_dict = {}
    for reaction in all_reactions:
        yaml_dict = yaml_conversion.reaction_to_yaml_dict(yaml_dict, reaction)

    yaml_json = json.dumps(yaml_dict)
    yaml_dict = json.loads(yaml_json)
    resp = make_response(yaml.dump(yaml_dict))
    resp.headers["Content-Disposition"] = "attachment; filename=rxns_yaml.yaml"
    return resp

@bp.route('/download_biocatdb', methods=['GET'])
@roles_required('admin')
def download_biocatdb():
    spec_df = specificity_data_query.query_specificity_data(['All'], ['All'])
    resp = make_response(spec_df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=biocatdb_2.csv"
    return resp

@bp.route('/download_papers', methods=['GET'])
@roles_required('admin')
def download_papers():
    papers = paper_queries.get_all_papers()
    papers_list = process_papers_to_table(papers)
    papers_df = pd.DataFrame(papers_list)
    resp = make_response(papers_df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=papers.csv"
    return resp

@bp.route('/download_sequences', methods=['GET'])
@roles_required('admin')
def download_sequences():
    list_cols = ["enzyme_type", "enzyme_name", "other_names", "other_names_data", "n_tag", "sequence",
                 "c_tag", "sequence_unavailable", "accession", "other_identifiers", "structure", "pdb", "mutant_of",
                 "notes", "reviewed"]
    sequences = Sequence.objects().only(*list_cols).as_pymongo()
    df = pd.DataFrame(list(sequences))
    resp = make_response(df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=sequences.csv"
    return resp


@bp.route('/download_the_mongo_dump', methods=['GET'])
@roles_required('admin')
def download_the_mongo_dump():
    print('download_mongo_dump')
    output_path = execute_mongo_dump()
    return send_file(output_path, as_attachment=True)

