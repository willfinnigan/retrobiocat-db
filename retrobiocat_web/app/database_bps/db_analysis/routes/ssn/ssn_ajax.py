from retrobiocat_web.analysis.uniprot.uniref_data_extraction import get_uniprot_id
from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, request, jsonify, session, current_app, redirect, make_response, send_file
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, UniRef50, SSN_record, Sequence
import mongoengine as db
from flask_security import roles_required, current_user

from retrobiocat_web.analysis.ssn.ssn_main import SSN
from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser
from retrobiocat_web.analysis.ssn.ssn_quickload import SSN_quickload
from retrobiocat_web.analysis.uniprot_and_web import rhea_api
import json
from distutils.util import strtobool
import yaml


@bp.route("/_load_uniref_data", methods=["POST"])
def load_uniref_data():
    name = request.form['name']
    enzyme_type = request.form['enzyme_type']
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    et = db.Q(enzyme_type=enzyme_type_obj)
    nq = db.Q(enzyme_name=name)

    uniref = UniRef50.objects(et & nq).first()

    result = {'rep_seq_name': uniref.protein_name,
              'rep_seq_organism': uniref.tax,
              'rep_seq_uniprot_id': get_uniprot_id(uniref.cluster_name),
              'cluster_id': uniref.cluster_name,
              'sp': uniref.sp_annotated,
              'num_uni90': uniref.num_uni90,
              'num_uni100': uniref.num_uni100,
              'num_uniprot': uniref.num_uniprot,
              'pfam_object': uniref.pfams,
              'pdbs': uniref.pdbs,
              'rheas': uniref.rhea}

    return jsonify(result=result)

@bp.route("/_edge_ajax", methods=["POST"])
def edge_ajax():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    selected_node = request.form['selected_node']

    ql = SSN_quickload(enzyme_type, log_level=0)
    ql.load_df()
    edges = ql.get_edges(selected_node, alignment_score)

    result = {'edges': edges}
    return jsonify(result=result)

@bp.route("/_connected_nodes_ajax", methods=["POST"])
def connected_nodes_ajax():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    selected_nodes = json.loads(request.form['selected_nodes'])

    ql = SSN_quickload(enzyme_type, log_level=1)
    ql.load_df()
    nodes = ql.get_connected_nodes(selected_nodes, alignment_score)

    result = {'nodes': nodes}
    return jsonify(result=result)

@bp.route("/_get_clusters", methods=["POST"])
def get_clusters():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    only_biocatdb = bool(strtobool(request.form['only_biocatdb']))
    hide_mutants = bool(strtobool(request.form['hide_mutants']))

    ssn = SSN(enzyme_type)
    ssn.load_sqlite(include_mutants=not hide_mutants, only_biocatdb=only_biocatdb)
    without_uniref, with_uniref = ssn.get_clusters(alignment_score)

    result = {'with_uniref': with_uniref,
              'without_uniref': without_uniref}
    return jsonify(result=result)

@bp.route("/_load_all_edges_ajax", methods=["POST"])
def load_all_edges_ajax():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    nodes = json.loads(request.form['list_nodes'])

    ql = SSN_quickload(enzyme_type, log_level=0)
    ql.load_df()
    edges = ql.get_multiple_edges(nodes, alignment_score)

    result = {'edges': edges}
    return jsonify(result=result)

@bp.route("/_download_cytoscape", methods=["POST"])
@roles_required('experimental')
def download_cytoscape():
    enzyme_type = request.form['enzyme_type']
    nodes = json.loads(request.form['list_nodes'])
    alignment_score = int(request.form['alignment_score'])
    only_biocatdb = bool(strtobool(request.form['only_biocatdb']))
    hide_mutants = bool(strtobool(request.form['hide_mutants']))

    ssn = SSN(enzyme_type)
    ssn.load_sqlite(include_mutants=not hide_mutants, only_biocatdb=only_biocatdb)
    cytoscape = ssn.send_to_cytoscape(alignment_score, selected_nodes=nodes)

    #resp = make_response(yaml.dump(cytoscape))
    #resp.headers["Content-Disposition"] = "attachment; filename=ssn_file.cyjs"
    return jsonify(result={'cytoscape': cytoscape})

@bp.route("/_get_rhea", methods=["POST"])
def get_rhea():
    rhea_code = request.form['rhea_code']
    lhs_dict, rhs_dict = rhea_api.get_rhea_imgs(rhea_code, molSize=(125,90))
    return jsonify(result={'lhs': lhs_dict,
                           'rhs': rhs_dict})

@bp.route("/_get_rhea_equation", methods=["POST"])
def get_rhea_equation():
    rhea_code = request.form['rhea_code']
    equation, chebis, ec = rhea_api.get_rhea(rhea_code)
    return jsonify(result={'equation': equation})


