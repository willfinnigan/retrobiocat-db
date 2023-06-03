from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, request, jsonify, session, current_app, redirect, make_response, send_file
from retrobiocat_web.mongo.models.biocatdb_models import UniRef50
from flask_security import roles_required
import pandas as pd
import json


@roles_required('admin')
@bp.route("/_download_selected_unireps", methods=["POST"])
def download_selected_unireps():
    unirep_ids = json.loads(request.form['unirep_ids'])
    unirefs = UniRef50.objects(enzyme_name__in=unirep_ids)
    seqs = []
    for uniref in unirefs:
        seqs.append(uniref.sequence)

    csv = "name, sequence \n"
    fasta = ""
    for seq_id, seq in zip(unirep_ids, seqs):
        csv += f"{seq_id}, {seq} \n"
        fasta += f">{seq_id} \n"
        fasta += f"{seq} \n"

    return jsonify(result={'csv': csv, 'fasta': fasta})

