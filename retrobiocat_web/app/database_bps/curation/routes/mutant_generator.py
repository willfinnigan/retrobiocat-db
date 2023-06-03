from flask import render_template
from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.mongo.model_queries.sequence_queries import seq_obj_from_name


@bp.route('/mutant_generator', methods=['GET'])
def mutant_generator():
    return render_template('mutant_generator.html', enzyme_name="", protein_sequence="")

@bp.route('/mutant_generator_load/<parent_enzyme_name>', methods=['GET'])
def mutant_generator_load(parent_enzyme_name):
    seq_obj = seq_obj_from_name(parent_enzyme_name)
    if seq_obj is None:
        protein_sequence = ""
    elif seq_obj.sequence is None:
        protein_sequence = ""
    else:
        protein_sequence = seq_obj.sequence

    return render_template('mutant_generator.html', enzyme_name=parent_enzyme_name, protein_sequence=protein_sequence)