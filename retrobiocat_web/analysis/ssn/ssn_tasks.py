from retrobiocat_web.analysis.sequence_analysis.all_by_all_blast import AllByAllBlaster
from flask import current_app
from retrobiocat_web.analysis.ssn.ssn_main import SSN
from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser
from retrobiocat_web.analysis.ssn.ssn_quickload import SSN_quickload
from retrobiocat_web.analysis.ssn.ssn_precalculate import SSN_Cluster_Precalculator
from pathlib import Path
import json
from retrobiocat_web.analysis import analysis_paths


def task_expand_ssn(enzyme_type, log_level=1, max_num=200):
    current_app.app_context().push()

    aba_blaster = AllByAllBlaster(enzyme_type, log_level=log_level)
    aba_blaster.make_blast_db()

    ssn = SSN(enzyme_type, aba_blaster=aba_blaster, log_level=log_level)
    ssn.load_sqlite()
    ssn.set_status('Checking SSN')
    ssn.remove_nonexisting_seqs()
    ssn.remove_seqs_marked_with_no_alignments()

    biocatdb_seqs = ssn.nodes_not_present(only_biocatdb=True, max_num=max_num)
    if len(biocatdb_seqs) != 0:
        ssn.clear_position_information()
        ssn.set_status('Adding and aligning BioCatDB sequences')
        ssn.add_multiple_proteins(biocatdb_seqs)
        ssn.save_sqlite()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    need_alignments = ssn.nodes_need_alignments(max_num=max_num)
    if len(need_alignments) != 0:
        ssn.clear_position_information()
        ssn.set_status('Aligning sequences in SSN')
        ssn.add_multiple_proteins(need_alignments)
        ssn.save_sqlite()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    not_present = ssn.nodes_not_present(max_num=max_num)
    if len(not_present) != 0:
        ssn.clear_position_information()
        ssn.set_status('Adding UniRef sequences which are not yet present')
        ssn.add_multiple_proteins(not_present)
        ssn.save_sqlite()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)

        return

    if ssn.db_object.precalc_status != 'Complete':
        print('SSN prec calc status is not complete - running positioning')
        ssn.db_object.precalc_status = 'Running'
        ssn.save_sqlite()
        ssn.set_status('Precalculating visualisations')
        current_app.preprocess_queue.enqueue(precalculate_job, enzyme_type)

    else:
        print(f'{ssn.db_object.precalc_status}')
        ssn.set_status('Complete')
        print(f'- SSN CONSTRUCTION FOR {enzyme_type} IS COMPLETE -')
        ssn.db_object.save()

def new_expand_ssn_job(enzyme_type):
    ssn = SSN(enzyme_type)
    if ssn.db_object.status != 'Complete':
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)

def precalculate_job(enzyme_type):
    ssn = SSN(enzyme_type)
    ssn.load_sqlite(include_mutants=False, only_biocatdb=False)

    ssn_precalc = SSN_Cluster_Precalculator(ssn)

    ssn.db_object.precalc_status = 'Precalculating SSN positions'
    ssn.db_object.status = 'Precalculating SSN positions'
    ssn.db_object.save()

    num_nodes = len(list(ssn.graph.nodes))
    if num_nodes > 3000:
        num = 1
    elif num_nodes > 1000:
        num = 5
    else:
        num = 20

    if len(list(ssn.db_object.num_at_alignment_score.values())) == 0:
        print('No existing % identity data, starting at alignment score 40')
        ssn.clear_position_information()
        ssn_precalc.start = 40
        current_num_clusters = 0
    else:
        start_list = [int(s) for s in list(ssn.db_object.identity_at_alignment_score.keys())]
        current_num_clusters = max(list(ssn.db_object.num_at_alignment_score.values()))
        ssn_precalc.start = max(start_list) + 5

    print(f"Start = {ssn_precalc.start}")
    precalculated_nodes, cluster_numbers, identity_at_score = ssn_precalc.precalulate(num=num, current_num_clusters=current_num_clusters)

    if len(precalculated_nodes) == 0:
        ssn.db_object.precalc_status = 'Complete'
        ssn.db_object.status = 'Complete'
        ssn.db_object.save()
        print('Precalc complete - checking SSN again')
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)
    else:
        #ssn.db_object.precalculated_vis.update(precalculated_nodes)
        ssn.db_object.num_at_alignment_score.update(cluster_numbers)
        ssn.db_object.identity_at_alignment_score.update(identity_at_score)
        ssn.db_object.save()
        save_precalculated_vis(enzyme_type, precalculated_nodes)
        current_app.preprocess_queue.enqueue(precalculate_job, enzyme_type)

def save_precalculated_vis(enzyme_type, precalculated_nodes):
    folder = f'{analysis_paths.ssn_folder}/{enzyme_type}/positions'
    Path(folder).mkdir(parents=True, exist_ok=True)
    for alignment_score in precalculated_nodes:
        filepath = f'{folder}/{alignment_score}_positions.json'

        with open(filepath, 'wb') as outfile:
            outfile.write(json.dumps(precalculated_nodes[alignment_score]).encode("utf-8"))

def load_precalculated_vis(enzyme_type, alignment_score):
    #try:
    folder = f'{analysis_paths.ssn_folder}/{enzyme_type}/positions'
    filepath = f'{folder}/{alignment_score}_positions.json'
    nodes = json.load(open(filepath))
    return nodes

    #except:
        #print('Error loading positions')

    #return None

def load_all_edges(enzyme_type, alignment_score):
    ql = SSN_quickload(enzyme_type, log_level=0)
    ql.load_df()
    edges = ql.get_all_edges(alignment_score)
    return edges

def clear_position_data_task_and_recalc(enzyme_type):
    ssn = SSN(enzyme_type)
    ssn.clear_position_information()
    current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)

def remove_sequence(enzyme_type, enzyme_name):
    ssn = SSN(enzyme_type)
    ssn.load_sqlite()

    if len(list(ssn.graph.nodes)) != 0:
        if enzyme_name in list(ssn.graph.nodes):
            ssn.graph.nodes.remove(enzyme_name)
            ssn.save_sqlite()
            current_app.preprocess_queue.enqueue(precalculate_job, enzyme_type)