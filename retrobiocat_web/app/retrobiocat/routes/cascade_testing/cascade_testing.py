from retrobiocat_web.app.retrobiocat import bp
from flask import render_template, request
from flask_security import roles_required
from retrobiocat_web.mongo.models.retro_tests import TestCascadeRun
from rdkit import Chem
from retrobiocat_web.analysis.drawing import images


def get_cascade_tests_table():
    query = TestCascadeRun.objects()
    table = []
    for run in query:
        data = {'id': str(run.id),
                'run_name': run.run_name,
                'date': str(run.date),
                'mode': run.settings['pathway_settings']['mode'],
                'num_results': len(run.results),
                'top_5': run.top_5_score,
                'top_25': run.top_25_score,
                'top_200': run.top_200_score}
        table.append(data)

    return table

def smiles_to_svg(activity_data):
    for i, record in enumerate(activity_data):
        if activity_data[i]["substrate_1_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["substrate_1_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["substrate_1_smiles"] = url

        if activity_data[i]["substrate_2_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["substrate_2_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["substrate_2_smiles"] = url

        if activity_data[i]["product_1_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["product_1_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["product_1_smiles"] = url

    return activity_data

def get_results_table(run):
    table = []
    for result in run.results:
        cascade = result.test_cascade

        target_mol = Chem.MolFromSmiles(cascade.target_smiles)
        target_img = images.moltosvg_url(target_mol)

        substrates_smi = ""
        for smi in list(cascade.starting_smiles):
            substrates_smi += smi
            substrates_smi += '.'
        substrates_smi = substrates_smi[:-1]
        substrates_mol = Chem.MolFromSmiles(substrates_smi)
        substrates_img = images.moltosvg_url(substrates_mol)

        data = {'id': str(result.id),
                'test_id': str(cascade.id),
                'test_name': cascade.name,
                'test_type': cascade.test_type,
                'doi': cascade.doi,
                'target_smiles': target_img,
                'substrates': substrates_img,
                'enzymes': cascade.enzymes,
                'ranking_any': result.ranking_any or '',
                'ranking_specific': result.ranking_specific or '',
                'time_taken': result.time_taken,
                'total_pathways': result.total_pathways}
        table.append(data)

    return table

@roles_required('admin')
@bp.route('/cascade_tests', methods=['GET'])
def cascade_tests():
    table = get_cascade_tests_table()
    return render_template('cascade_testing/cascade_tests.html', cascade_tests_table=table)

@roles_required('admin')
@bp.route('/cascade_result/', methods=['GET'])
def cascade_result():
    run_id = request.args.get("run_id")
    run = TestCascadeRun.objects(id=run_id).first().select_related()
    table = get_results_table(run)
    run_name = run.run_name
    run_date = run.date
    return render_template('cascade_testing/cascade_results.html',
                           cascade_results_table=table,
                           run_name=run_name, run_date=run_date)



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    table = get_cascade_tests_table()
    print(table)
