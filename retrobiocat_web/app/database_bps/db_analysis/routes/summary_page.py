from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.analysis.reaction_summary.enzyme_count_graph import EnzymeCountGraphCreator
from retrobiocat_web.analysis.reaction_summary.formulation_counts import FormulationCountGraphCreator
from retrobiocat_web.analysis.reaction_summary.mol_count_graph import MolCountGraphCreator
from retrobiocat_web.analysis.reaction_summary.mol_descriptor_graphs import MolDescriptorHistrogramCreator
from retrobiocat_web.analysis.reaction_summary.ph_histogram import pHHistrogramCreator
from retrobiocat_web.analysis.reaction_summary.solvent_counts import SolventCountGraphCreator
from retrobiocat_web.analysis.reaction_summary.temperature_histogram import TemperatureHistrogramCreator
from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, request
from bokeh.embed import components
from retrobiocat_web.mongo.model_queries import reaction_queries, activity_queries, enzyme_type_queries

@bp.route('/summary_page/<enzyme_type>/<reaction>', methods=['GET'])
def summary_page(enzyme_type, reaction):

    enzyme_types = activity_queries.enzyme_types_from_activity(reaction=reaction, include_chemical=False)
    enzyme_type_full_names = enzyme_type_queries.get_enzyme_types_with_full_names(enzyme_types)
    reactions = activity_queries.reactions_from_activity(enzyme_type=enzyme_type, include_chemical=False)

    if enzyme_type not in enzyme_types:
        enzyme_type = 'All'
    if reaction not in reactions:
        reaction = 'All'

    remove_negative = None
    only_positive = ''
    args = request.args.to_dict()
    print(args)
    if 'only_positive'in args:
        remove_negative = 'all'
        only_positive = 'checked'

    data_query = get_data.DataQuery(reaction=reaction, enzyme_type=enzyme_type, log_level=1, remove_negative=remove_negative)

    enzyme_count_graph = EnzymeCountGraphCreator().create_graph(data_query)
    enzyme_script, enzyme_div = components(enzyme_count_graph)

    mol_count_graph = MolCountGraphCreator().create_graph(data_query)
    mol_script, mol_div = components(mol_count_graph)

    temp_hist_graph = TemperatureHistrogramCreator().create_graph(data_query)
    temp_script, temp_div = components(temp_hist_graph)

    ph_hist_graph = pHHistrogramCreator().create_graph(data_query)
    ph_script, ph_div = components(ph_hist_graph)

    solvent_count_graph = SolventCountGraphCreator().create_graph(data_query)
    solvent_script, solvent_div = components(solvent_count_graph)

    formulation_count_graph = FormulationCountGraphCreator().create_graph(data_query)
    formulation_script, formulation_div = components(formulation_count_graph)

    descriptor_count_graphs = MolDescriptorHistrogramCreator().create_graphs(data_query)
    descriptor_count_scripts_and_divs = []
    for graph in descriptor_count_graphs:
        script, div = components(graph)
        descriptor_count_scripts_and_divs.append([script, div])

    graphs = {'enzyme_graph': {'div': enzyme_div, 'script': enzyme_script},
              'mol_graph': {'div': mol_div, 'script': mol_script},
              'ph_graph': {'div': ph_div, 'script': ph_script},
              'temp_graph': {'div': temp_div, 'script': temp_script},
              'solvent_graph': {'div': solvent_div, 'script': solvent_script},
              'formulation_graph': {'div': formulation_div, 'script': formulation_script},
              'descriptor_graphs': descriptor_count_scripts_and_divs}

    return render_template('summary_page/summary_page.html',
                           graphs=graphs,
                           enzyme_type=enzyme_type, enzyme_types=enzyme_type_full_names,
                           reaction=reaction, reactions=reactions, only_positive=only_positive)

@bp.route('/reaction_conditions/<enzyme_type>/<reaction>', methods=['GET'])
def reaction_conditions(enzyme_type, reaction):

    enzyme_types = activity_queries.enzyme_types_from_activity(reaction=reaction, include_chemical=False)
    enzyme_type_full_names = enzyme_type_queries.get_enzyme_types_with_full_names(enzyme_types)
    reactions = activity_queries.reactions_from_activity(enzyme_type=enzyme_type, include_chemical=False)

    if enzyme_type not in enzyme_types:
        enzyme_type = 'All'
    if reaction not in reactions:
        reaction = 'All'

    remove_negative = None
    only_positive = ''
    args = request.args.to_dict()
    print(args)
    if 'only_positive'in args:
        remove_negative = 'all'
        only_positive = 'checked'

    data_query = get_data.DataQuery(reaction=reaction, enzyme_type=enzyme_type, log_level=1, remove_negative=remove_negative)

    temp_hist_graph = TemperatureHistrogramCreator().create_graph(data_query)
    temp_script, temp_div = components(temp_hist_graph)

    ph_hist_graph = pHHistrogramCreator().create_graph(data_query)
    ph_script, ph_div = components(ph_hist_graph)

    solvent_count_graph = SolventCountGraphCreator().create_graph(data_query)
    solvent_script, solvent_div = components(solvent_count_graph)

    formulation_count_graph = FormulationCountGraphCreator().create_graph(data_query)
    formulation_script, formulation_div = components(formulation_count_graph)

    graphs = {
              'ph_graph': {'div': ph_div, 'script': ph_script},
              'temp_graph': {'div': temp_div, 'script': temp_script},
              'solvent_graph': {'div': solvent_div, 'script': solvent_script},
              'formulation_graph': {'div': formulation_div, 'script': formulation_script},
              }

    return render_template('summary_page/reaction_conditions_page.html',
                           graphs=graphs,
                           enzyme_type=enzyme_type, enzyme_types=enzyme_type_full_names,
                           reaction=reaction, reactions=reactions, only_positive=only_positive)

@bp.route('/top_enzymes_and_products/<enzyme_type>/<reaction>', methods=['GET'])
def top_enzymes_and_products(enzyme_type, reaction):

    enzyme_types = activity_queries.enzyme_types_from_activity(reaction=reaction, include_chemical=False)
    enzyme_type_full_names = enzyme_type_queries.get_enzyme_types_with_full_names(enzyme_types)
    reactions = activity_queries.reactions_from_activity(enzyme_type=enzyme_type, include_chemical=False)

    if enzyme_type not in enzyme_types:
        enzyme_type = 'All'
    if reaction not in reactions:
        reaction = 'All'

    remove_negative = None
    only_positive = ''
    args = request.args.to_dict()
    print(args)
    if 'only_positive'in args:
        remove_negative = 'all'
        only_positive = 'checked'

    data_query = get_data.DataQuery(reaction=reaction, enzyme_type=enzyme_type, log_level=1, remove_negative=remove_negative)


    enzyme_count_graph = EnzymeCountGraphCreator().create_graph(data_query)
    enzyme_script, enzyme_div = components(enzyme_count_graph)

    mol_count_graph = MolCountGraphCreator().create_graph(data_query)
    mol_script, mol_div = components(mol_count_graph)


    descriptor_count_graphs = MolDescriptorHistrogramCreator().create_graphs(data_query)
    descriptor_count_scripts_and_divs = []
    for graph in descriptor_count_graphs:
        script, div = components(graph)
        descriptor_count_scripts_and_divs.append([script, div])

    graphs = {'enzyme_graph': {'div': enzyme_div, 'script': enzyme_script},
              'mol_graph': {'div': mol_div, 'script': mol_script},
              'descriptor_graphs': descriptor_count_scripts_and_divs}

    return render_template('summary_page/top_enzymes_and_products.html',
                           graphs=graphs,
                           enzyme_type=enzyme_type, enzyme_types=enzyme_type_full_names,
                           reaction=reaction, reactions=reactions, only_positive=only_positive)
