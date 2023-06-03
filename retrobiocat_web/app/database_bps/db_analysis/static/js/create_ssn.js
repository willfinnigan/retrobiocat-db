function create_network_data(node_list, edge_list) {
    const nodes = new vis.DataSet(node_list)
    const edges = new vis.DataSet(edge_list)

    return {'nodes': nodes, 'edges': edges};
}

function setup_network_options(options_div_id) {
    var edges_options = {'smooth': false}

    var physics_options = {'enabled': false,
                            'stabilization': {'enabled': false,
                                             'iterations': 50},
                            "repulsion": {
                                  "centralGravity": 0.1,
                                  "springLength": 860,
                                  "nodeDistance": 860,
                                  "damping": 0.5
                            },
                            "maxVelocity": 100,
                            "minVelocity": 0.75,
                            "solver": "repulsion",
                            "timestep": 0.75}

    var interaction_options = {'tooltipDelay': 0,
                               'hover': false,
                               'hoverConnectedEdges': false,
                               'selectConnectedEdges': false,
                               'zoomSpeed': 0.5}

    var configure_options = {filter: function (option, path) {
                                        if (path.indexOf("physics") !== -1) {
                                            return true;}
                                        if (path.indexOf("smooth") !== -1 || option === "smooth") {
                                            return true;
                                        }
                                        return false;},
                                        container: document.getElementById(options_div_id)}

    var options = {edges: edges_options,
                   physics: physics_options,
                   interaction: interaction_options,
                   configure: configure_options,
                   nodes: {shapeProperties: {interpolation: false}}
                  };

    return options

}

function create_ssn(network_data, options, ssn_div_id) {
    var container = document.getElementById(ssn_div_id);

    return new vis.Network(container, network_data, options);

}



