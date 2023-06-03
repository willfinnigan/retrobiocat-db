function move_view(selected) {
    if (selected.length !== 0) {
        var node = data['nodes'].get(selected)[0]
        var pos = {'x': node['x'], 'y': node['y']}
        network.moveTo({position: pos,
                        scale: 0.5})
        network.fit({nodes: selected})
    }
}

function move_view_offset(selected) {
    if (selected.length !== 0) {
        network.fit({nodes: selected})
        var pos = network.getViewPosition()
        pos['y'] -= 250
        network.moveTo({position: pos,
                        scale: 0.7})
    }
}

function select_cluster() {
    var selected = data['nodes'].get(network.getSelectedNodes())
    var selected_clusters = []
    selected.forEach(function (node_dict, i) {
        if (!(selected_clusters.includes(node_dict['cluster_label']))) {
            selected_clusters.push(node_dict['cluster_label'])
        }
    })

    var to_select = []
    data['nodes'].forEach(function (node_dict, i) {
        if (selected_clusters.includes(node_dict['cluster_label'])) {
            to_select.push(node_dict['id'])
        }
    })
    network.selectNodes(to_select)
}

function select_cluster_by_number(number) {
    var to_select = []
    data['nodes'].forEach(function (node_dict, i) {
        if (number === node_dict['cluster_label']) {
            to_select.push(node_dict['id'])
        }
    })
    network.selectNodes(to_select)
    return to_select
}

function stabilise_network() {
            response_msg('Stabilising please wait...', 'success', [], "live_ssn_response_bar")
            network.stabilize(500)
        }