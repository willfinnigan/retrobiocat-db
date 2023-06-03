function add_info_cluster_dict(cluster_dict, node_dict) {
    // this function takes a dict for a single node and cumulatively adds its info to cluster dict

    if (node_dict['node_type'] !== 'biocatdb') {
        // node vars
        var cluster_num = node_dict['cluster_label']
        var name = node_dict['protein_name']
        var rheas = node_dict['rhea']
        var pfams = node_dict['pfam_codes']

        // add cluster labels or increase count
        if (!cluster_dict['cluster_labels'].hasOwnProperty(cluster_num)) {
            cluster_dict['cluster_labels'][cluster_num] = 1
        } else {
            cluster_dict['cluster_labels'][cluster_num] += 1
        }

        // add name or increase count
        if (!cluster_dict['names'].hasOwnProperty(name)) {
            cluster_dict['names'][name] = 1
        } else {
            cluster_dict['names'][name] += 1
        }

        // add rheas
        rheas.forEach(function (rhea, index) {
            if (!cluster_dict['rheas'].hasOwnProperty(rhea)) {
                cluster_dict['rheas'][rhea] = 1
            } else {
                cluster_dict['rheas'][rhea] += 1
            }
        });

        // add pfams
        pfams.forEach(function (pfam, index) {
            if (!cluster_dict['pfams'].hasOwnProperty(pfam)) {
                cluster_dict['pfams'][pfam] = 1
            } else {
                cluster_dict['pfams'][pfam] += 1
            }
        });
    } else {
        cluster_dict['biocatdb'].push(node_dict['id'])
    }

    // add biocatdb seqs

    return cluster_dict
}

function get_selected_nodes_info(selected_nodes, data) {
    var node_data = data['nodes'].get(selected_nodes)

    var cluster_dict = {'names': {},
                        'rheas': {},
                        'pfams': {},
                        'cluster_labels': {},
                        'biocatdb': []
                         }

    // loop over nodes, adding info to dict
    node_data.forEach(function (node_dict, index) {
        cluster_dict = add_info_cluster_dict(cluster_dict, node_dict)
    })

    return cluster_dict
}

function get_sorted_array(dict) {
    keysSorted = Object.keys(dict).sort(function(a,b){return dict[b]-dict[a]})
    return keysSorted
}

function add_names_to_div(div_id, names_dict) {
    var name_div = document.getElementById(div_id);
    name_div.innerHTML = ""

    var sorted_keys = get_sorted_array(names_dict)

    sorted_keys.forEach(function (name, index) {
        var new_name = document.createElement("p")
        new_name.innerHTML = "<b>" + names_dict[name] + "</b> " + name
        name_div.appendChild(new_name);
    });
}

function add_biocatdb_to_div(div_id, biocatdb_list) {
    var name_div = document.getElementById(div_id);
    name_div.innerHTML = ""

    if (biocatdb_list.length === 0) {
        name_div.innerHTML = "<p><i>None</i></p>"
    } else {
        var new_seq_ele = document.createElement("p")
        name_div.appendChild(new_seq_ele);
        biocatdb_list.forEach(function (seq, index) {
            new_seq_ele.innerHTML += seq
            if (index !== biocatdb_list.length-1) {
                new_seq_ele.innerHTML += ", "
            }
        });
    }

}

function add_pfam_to_div(div_id, pfam_dict, pfam_id_dict) {
    var name_div = document.getElementById(div_id);
    name_div.innerHTML = ""

    var sorted_keys = get_sorted_array(pfam_dict)

    sorted_keys.forEach(function (pfam, index) {
        var new_pfam = document.createElement("p")

        var pfam_count = document.createElement('b')
        pfam_count.innerHTML = pfam_dict[pfam] + ' '
        new_pfam.append(pfam_count)

        var pfam_link = document.createElement('a')
        pfam_link.innerHTML = pfam + " - " + pfam_id_dict[pfam]
        pfam_link.href = 'https://pfam.xfam.org/family/' + pfam
        pfam_link.target = '_blank'
        new_pfam.append(pfam_link)

        name_div.appendChild(new_pfam);
    });
}

function add_rhea_to_div(div_id, rhea_dict, rhea_img_dict, rhea_equation_dict) {
    var rhea_div = document.getElementById(div_id);
    rhea_div.innerHTML = ""

    var sorted_keys = get_sorted_array(rhea_dict)

    parse_rheas(rhea_div, sorted_keys, rhea_img_dict, rhea_equation_dict, only_draw_first=true)

}

function create_list_clusters(nodes) {
    var cluster_list = []
    nodes.forEach(function (node_dict, index) {
        if (!cluster_list.includes(node_dict['cluster_label'])) {
            cluster_list.push(node_dict['cluster_label'])
        }
    })

    cluster_list.sort(function(a, b){return a-b})
    return cluster_list
}

