
    function addNodes(data, newNodes, newEdges) {
        for (let i = 0; i < newNodes.length; i++) {
            data['nodes'].update(newNodes[i])
        }
        for (let i = 0; i < newEdges.length; i++) {
            data['edges'].update(newEdges[i])
        }
    }

    function remove_all_nodes(data) {
        data['nodes'].clear()
        data['edges'].clear()
    }

    function removeNodes(data, nodesToRemove) {

        data['nodes'].remove(nodesToRemove)

        var edges = data['edges'].getIds({
          filter: function(item) {
            return (nodesToRemove.indexOf(item.to)   !== -1) ||
                   (nodesToRemove.indexOf(item.from) !== -1);
            }
          });

        data['edges'].remove(edges)
    }

    function get_biocathub_reaction_json(reaction, network) {
        var reaction_id = reaction['id']
        var products = network.getConnectedNodes(reaction_id, 'to')
        var substrates = network.getConnectedNodes(reaction_id, 'from')
        var enzyme = reaction['metadata']['selected_enzyme']
        var cofactors = reaction['metadata']['enzyme_cofactors'][enzyme]
        var reaction_name = reaction['label']

        return {'substrates': substrates,
                'products': products,
                'cofactors': cofactors,
                'enzyme': enzyme,
                'reaction': reaction_name}




    }
