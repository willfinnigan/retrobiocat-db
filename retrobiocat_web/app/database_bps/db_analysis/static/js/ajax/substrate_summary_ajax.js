var substrate_summary_ajax = (function () {
    /*
    Methods for making ajax calls for the substrate summary page
    */

    function get_examples(leafs, args, callback, response_div, on_finish) {
        // get the examples for a node

        $.post('/_load_substrate_summary_examples', {
            leafs: JSON.stringify(leafs),
            args: JSON.stringify(args),
        }).done(function (data) {
            callback(data.html)
        }).fail(function (xhr, status, error) {
            response_msg('Error loading examples', 'danger', [''+error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function send_heatmap_mol_list(leafs, callback, response_div, on_finish) {
        // send the list of smis to use in the heatmap

        $.post('/_send_heatmap_mol_list', {
            leafs: JSON.stringify(leafs),
        }).done(function (data) {
            callback(data.mol_list_uuid)
        }).fail(function (xhr, status, error) {
            response_msg('Error loading examples', 'danger', [''+error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    return {
        get_examples: get_examples,
        send_heatmap_mol_list: send_heatmap_mol_list
    }

})()