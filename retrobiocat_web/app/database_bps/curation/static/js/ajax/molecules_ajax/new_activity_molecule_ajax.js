var new_activity_molecule_ajax = (function () {

    function activity_mol_from_name(paper_id, mol_name, callback, response_div, on_finish) {
        $.post('/_activity_mol_from_name', {
            paper_id: paper_id,
            mol_name: mol_name
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying set issues status', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

   function new_activity_mol(paper_id, smi, callback, response_div, on_finish) {
        $.post('/_new_activity_mol', {
            paper_id: paper_id,
            smi: smi
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying set issues status', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }


    return {
        activity_mol_from_name: activity_mol_from_name,
        new_activity_mol: new_activity_mol,
    }

})()