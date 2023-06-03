var update_activity_molecule_ajax = (function () {

    function update_activity_molecule(paper_id, mol_id, smi, name, callback, response_div, on_finish) {
        $.post('/_update_activity_molecule', {
            paper_id: paper_id,
            mol_id: mol_id,
            mol_name: name,
            smi: smi
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to delete molecule ' + mol_id, 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function update_activity_molecule_name_only(paper_id, mol_id, name, callback, response_div, on_finish) {
        $.post('/_update_activity_molecule_name_only', {
            paper_id: paper_id,
            mol_id: mol_id,
            mol_name: name,
            }).done(function(data) {
                callback(data.result)
                if (data.result.status !== 'success') {
                    response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to delete molecule ' + mol_id, 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }


    return {
        update_activity_molecule: update_activity_molecule,
        update_activity_molecule_name_only: update_activity_molecule_name_only
    }

})()