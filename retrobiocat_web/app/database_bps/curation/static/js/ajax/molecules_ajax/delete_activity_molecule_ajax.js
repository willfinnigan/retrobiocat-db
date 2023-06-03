var delete_activity_molecule_ajax = (function () {

    function delete_activity_molecule(paper_id, mol_id, callback, response_div) {
        $.post('/_delete_activity_molecule', {
            paper_id: paper_id,
            mol_id: mol_id
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to delete molecule ' + mol_id, 'danger', [error], response_div)
                console.log(error)
            })
        }

   function delete_many_paper_molecules(paper_id, mode, callback, response_div) {
        $.post('/_delete_many_paper_molecules', {
            paper_id: paper_id,
            mode: mode
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to delete multiple molecules', 'danger', [error], response_div)
                console.log(error)
            })
        }


    return {
        delete_activity_molecule: delete_activity_molecule,
        delete_many_paper_molecules: delete_many_paper_molecules
    }

})()