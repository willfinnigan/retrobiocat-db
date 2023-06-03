var curation_admin_ajax = (function () {
    /*
    Methods for making ajax calls when using admin tab for changing paper assignments
    All functions have a response div to show if operation successful or not.
    All functions also have an optional on_finish argument for a function to run on completion
    */


    function admin_set_owner(paper_id, new_owner_id, callback, response_div, on_finish) {
        $.post('/_admin_set_owner', {
            paper_id: paper_id,
            new_owner_id: new_owner_id
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback()
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change set ownership', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function admin_activity_to_owner(paper_id, response_div) {
        $.post('/_admin_activity_to_owner', {
            paper_id: paper_id,
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change set activity ownership', 'danger', [error], response_div)
                console.log(error)
            })
        }

    function admin_unassigned_seqs_to_owner(paper_id, response_div) {
        $.post('/_admin_unassigned_seqs_to_owner', {
            paper_id: paper_id,
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change set sequence ownership', 'danger', [error], response_div)
                console.log(error)
            })
        }


    function admin_all_seqs_to_owner(paper_id, response_div) {
        $.post('/_admin_all_seqs_to_owner', {
            paper_id: paper_id,
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change set sequence ownership', 'danger', [error], response_div)
                console.log(error)
            })
        }


    return {
        admin_set_owner: admin_set_owner,
        admin_activity_to_owner: admin_activity_to_owner,
        admin_unassigned_seqs_to_owner: admin_unassigned_seqs_to_owner,
        admin_all_seqs_to_owner: admin_all_seqs_to_owner,
    }

})()