var paper_assignment_ajax = (function () {
    /*
    Methods for making ajax calls when changing paper assignment
    All functions have a response div to show if operation successful or not.
    All functions also have an optional on_finish argument for a function to run on completion
    */

    function self_assign_paper(paper_id, callback, response_div, on_finish) {
        $.post('/_self_assign', {
            paper_id: paper_id
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback()
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change paper ownership', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function unassign_paper(paper_id, callback, response_div, on_finish) {
        $.post('/_un_self_assign', {
            paper_id: paper_id
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                callback(data.result.status)

            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change paper ownership', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    return {
        unassign_paper: unassign_paper,
        self_assign_paper: self_assign_paper,
    }

})()