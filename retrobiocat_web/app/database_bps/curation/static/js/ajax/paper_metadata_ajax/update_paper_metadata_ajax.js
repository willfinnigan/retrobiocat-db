var update_paper_metadata_ajax = (function () {
    /*
    Methods for making ajax calls to update paper metadata fields
    */


    function save_updated_paper_metadata(paper_id, update_dict, callback, response_div, on_finish) {
        $.post('/_save_updated_paper_metadata', {
            paper_id: paper_id,
            update_dict: JSON.stringify(update_dict)
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback()
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to update paper metadata', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function update_paper_importance(paper_id, important_bool, response_div) {
        $.post('/_update_paper_importance', {
            paper_id: paper_id,
            importance: important_bool
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to update paper metadata', 'danger', [error], response_div)
                console.log(error)
            })
        }


    return {
        save_updated_paper_metadata: save_updated_paper_metadata,
        update_paper_importance: update_paper_importance
    }

})()