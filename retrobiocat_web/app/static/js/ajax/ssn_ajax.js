var ssn_ajax = (function () {
    /*
    Methods for making ajax calls for working with ssns
    */

    function download_selected_unireps(unirep_ids, callback, on_finish) {
        // for a set of unirep ids, get the sequences

        $.post('/_download_selected_unireps', {
            unirep_ids: JSON.stringify(unirep_ids)
        }).done(function (data) {
            callback(data.result)
        }).fail(function (xhr, status, error) {
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    return {
        download_selected_unireps: download_selected_unireps,
    }

})()