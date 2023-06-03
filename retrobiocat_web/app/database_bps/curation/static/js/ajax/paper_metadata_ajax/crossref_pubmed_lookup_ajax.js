var crossref_pubmed_lookup_ajax = (function () {
    /*
    Methods for making ajax calls to get paper metadata
    */

    function query_pubmed_or_crossref(paper_id, service_name, callback, response_div, on_finish) {
        $.post('/_query_pubmed_or_crossref', {
            paper_id: paper_id,
            service_name: service_name
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying get metadata from pubmed/crossref', 'danger', [error], response_div)
                console.log(error)
            })
        }

    return {
        query_pubmed_or_crossref: query_pubmed_or_crossref,
    }

})()