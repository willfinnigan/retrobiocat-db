var self_assign_ajax = (function () {

    function self_assign(paper_id, response_div) {
        $.post('/_self_assign', {
            paper_id: paper_id
        }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            if (data.result.status === 'success') {
                window.location.href = Flask.url_for("curation.curation_overview", {"paper_id": paper_id})
            }
        })
    }

    return {
        self_assign: self_assign,
    }

})()

