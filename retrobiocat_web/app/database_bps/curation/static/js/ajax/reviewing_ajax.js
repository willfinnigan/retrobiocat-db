var reviewing_ajax = (function () {
    /*
    Methods for making ajax calls for reviewing parts of the curation process
    */


    function review_paper_metadata(paper_id, reviewed, callback, response_div, on_finish) {
        $.post('/_review_paper_metadata', {
            paper_id: paper_id,
            reviewed: reviewed
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                callback()
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying update paper metadata review status', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function sequences_ready_review(paper_id, ready, callback, response_div, on_finish) {
        // Mark sequences are ready for review

        $.post('/_sequences_ready_for_review', {
            paper_id: paper_id,
            ready: ready
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                callback(data.result.status)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying update paper sequences ready review', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function sequences_review(paper_id, reviewed, callback, response_div, on_finish) {
        $.post('/_sequences_review', {
            paper_id: paper_id,
            reviewed: reviewed
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                callback(data.result.status)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying update paper sequences review status', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

        function activity_ready_review(paper_id, ready, callback, response_div, on_finish) {
        $.post('/_activity_ready_for_review', {
            paper_id: paper_id,
            ready: ready
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                callback(data.result.status)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying update paper activity ready review', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function activity_review(paper_id, reviewed, callback, response_div, on_finish) {
        $.post('/_activity_review', {
            paper_id: paper_id,
            reviewed: reviewed
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                callback(data.result.status)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying update paper activity review status', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function review_all(paper_id, callback, response_div, on_finish) {
        // run all review routes sequentially

        function paper_callback() {
            paper_sequences_reviewed(paper_id, true, seq_callback, response_div, on_finish)
        }

        function seq_callback(status) {
            paper_activity_reviewed(paper_id, true, callback, response_div, on_finish)
        }

         paper_metadata_update(paper_id, true, paper_callback, response_div, on_finish)
    }


    return {
        review_paper_metadata: review_paper_metadata,
        sequences_ready_review: sequences_ready_review,
        sequences_review: sequences_review,
        activity_ready_review: activity_ready_review,
        activity_review: activity_review,
        review_all: review_all
    }

})()