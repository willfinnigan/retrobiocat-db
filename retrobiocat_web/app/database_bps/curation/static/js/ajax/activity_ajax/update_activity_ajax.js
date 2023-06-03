var update_activity_ajax = (function () {

    function save_activity_data(data, paper_id, edited_rows, response_div, on_finish) {
        $.post('/_save_activity_data', {
            paper_id: paper_id,
            data: JSON.stringify(data),
            edited_rows: JSON.stringify(edited_rows)
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    window.onbeforeunload = null;  // Remove navigation prompt
                    location.reload();
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
        save_activity_data: save_activity_data,
    }



})()