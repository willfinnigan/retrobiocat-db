var activity_excel_ajax = (function () {

    function upload_activity_excel(form_data, table, response_div, callback, on_fail) {

        $.ajax({
                type: 'POST',
                url: '/_upload_activity_excel',
                data: form_data,
                contentType: false,
                cache: false,
                processData: false,
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                } else {
                    onfail()
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to upload excel', 'danger', [error], response_div)
                console.log(error)
                on_fail()
            })
        }

    return {
        upload_activity_excel: upload_activity_excel,
    }

})()
