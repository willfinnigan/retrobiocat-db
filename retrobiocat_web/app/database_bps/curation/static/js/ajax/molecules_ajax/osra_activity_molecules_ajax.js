var osra_activity_molecule_ajax = (function () {

    function add_smi_to_paper(paper_id, smi, name, response_div, callback) {
        $.post('/_add_smi_to_paper', {
            paper_id: paper_id,
            smi: smi,
            name: name
            }).done(function(data) {
                if (data.result.status === 'success') {
                    callback()
                }
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to delete molecule ' + smi, 'danger', [error], response_div)
            })
        }

    function remove_osr_molecule(paper_id, smi, response_div, callback) {
        $.post('/_remove_osr_molecule', {
            paper_id: paper_id,
            smi: smi
            }).done(function(data) {
                if (data.result.status === 'success') {
                    callback()
                }
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying to delete molecule ' + smi, 'danger', [error], response_div)
            })
        }

    function submit_images_for_processing(form_data, callback, response_div, on_finish) {
        $.ajax({
            type: 'POST',
            url: '/_upload_molecule_images',
            data: form_data,
            dataType: 'json',
            contentType: false,
            cache: false,
            processData: false,
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)  //"mol_upload_response"
            if (data.result.status === 'success') {
                callback(data.result)  //getProcessingStatus();
            }
        }).fail(function (xhr, status, error) {
            response_msg('JS error while trying to send images for osra', 'danger', [error], response_div)
            console.log(error)
        }).always(function (data) {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        });
    }

    return {
        add_smi_to_paper: add_smi_to_paper,
        remove_osr_molecule: remove_osr_molecule,
        submit_images_for_processing: submit_images_for_processing,
        }
})()



