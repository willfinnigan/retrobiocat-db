var upload_molecules_excel_ajax = (function () {


    function set_progress_bar(progress_bar_id, style, innerText, className) {
        var progress_bar = document.getElementById(progress_bar_id)
        progress_bar.style = style
        progress_bar.innerText = innerText
        progress_bar.className = className
    }

    function process_result_msgs(result, response_div) {
        if ("new_mol_msgs" in result) {
            if (result.new_mol_msgs.length !== 0) {
                response_msg_stay(result.msg, 'success', result.new_seq_msgs, response_div)
            }
        }
        if (result.issues.length !== 0) {
            response_msg_stay(result.msg, 'warning', result.issues, response_div)
        }
    }

    function add_data_to_table(result, table) {
        let new_data = result.new_mols_data
        for (const [key, value] of Object.entries(new_data)) {
            let row = table.getRow('' + key)
            if (row === false) {
                table.addRow(value)
            } else {
                table.updateRow(row, value);
            }
        }
    }

    function upload_mol_excel(form_data, mol_table, progressbar_id, response_div, on_finish) {
        set_progress_bar(progressbar_id, "width: 50%", "Uploading...", "progress-bar")
        $.ajax({type: 'POST',
                url: '/_upload_molecules_excel',
                data: form_data,
                contentType: false,
                cache: false,
                processData: false,
            }).done(function(data) {
                process_result_msgs(data.result, response_div)
                if (data.result.status !== 'danger') {
                    add_data_to_table(data.result, mol_table)
                    set_progress_bar(progressbar_id, "width: 100%", "Complete", "progress-bar bg-success")
                } else {
                    response_msg_stay(data.result.msg, data.result.status, data.result.issues, response_div)
                    set_progress_bar(progressbar_id, "width: 50%", "Error", "progress-bar bg-danger")
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error', 'danger', [error], response_div)
                set_progress_bar(progressbar_id, "width: 50%", "Error", "progress-bar bg-danger")
                console.log(error)
            }).always(function(data) {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
    }

    return {
        upload_mol_excel: upload_mol_excel
      }

})()


