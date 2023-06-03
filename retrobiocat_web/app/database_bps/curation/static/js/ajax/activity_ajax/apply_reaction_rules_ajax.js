var apply_reaction_rules_ajax = (function () {


    function apply_reaction_rules(selectedData, selectedRows, response_div, callback, on_finish) {
        var interconvert = false
        $.post('/_apply_reaction_rules_table', {
            rows: JSON.stringify(selectedData)
        }).done(function (data) {
            if (data.result.status === 'success') {
                callback(selectedRows, data.result, interconvert)
            }
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
        }). fail(function (xhr, status, error) {
            response_msg('JS error while trying set issues status', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function autodetect_reaction_names(selectedData, selectedRows, response_div, callback, on_finish) {
        var interconvert = false
        $.post('/_autodetect_reaction_names', {
            rows: JSON.stringify(selectedData)
        }).done(function (data) {
            if (data.result.status === 'success') {
                callback(selectedRows, data.result, interconvert)
            }
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
        }). fail(function (xhr, status, error) {
            response_msg('JS error while trying set issues status', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }






return {
    apply_reaction_rules: apply_reaction_rules,
    autodetect_reaction_names: autodetect_reaction_names
}

})
()