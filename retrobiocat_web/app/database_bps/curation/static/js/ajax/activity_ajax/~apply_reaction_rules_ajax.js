var apply_reaction_rules_ajax = (function () {


    function apply_reaction_rules(selectedData) {
        return new Promise((resolve, reject) => {
            $.post('/_apply_reaction_rules_table', {
                rows: JSON.stringify(selectedData)
            }).done(function (data) {
                resolve(data)
            }).fail(function (xhr, status, error) {
                reject(error)
            })
        })
    }

    function autodetect_reaction_names(selectedData) {
        return new Promise((resolve, reject) => {
            $.post('/_autodetect_reaction_names', {
                rows: JSON.stringify(selectedData)
            }).done(function (data) {
                resolve(data)
            }). fail(function (xhr, status, error) {
                reject(error)
            })
        })
    }



return {
    apply_reaction_rules: apply_reaction_rules,
    autodetect_reaction_names: autodetect_reaction_names
}

})
()