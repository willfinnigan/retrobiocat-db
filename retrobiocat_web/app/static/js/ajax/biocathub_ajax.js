var biocathub_ajax = (function () {
    /*
    Methods for making ajax calls for working with biocathub
    */

    function get_network_json(network_id, callback, response_div, on_finish) {
        // get the json for a network, ready to send to biocathub

        $.post('/_get_biocathub_network_json', {
            network_id: network_id
        }).done(function (data) {
            callback(data.result.biocathub_json)
        }).fail(function (xhr, status, error) {
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function send_and_redirect(biocathub_json, callback, response_div, on_finish) {
        // send json to biocathub and redirect if ok

        $.post('/_send_data_to_retrobiohub', {
            biocathub_json: biocathub_json
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            if (data.result.status_code === 200) {
                callback(data.result.redirect_url)
            }
        }).fail(function (xhr, status, error) {
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    return {
        get_network_json: get_network_json,
        send_and_redirect: send_and_redirect
    }

})()