var get_status_ajax = (function () {

    function progress_bar_complete(progressbar_id) {
        document.getElementById(progressbar_id).style =  "width: 100%"  //"width: 50%"
        document.getElementById(progressbar_id).innerText = 'Complete' // "Loading heatmap"
        document.getElementById(progressbar_id).className = "progress-bar bg-success"  //"progress-bar"
    }

    function progress_bar_failed(progressbar_id) {
        document.getElementById(progressbar_id).style =  "width: 75%"  //"width: 50%"
        document.getElementById(progressbar_id).innerText = 'Failed' // "Loading heatmap"
        document.getElementById(progressbar_id).className = "progress-bar bg-danger"  //"progress-bar"
    }

    function progress_bar_from_response(progressbar_id, response) {
        document.getElementById(progressbar_id).style =  response.data.task_progress.progressbar_style  //"width: 50%"
        document.getElementById(progressbar_id).innerText = response.data.task_progress.progressbar_text // "Loading heatmap"
        document.getElementById(progressbar_id).className = response.data.task_progress.progressbar_classname  //"progress-bar"
    }

    function get_status(taskID, taskQueue, progressbar_id, retry_time, callback) {
        var url = "/queue_status/"+taskID+"/"+taskQueue
        $.get(url, {
        }).done(function (response) {
            progress_bar_from_response(progressbar_id, response)
            const taskStatus = response.data.task_status;
            if (taskStatus === 'finished') {
                progress_bar_complete(progressbar_id)
                callback(response.data)
                return false;
            } else if (taskStatus === 'failed') {
                progress_bar_failed(progressbar_id)
                return false;
            } else {
                setTimeout(function() {
                    get_status(taskID, taskQueue, progressbar_id, retry_time, callback);
                    }, retry_time);
            }
        }).fail(function (xhr, status, error) {
            console.log(error)
        })
    }

    return {
        get_status: get_status,
    }

})()