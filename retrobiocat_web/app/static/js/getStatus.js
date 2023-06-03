

function general_getStatus(status_url, taskID, progress_bar_obj, progressDict, timeout, completion_func) {
    $.get(status_url + '/' + taskID, {}).done(function(response) {
        const taskStatus = response.data.task_status;
        const taskProgress = response.data.task_progress;
        progress_bar_obj.style = progressDict[taskProgress][0]
        progress_bar_obj.innerText = progressDict[taskProgress][1]
        progress_bar_obj.className = progressDict[taskProgress][2]

        if (taskStatus === 'finished') {
            progress_bar_obj.style = "width: 100%"
            progress_bar_obj.innerText = "Done"
            completion_func(task_id)
            return false;
        } else if (taskStatus === 'failed') {
            progress_bar_obj.style = "width: 50%"
            progress_bar_obj.innerText = "Failed"
            progress_bar_obj.className = "progress-bar bg-danger"
            return false;
        } else {
            setTimeout(function() {
                general_getStatus(status_url, taskID, progress_bar_obj, progressDict, timeout, completion_func);
                }, timeout);
        }
    })
}

export {general_getStatus}
