function response_msg(msg, type, issues, parent_div) {
    if (typeof parent_div === 'undefined') {
        return
    }


    var newDiv = document.createElement('div');
    newDiv.className = "alert alert-" + type;
    newDiv.setAttribute("role", "alert");
    newDiv.innerHTML = msg
    var id = 'ID_' + Math.random().toString(36).substr(2, 9);
    newDiv.setAttribute("id", id);

    var close_button = document.createElement('button');
    close_button.className = "close"
    close_button.setAttribute("type", "button");
    close_button.setAttribute("data-dismiss", "alert");
    close_button.innerHTML = "<span aria-hidden=\"true\">&times;</span>"
    newDiv.appendChild(close_button)

    document.getElementById(parent_div).appendChild(newDiv)

    if (issues.length !== 0) {
        newList = document.createElement('ul');
        issues.forEach(function (item, index) {
            var newListItem = document.createElement('li');
            newListItem.innerHTML = item
            newList.appendChild(newListItem)
        });
        newDiv.appendChild(newList)
    }

    var timeout_time = 20000 + (issues.length*5000)

    setTimeout(function() {
        $("#" + id).alert('close');
    }, timeout_time);
}

function close_alerts() {
    $(".alert").alert('close')
}

function response_msg_stay(msg, type, issues, parent_div) {
    var newDiv = document.createElement('div');
    newDiv.className = "alert alert-" + type;
    newDiv.setAttribute("role", "alert");
    newDiv.innerHTML = msg
    var id = 'ID_' + Math.random().toString(36).substr(2, 9);
    newDiv.setAttribute("id", id);

    var close_button = document.createElement('button');
    close_button.className = "close"
    close_button.setAttribute("type", "button");
    close_button.setAttribute("data-dismiss", "alert");
    close_button.innerHTML = "<span aria-hidden=\"true\">&times;</span>"
    newDiv.appendChild(close_button)

    document.getElementById(parent_div).appendChild(newDiv)

    if (issues.length !== 0) {
        newList = document.createElement('ul');
        issues.forEach(function (item, index) {
            var newListItem = document.createElement('li');
            newListItem.innerHTML = item
            newList.appendChild(newListItem)
        });
        newDiv.appendChild(newList)
    }
}

