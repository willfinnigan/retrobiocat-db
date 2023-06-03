function set_select(select_id, arraySeq) {
    var select = document.getElementById(select_id);
    select.options.length = 1;
    for (index in arraySeq) {
        select.options[select.options.length] = new Option(arraySeq[index], index);
    }
}

