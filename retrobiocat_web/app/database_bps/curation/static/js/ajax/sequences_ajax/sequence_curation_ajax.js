var seq_curate_ajax = (function () {
    /*
    Methods for making ajax calls when curating sequence information.
    All functions have a response div to show if operation successful or not.
    All functions also have an optional on_finish argument for a function to run on completion
    */

    function uniprot_lookup(accession, callback, response_div, on_finish) {
        // Lookup a uniprot or ncbi accession number and result the sequence

        $.post('_get_sequence_from_uniprot', {
            accession: accession
        }).done(function (data) {
            if (data.result.status === 'success' || data.result.status === 'warning') {
                callback(data.result.seq)
            }
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
        }).fail(function (xhr, status, error) {
            response_msg('JS error', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function delete_seq(seq_name, callback, response_div, on_finish) {
        // Try to delete a sequence

        $.post('/_delete_sequence', {
            to_delete: seq_name
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success' || data.result.status === 'warning') {
                    callback(seq_name)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while deleting seq', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }


    function merge_seqs(merge_with, to_merge, save_other_name, callback, response_div, on_finish) {
        $.post('/_merge_seq', {
            merge_with: merge_with,
            to_merge: to_merge,
            save_other_name: save_other_name
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success' || data.result.status === 'warning') {
                    callback(merge_with, to_merge)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while merging seqs', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function save_seq_edit(update_dict, callback, response_div, on_finish) {
        // save an update to an existing sequence entry

        $.post('/_save_edited_sequence', update_dict).done(function(data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            if (data.result.status === 'success' || data.result.status === 'warning') {
                callback(data.result) //for updating tabulator
            }
        }).fail(function (xhr, status, error) {
            response_msg('JS error while saving seq', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
    }


    function add_new_sequence(new_enzyme_name, enzyme_type, paper_id, callback, response_div, on_finish) {
        $.post('/_add_new_sequence', {
            enzyme_type: document.getElementById("new_enzyme_type").value,
            new_name: document.getElementById("new_enzyme_name").value,
            paper_id: paper_id

            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success' || data.result.status === 'warning') {
                    callback(data.result.seq_table_entry) //for updating tabulator)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while creating new seq', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function add_existing_sequence(existing_name, paper_id, callback, response_div, on_finish) {
        $.post('/_add_existing_sequence', {
            existing_name: existing_name,
            paper_id: paper_id

            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success' || data.result.status === 'warning') {
                    callback(data.result.seq_table_entry) //for updating tabulator)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying add existing sequence', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function change_sequence_assign(enzyme_name, self_assigned, callback, response_div, on_finish) {
        $.post('/_change_sequence_assign', {
            enzyme_name: enzyme_name,
            self_assigned: self_assigned,
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            callback(enzyme_name) // reload_enzyme_data to update any disabled fields
        }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change sequence ownership', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
    }

    function change_sequence_reviewed(enzyme_name, reviewed, callback, response_div, on_finish) {
        $.post('/_mark_sequence_reviewed', {
            enzyme_name: enzyme_name,
            reviewed: reviewed,
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            callback(enzyme_name) // reload_enzyme_data to update any disabled fields
        }).fail(function (xhr, status, error) {
                response_msg('JS error while trying change sequence ownership', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
    }

    function update_other_names(enzyme_name, other_names_dict, callback, response_div, on_finish) {
        $.post('/_update_other_names', {
            enzyme_name: enzyme_name,
            other_names_dict: JSON.stringify(other_names_dict),
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            callback(data.result)
        }).fail(function (xhr, status, error) {
            response_msg('JS error while trying update other names', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function delete_other_names(enzyme_name, other_name, callback, response_div, on_finish) {
        $.post('/_delete_other_names', {
            enzyme_name: enzyme_name,
            other_name: other_name,
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            callback(data.result)
        }).fail(function (xhr, status, error) {
            response_msg('JS error while trying delete other name', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function save_alt_naming_selection(alt_names_to_set, paper_id, callback, response_div, on_finish) {
        $.post('/_save_alt_naming_selection', {
            alt_names_to_set: JSON.stringify(alt_names_to_set),
            paper_id: paper_id,
        }).done(function (data) {
            response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            callback(data.result)
        }).fail(function (xhr, status, error) {
            response_msg('JS error while trying delete other name', 'danger', [error], response_div)
            console.log(error)
        }).always(function () {
            if (typeof on_finish !== 'undefined') {
                on_finish()
            }
        })
    }

    function remove_sequence_from_paper(name, paper_id, callback, response_div, on_finish) {
        $.post('/_remove_seq_from_paper', {
            paper_id: paper_id,
            enzyme_name: name,

            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
                if (data.result.status === 'success') {
                    callback(data.result)
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying delete other name', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    function remove_all_sequence_from_paper(paper_id, callback, response_div, on_finish) {
        $.post('/_remove_all_seq_from_paper', {
            paper_id: paper_id,

            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, "seq_response")
                if (data.result.status === 'success') {
                    callback()
                }
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying delete other name', 'danger', [error], response_div)
                console.log(error)
            }).always(function () {
                if (typeof on_finish !== 'undefined') {
                    on_finish()
                }
            })
        }

    return {
        uniprot_lookup: uniprot_lookup,
        delete_seq: delete_seq,
        merge_seqs: merge_seqs,
        save_seq_edit: save_seq_edit,
        add_existing_sequence: add_existing_sequence,
        add_new_sequence: add_new_sequence,
        change_sequence_assign: change_sequence_assign,
        change_sequence_reviewed: change_sequence_reviewed,
        update_other_names: update_other_names,
        delete_other_names: delete_other_names,
        save_alt_naming_selection: save_alt_naming_selection,
        remove_sequence_from_paper: remove_sequence_from_paper,
        remove_all_sequence_from_paper: remove_all_sequence_from_paper
      }

})()


