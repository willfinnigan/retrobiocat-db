var seq_load_ajax = (function () {
    /*
    Methods for making ajax calls to load sequence information
    */

    function get_seq_papers(enzyme_name, callback) {
        // get the papers which a sequence is attached to

        $.post('/_load_sequence_papers', {
            enzyme_name: enzyme_name,
            }).done(function(data) {
                callback(data.result.papers) //a dict containing paper rows
            })
    }

    function get_sequences_of_same_type(enzyme_name, callback) {
        // get the sequences which are the same type

        $.post('/_sequences_of_same_type', {
            enzyme_name: enzyme_name,
            }).done(function(data) {
                callback(data.result.sequences) // a dict where both key and value are seq name
            })
        }

    function get_sequences_of_type(enzyme_type, callback) {
        // get the sequences of a type

        $.post('/_sequences_of_type', {
            enzyme_type: enzyme_type,
            }).done(function(data) {
                callback(data.result.sequences) // a dict where both key and value are seq name
            })
        }

    function get_sequences_of_type_with_other_names(enzyme_type, callback) {
        // get the sequences of a type and the other_names

        $.post('/_sequences_of_type_with_other_names', {
            enzyme_type: enzyme_type,
            }).done(function(data) {
                callback(data.result.sequences) // a dict where both key and value are seq name
            })
        }


    function load_seq_data(enzyme_name, callback) {
        // load all the data for a given sequence

        $.post('/_load_sequence_data', {
            enzyme_name: enzyme_name,
            }).done(function(data) {
                callback(data.result)
            })
        }

    function load_seq_other_names_data(enzyme_name, existing_name, callback) {
        // load the other_names data for an enzyme

        $.post('/_load_sequence_other_names_data', {
            enzyme_name: enzyme_name,
            existing_name: existing_name
            }).done(function(data) {
                callback(data.result)
            })
        }

    function get_possible_alternative_naming_for_paper(paper_id, callback) {
        // get all the other names for the sequences in a paper

        $.post('/_get_possible_alternative_naming_for_paper', {
            paper_id: paper_id,
            }).done(function(data) {
                callback(data.result.alt_names)
            })
        }

    return {
        get_seq_papers: get_seq_papers,
        get_sequences_of_same_type: get_sequences_of_same_type,
        get_sequences_of_type: get_sequences_of_type,
        get_sequences_of_type_with_other_names: get_sequences_of_type_with_other_names,
        load_seq_data: load_seq_data,
        load_seq_other_names_data: load_seq_other_names_data,
        get_possible_alternative_naming_for_paper: get_possible_alternative_naming_for_paper
    }

})()