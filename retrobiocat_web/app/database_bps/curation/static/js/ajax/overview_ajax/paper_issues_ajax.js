var paper_issues_ajax = (function () {
    /*
    Methods for making ajax calls when changing paper assignment
    All functions have a response div to show if operation successful or not.
    All functions also have an optional on_finish argument for a function to run on completion
    */

    function set_paper_has_issues(paper_id, issues_bool, response_div) {
        $.post('/_paper_issues', {
            paper_id: paper_id,
            issues: issues_bool
            }).done(function(data) {
                response_msg(data.result.msg, data.result.status, data.result.issues, response_div)
            }).fail(function (xhr, status, error) {
                response_msg('JS error while trying set issues status', 'danger', [error], response_div)
                console.log(error)
            })
        }



    return {
        set_paper_has_issues: set_paper_has_issues,
    }

})()