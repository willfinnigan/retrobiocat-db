function parse_rheas(dom_div, rheas, rhea_img_dict, rhea_equation_dict, only_draw_first=true) {
    dom_div.innerHTML = ""
    var drawn_first = false
    rheas.forEach(function (item, index) {
        var rhea_code = item.replaceAll('RHEA:', '')
        var rhea_title = document.createElement("div");
        dom_div.appendChild(rhea_title);
        rhea_title.innerHTML = "<a href='https://www.rhea-db.org/rhea/" + rhea_code + "' target='_blank'>" + item + "</a>"

        var rhea_draw = document.createElement("div");
        dom_div.appendChild(rhea_draw);
        rhea_draw.className = "row my-auto text-center"
        rhea_draw.innerHTML = '<p class="ml-5">Loading rhea reaction..</p>'

        if (drawn_first === false || only_draw_first === false) {
            //rhea_ajax(item, rhea_draw)
            draw_rhea(item, rhea_img_dict, rhea_draw, 125, 90)
            drawn_first = true
        } else {
            write_rhea(item, rhea_equation_dict, rhea_draw)
        }
    });
}

function draw_rhea(rhea_code, rhea_img_dict, draw_div, width, height) {

        var lhs = rhea_img_dict[rhea_code][0]
        var rhs = rhea_img_dict[rhea_code][1]

        draw_div.innerHTML = ''
        for (const [key, value] of Object.entries(lhs)) {
            draw_div.innerHTML += "<div class='col padding-0'><div class='card border-0' style='height: 200px'>" +
                "<img class='card-img-top padding-0' width='"+width+"' height='"+height+"' src='" + value[2] +"' alt='" + value[0] + "'/>" + value[0] +
                "</div></div>"
        }
        draw_div.innerHTML += "<div class='col padding-0 my-auto text-center'><p class='no_margin padding-0 my-auto text-center'><b>=</b></p></div>"
        for (const [key, value] of Object.entries(rhs)) {
            draw_div.innerHTML += "<div class='col padding-0'><div class='card border-0' style='height: 200px'>" +
                "<img class='card-img-top padding-0' width='"+width+"' height='"+height+"' src='" + value[2] +"' alt='" + value[0] + "'/>" + value[0] +
                "</div></div>"
        }
    }

function write_rhea(rhea_code, rhea_equation_dict,  draw_div) {
    var equation = rhea_equation_dict[rhea_code]
    draw_div.innerHTML = "<p class='ml-3'>" + equation + "</p>"
}
