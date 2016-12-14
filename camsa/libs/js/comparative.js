/**
 * Created by aganezov on 12/10/16.
 */

function format_detailed_child_row(ap, sources, conflicts, assemblies_to_ids) {
    var table = '<table class="table table-condensed" cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">';
    table += '<thead><tr>';
    table += '<th class="text-center col-lg-1">Id</th>';
    // table += '<th>Sources</th>';
    table += '<th class="text-center col-lg-2">ISC</th>';
    table += '<th class="text-center col-lg-3">IC</th>';
    table += '<th class="text-center col-lg-3">OSC</th>';
    table += '<th class="text-center col-lg-3">OC</th>';
    table += '</tr></thead>';
    table += '<tbody>';

    table += '<tr>';
    table += '<td class="text-center">' + format_id(ap.self_id) + '</td>';
    table += '<td class="text-center">' + format_conflicting(conflicts[ap.self_id]['ISC'], assemblies_to_ids) + '</td>';
    table += '<td class="text-center">' + format_conflicting(conflicts[ap.self_id]['IC'], assemblies_to_ids) + '</td>';
    table += '<td class="text-center">' + format_conflicting(conflicts[ap.self_id]['OSC'], assemblies_to_ids) + '</td>';
    table += '<td class="text-center">' + format_conflicting(conflicts[ap.self_id]['OC'], assemblies_to_ids) + '</td>';
    table += '</tr>';

    table += '</tbody></table>';
    return table;
}

function format_id(ap_id) {
    return '<span ap_id="' + ap_id + '" class="ap_id_badge label label-info">' + ap_id + '</span>';
}

function get_sources_from_aps(aps) {
    var result = [];
    $.each(aps, function (i, ap) {
        result.push(ap.source);
    });
    return result;
}

function format_sources(sources) {
    var result = "";
    result += '[';
    $.each(sources, function (i, source) {

    });
    return result;
}

function format_conflicting(conflicts, assemblies_to_ids) {
    var result = '';
    if (Object.keys(conflicts).length == 0){
        return '<div class="alert alert-success">No conflicts!</div>'
    }
    for (var key in conflicts) {
        if (conflicts.hasOwnProperty(key)) {
            result += '<span>' + assemblies_to_ids[key] + '</span><br>';
            $.each(conflicts[key], function (i, a_conflict) {
                result += format_id(a_conflict) + ' ';
            });
            result += '<hr style="margin: 5px;">';
        }
    }
    return result;
}
