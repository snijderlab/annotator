"use strict";

const { invoke } = window.__TAURI__.core
var assets_folder = undefined;

function sortTable(id, column_number, type) {
    var table, rows, switching, i, x, y, shouldSwitch, dir = 0;
    table = document.getElementById(id);
    switching = true;
    dir = "desc";
    let headers = table.getElementsByTagName("tr")[0].getElementsByTagName("th");

    var sorted = false;
    if (headers[column_number].getAttribute('data-sort-order') == "asc") {
        dir = "desc";
        sorted = true
    }
    if (headers[column_number].getAttribute('data-sort-order') == "desc") {
        dir = "asc";
        sorted = true
    }

    for (var j = 0; j < headers.length; j++) {
        headers[j].setAttribute('data-sort-order', "");
    }
    headers[column_number].setAttribute('data-sort-order', dir);

    rows = Array.from(table.getElementsByTagName("tr"));
    var values = [null];
    for (i = 1; i < rows.length; i++) {
        x = rows[i].getElementsByTagName("td")[column_number]
        switch (type) {
            case "string":
                values.push(x.innerHTML.toLowerCase())
                break
            case "id":
                values.push(x.firstChild.innerHTML.toLowerCase())
                break
            case "number":
                values.push(Number(x.innerHTML))
                break
        }
    }

    if (sorted) {
        rows = rows.reverse()
        rows.splice(0, 0, rows[rows.length - 1])
        rows.pop()
        switching = false
    } else {
        switching = true
    }

    while (switching) {
        switching = false;
        /* Loop through all table rows (except the
        first, which contains table headers): */
        for (i = 1; i < (rows.length - 1); i++) {
            shouldSwitch = false;
            if (typeof values[i] == "string") {
                if (values[i].localeCompare(values[i + 1], undefined, { numeric: true, sensitivity: "base" }) == 1) {
                    shouldSwitch = true;
                    break;
                }
            } else if (values[i] < values[i + 1]) {
                shouldSwitch = true;
                break;
            }
        }
        if (shouldSwitch) {
            var el = rows[i + 1]
            rows.splice(i + 1, 1)
            rows.splice(i, 0, el)
            el = values[i + 1]
            values.splice(i + 1, 1)
            values.splice(i, 0, el)
            //rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
            switching = true;
        }
    }
    // Remove old rows and append all new rows
    while (table.lastChild) {
        table.removeChild(table.lastChild);
    }
    var frag = document.createDocumentFragment();
    for (var i = 0; i < rows.length; ++i) {
        frag.appendChild(rows[i]);
    }
    table.appendChild(frag);
}

function Select(id) {
    window.location.href = assets_folder + "/contigs/" + id.replace(':', '-') + ".html";
}

var dragging = false;

function get(name) {
    if (name = (new RegExp('[?&]' + encodeURIComponent(name) + '=([^&]*)')).exec(location.search))
        return decodeURIComponent(name[1]);
}

var t;
function Setup() {
    if (assets_folder != undefined) {
        var back_button = document.getElementById("back-button");
        var ref = window.name.split('|').pop();
        if (back_button && ref) {
            back_button.href = ref;
            var id = ref.split('/').pop().split('.')[0].replace(/-/g, ':');
            back_button.innerText = id;
            if (id != assets_folder.replace(/-/g, ':')) back_button.style = "display: inline-block;"
        }

        window.name = window.name + "|" + window.location.href

        if (window.name.split('|').pop().split('/').pop().split('.')[0].replace(/-/g, ':') == assets_folder.replace(/-/g, ':')) window.name = window.location.href;
    }

    if (get('pos')) {
        var pos = parseInt(get('pos'));
        var text = document.getElementsByClassName('aside-seq')[0].innerText;
        document.getElementsByClassName('aside-seq')[0].innerHTML = text.substring(0, pos) + "<span class='highlight'>" + text.substring(pos, pos + 1) + "</span>" + text.substring(pos + 1);
    }

    var elements = document.getElementsByClassName("copy-data")
    for (var i = 0; i < elements.length; i++) {
        elements[i].addEventListener("click", CopyGraphData)
        var example = elements[i].getElementsByClassName("example")[0]
        if (example)
            example.innerHTML = GetGraphDataExample(elements[i]);
    }

    var elements = document.getElementsByClassName("user-help")
    for (var i = 0; i < elements.length; i++) {
        if (!elements[i].classList.contains("copy-data")) {
            elements[i].addEventListener("click", (event) => { event.stopPropagation() })
        }
        elements[i].addEventListener("focusin", openHelp)
        elements[i].addEventListener("focusout", closeHelp)
        elements[i].addEventListener("mouseenter", openHelp)
        elements[i].addEventListener("mouseleave", closeHelp)
    }

    // Highlight the correct nodes and table row on the homepage after clicking on the point graph
    var split = window.location.href.split('#')
    if (split.length == 2) {
        var parts = split[1].split('_');
        var target = parts[parts.length - 1];
        var elements = document.querySelectorAll("[id$=\"" + target + "\"]")
        for (var i = 0; i < elements.length; i++) {
            elements[i].classList.add("highlight");
        }
    }

    SpectrumSetUp()
}

function openHelp(event) {
    var content = event.target.children[1]
    content.classList.add("focus");
    var box = content.getBoundingClientRect()
    if (box.x < 0) {
        content.style.marginLeft = -box.x + "px";
    }
    if (box.right > window.innerWidth) {
        content.style.marginLeft = (window.innerWidth - box.right - 20) + "px";
    }
    if (box.bottom > window.innerHeight) {
        content.style.marginTop = (window.innerHeight - box.bottom) + "px";
    }
}

function closeHelp(event) {
    var content = event.target.children[1]
    content.classList.remove("focus")
}

function pauseEvent(e) {
    if (e.stopPropagation) e.stopPropagation();
    if (e.preventDefault) e.preventDefault();
    e.cancelBubble = true;
    e.returnValue = false;
    return false;
}

function GoBack() {
    var history = window.name.split('|');
    history.pop(); //remove current
    var url = history.pop();
    window.name = history.join('|');
    window.location.href = url;
}

function AlignmentDetails(number) {
    AlignmentDetailsClear();
    document.getElementById("alignment-details-" + number).className += " active";
}

function AlignmentDetailsClear() {
    document.getElementById("index-menus").childNodes.forEach(a => a.className = "alignment-details");
}

// Copy the graph data for the clicked on graph
function CopyGraphData(event) {
    var t = event.target;
    if (event.target.classList.contains("mark"))
        t = event.target.parentElement.parentElement;
    if (event.target.classList.contains("copy-data"))
        t = event.target.parentElement;

    var elements = t.getElementsByClassName("graph-data");
    for (var i = 0; i < elements.length; i++) {
        elements[i].select();
        elements[i].setSelectionRange(0, 99999);
        navigator.clipboard.writeText(elements[i].value);
        return;
    }
}

/** Get the first four lines of the data to work as an example of how the data looks.
 * @param {Element} element the .copy-data element to get the example for
*/
function GetGraphDataExample(element) {
    var children = element.children;
    if (element.classList.contains("mark"))
        children = element.parentElement.parentElement.children;
    if (element.classList.contains("copy-data"))
        children = element.parentElement.children;

    var example = undefined;
    for (var i = 0; i < children.length; i++) {
        if (children[i].classList.contains("graph-data")) {
            example = children[i].value.split("\n").slice(0, 4);
            break;
        }
    }

    if (example == undefined) {
        return "No data found";
    } else {
        var output = "";
        for (var line in example) {
            output += example[line].substring(0, Math.min(30, example[line].length)).replace(/ /gi, "<span class='sp'> </span>").replace(/\t/gi, "<span class='tab'>\t</span>")
            if (example[line].length <= 30) {
                output += "<span class='nl'></span>\n";
            } else {
                output += "...\n";
            }
        }
        return output;
    }
}

function ToggleCDRReads() {
    document.querySelector(".alignment-body").classList.toggle("cdr");
}

function ToggleAlignmentComic() {
    document.querySelector(".alignment-body").classList.toggle("comic");
}

/** Update the reads alignment to only show the read supporting the given position. Or if only this position is shown remove the highlight.
 * @param {String} position the class of the position eg for position 42 this would be "a42".
 * @param {Number} alignment_position the position of this node in the template alignment.
*/
function HighlightAmbiguous(position, alignment_position) {
    var body = document.querySelector(".alignment-body")
    if (body.classList.contains("highlight") && body.dataset.position == position) {
        body.style.setProperty("--highlight-pos", "-1ch");
        body.dataset.position = undefined;
        document.querySelectorAll(".highlight").forEach(e => e.classList.remove("highlight")); // Removes highlight from all reads, the .node, and .alignment-body
    } else {
        document.querySelectorAll(".highlight").forEach(e => e.classList.remove("highlight"));
        body.style.setProperty("--highlight-pos", String(alignment_position) + "ch");
        body.dataset.position = position;
        body.classList.add("highlight")
        document.querySelector(".higher-order-graphs").classList.add("highlight")
        document.querySelectorAll("." + position).forEach(e => e.classList.add("highlight")); // Highlights all reads and the .node
    }
}

/* Spectrum viewer code */
function SpectrumSetUp() {
    // Settings
    var elements = document.querySelectorAll("#spectrum-wrapper .graphics-settings input");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("change", GraphicsSettings);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .graphics-settings button");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("click", GraphicsSettings);
    }
    document.querySelectorAll("#spectrum-wrapper .error-graph-settings input")
        .forEach(ele => ele.addEventListener("change", ErrorGraphSettings));
    var elements = document.querySelectorAll("#spectrum-wrapper .error-graph-settings .manual-zoom input");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("change", ManualZoomSpectrumGraph);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .peptide-settings input");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("change", PeptideSettings);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .peptide-settings button");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("click", PeptideSettings);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .spectrum-settings input");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("change", SpectrumSettings);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .spectrum-settings button");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("click", SpectrumSettings);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .spectrum-settings .manual-zoom input");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("change", ManualZoom);
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .spectrum-settings .reset-zoom");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("mouseup", spectrumZoomOut)
        elements[i].addEventListener("keyup", e => { if (e.key == "Enter") spectrumZoomOut(e) })
    }

    // Legend
    var elements = document.querySelectorAll("#spectrum-wrapper .legend .ion");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("mouseenter", e => { ToggleHighlight(e.target, false, true) });
        elements[i].addEventListener("mouseleave", e => { ToggleHighlight(e.target, false, false) });
        elements[i].addEventListener("focusin", e => { ToggleHighlight(e.target, false, true) });
        elements[i].addEventListener("focusout", e => { ToggleHighlight(e.target, false, false) });
        elements[i].addEventListener("mouseup", e => {
            ToggleHighlight(e.target, true, !e.target.classList.contains("permanent"), null)
        });
        elements[i].addEventListener("keyup", e => { if (e.key == "Enter") e.target.click() });
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .legend .unassigned");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("keyup", e => { if (e.key == "Enter") e.target.click() });
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .legend .legend-peptide");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("keyup", e => { if (e.key == "Enter") e.target.click() });
    }
    var elements = document.querySelectorAll("#spectrum-wrapper .legend input");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("change", SpectrumInputChange);
    }

    // Peptide
    var elements = document.querySelectorAll("#spectrum-wrapper .peptide span");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("mouseenter", e => { SequenceElementEvent(e, false, true) });
        elements[i].addEventListener("mouseleave", e => { SequenceElementEvent(e, false, false) });
        elements[i].addEventListener("focusin", e => { SequenceElementEvent(e, false, true) });
        elements[i].addEventListener("focusout", e => { SequenceElementEvent(e, false, false) });
        elements[i].addEventListener("mousedown", e => { SequenceElementEvent(e, false, true) });
        elements[i].addEventListener("mouseup", e => { SequenceElementEvent(e, true) });
        elements[i].addEventListener("keyup", e => { if (e.key == "Enter") e.target.click() });
    }

    // Spectrum
    document.querySelectorAll("#spectrum-wrapper .wrapper").forEach(element => {
        element.addEventListener("mousedown", spectrumDragStart)
        element.addEventListener("mousemove", spectrumDrag)
        element.addEventListener("mouseup", spectrumDragEnd)
        element.addEventListener("mouseout", spectrumDragOut)
    })
    document.querySelectorAll("#spectrum-wrapper .peak").forEach(element => {
        element.addEventListener("mousedown", spectrumDistanceStart)
        element.addEventListener("mouseup", spectrumDistanceEnd)
    })
    document.querySelectorAll("#spectrum-wrapper .canvas-wrapper").forEach(element => {
        var d = element.dataset;
        d.minMz = 0;
        d.maxMz = d.initialMaxMz;
        d.maxIntensity = d.initialMaxIntensity;
    })
    document.querySelectorAll("#spectrum-wrapper .peak:not(.fragment)").forEach(element => {
        element.addEventListener("click", ForceShowLabels)
    })
    document.querySelectorAll("#spectrum-wrapper .canvas-spectrum").forEach(element => {
        element.addEventListener("wheel", spectrumScroll)
    })
    document.addEventListener("keyup", e => {
        if (e.key == "0" && e.ctrlKey) {
            spectrumZoomOut(e)
        } else if ((e.key == "=" || e.key == "-") && e.ctrlKey) {
            const wrapper = document.querySelector(".canvas-wrapper");
            const minMz = Number(wrapper.dataset.minMz);
            const maxMz = Number(wrapper.dataset.maxMz);
            const delta = 0.05 * (maxMz - minMz);
            const direction = (e.key == "-") ? -1 : 1;
            Zoom(wrapper, Math.max(0, minMz + direction * delta * 0.5), maxMz - direction * delta * 0.5, Number(wrapper.dataset.maxIntensity));
        }
    })

    // Spectrum graph
    var elements = document.querySelectorAll("#spectrum-wrapper .error-graph");
    for (let i = 0; i < elements.length; i++) {
        elements[i].addEventListener("mousemove", SpectrumGraphMouseMove)
    }

    // Make the axes look nice
    UpdateSpectrumAxes(document.querySelector("#spectrum-wrapper .canvas-wrapper"));

    // Update error graph
    UpdateErrorGraph()
}

/**
 * Force show the labels of a peak (click event)
 * @param {MouseEvent} e The event
 */
function ForceShowLabels(e) {
    let spectrum = document.getElementById("spectrum-wrapper");
    if (spectrum.classList.contains("force-show-label")) {
        if ("showLabelManually" in e.target.dataset) {
            delete e.target.dataset.showLabelManually;
        } else {
            e.target.dataset.showLabelManually = "true";
        }
    } else if (spectrum.classList.contains("force-show-m-z")) {
        if ("showMZManually" in e.target.dataset) {
            delete e.target.dataset.showMZManually;
        } else {
            e.target.dataset.showMZManually = "true";
        }
    } else if (spectrum.classList.contains("force-show-hide")) {
        e.target.dataset.showLabelManually = "false";
        e.target.dataset.showMZManually = "false";
    }
}


/**
 * @type {EventTarget}
 */
let sequence_element_start = undefined;
/**
 * Do the event for a sequence element (one element in the peptide)
 * @param {MouseEvent} e The event
 * @param {boolean} permanent If it will be applied permanently (true if clicked, false if hovered)
 * @param {boolean?} turn_on If this is the turns on (true) or turns off (false) the event, the default is toggle (null)
 */
function SequenceElementEvent(e, permanent, turn_on = null) {
    // if (selected_colour != "highlight") {
    if (e.type == "mouseup" && sequence_element_start != undefined && sequence_element_start.dataset.pos != e.target.dataset.pos) {
        let state = sequence_element_start.classList.contains(selected_colour);
        let start = sequence_element_start.dataset.pos.split("-");
        let end = e.target.dataset.pos.split("-");
        if (start[0] == end[0] && start[1] == end[1]) {
            let range = [Math.min(Number(start[2]), Number(end[2])), Math.max(Number(start[2]), Number(end[2]))];
            document.querySelectorAll(".spectrum .peptide>span[title]").forEach(element => {
                let pos = element.dataset.pos.split("-");
                element.classList.remove("select");
                if (pos[0] == start[0] && pos[1] == start[1] && Number(pos[2]) >= range[0] && Number(pos[2]) <= range[1]) {
                    element.classList.remove("red", "green", "blue", "yellow", "purple", "default", "highlight");
                    ToggleHighlight(element, true, (selected_colour == "remove" ? false : !state), selected_colour);//, ".p" + pos[0] + "-" + pos[1]);
                    if (selected_colour != "remove") {
                        element.classList.toggle(selected_colour, !state);
                        element.classList.toggle("highlight", !state);
                    }
                }
            })
        }
        sequence_element_start = undefined;
    } else if (permanent) {
        let state = e.target.classList.contains(selected_colour);
        e.target.classList.remove("red", "green", "blue", "yellow", "purple", "default", "highlight");
        ToggleHighlight(e.target, true, (selected_colour == "remove" ? false : !state), selected_colour);
        if (selected_colour != "remove") {
            e.target.classList.toggle(selected_colour, !state);
            e.target.classList.toggle("highlight", !state);
        }
        sequence_element_start = undefined;
    } else if (e.type == "mousedown") {
        sequence_element_start = e.target;
    } else if (e.type == "mouseenter" && sequence_element_start != undefined) {
        let start = sequence_element_start.dataset.pos.split("-");
        let end = e.target.dataset.pos.split("-");
        if (start[0] == end[0] && start[1] == end[1]) {
            let range = [Math.min(Number(start[2]), Number(end[2])), Math.max(Number(start[2]), Number(end[2]))];
            document.querySelectorAll(".spectrum .peptide>span[title]").forEach(element => {
                let pos = element.dataset.pos.split("-");
                if (pos[0] == start[0] && pos[1] == start[1] && Number(pos[2]) >= range[0] && Number(pos[2]) <= range[1]) {
                    element.classList.add("select");
                } else {
                    element.classList.remove("select");
                }
            })
        }
    } else if (sequence_element_start == undefined && e.type == "mouseenter") {
        let state = e.target.classList.contains(selected_colour);
        ToggleHighlight(e.target, false, true, selected_colour);
    } else if (sequence_element_start == undefined && e.type == "mouseleave") {
        let state = e.target.classList.contains(selected_colour);
        ToggleHighlight(e.target, false, false, selected_colour);
    }
}

export function SetUpSpectrumInterface() {
    document.getElementById("spectrum-wrapper").classList.add("error-graph-relative");
    SpectrumSetUp();
    var event = new Event('change');
    document.getElementById("error-graph-y-max").dispatchEvent(event);
}

var highlight;
/**
 * Toggle a highlight
 * @param {Element} t The target
 * @param {boolean} permanent If it will be applied permanently (true if clicked, false if hovered)
 * @param {boolean?} state If this is the apply (true) or clear (false) operation
 * @param {string?} selected_colour For peptide highlights the colour selected or "default"
 */
function ToggleHighlight(t, permanent, state, selected_colour) {
    let current_state = t.classList.contains("permanent");
    if (permanent) t.classList.toggle("permanent", state);
    if (t.classList.contains("permanent") && !permanent) return;

    highlight = t;
    var container = document.getElementById("spectrum-wrapper");

    // Set the correct state for a single peak
    function toggle(element) {
        // Make sure the n is set, this tracks the number of reasons why a peak should be highlighted
        if (element.dataset.n == undefined) element.dataset.n = 0;
        if (element.dataset.n < 0) element.dataset.n = 0;

        // Clear any colour set for this peak
        if (selected_colour != null && permanent)
            element.classList.remove("red", "green", "blue", "yellow", "purple", "default");

        if (permanent) {
            // A permanent highlight, where this position or ion type is clicked
            if (state === true) {
                // If this position in the peptide was already highlighted do not add one to n
                if (current_state != state) element.dataset.n = Number(element.dataset.n) + 1;

                // Set up the classes
                if (element.dataset.n == 1) element.classList.add("highlight");
                if (selected_colour != null) element.classList.add(selected_colour);
            }
            else if (state === false) {
                if (current_state != state) element.dataset.n = Number(element.dataset.n) - 1;

                // Set up the classes
                if (element.dataset.n <= 0) element.classList.remove("highlight");
                if (selected_colour != null) element.classList.remove(selected_colour);
            } else {
                console.error("When using 'permanent' the state should be set");
            }
        } else {
            // A temporary highlight, while hovering over this position or ion type
            if (element.dataset.n == undefined || element.dataset.n == 0) {
                if (state === true) {
                    element.classList.add("highlight");
                    if (selected_colour != null) element.classList.add(selected_colour);
                } else if (state === null) {
                    element.classList.toggle("highlight");
                    if (selected_colour != null) element.classList.toggle(selected_colour);
                } else {
                    element.classList.remove("highlight");
                    if (selected_colour != null) element.classList.remove(selected_colour);
                }
            }
        }
    }

    // Select all peaks for this operation
    container.querySelectorAll(".canvas").forEach(canvas => {
        let selector = "";
        if (t.classList.contains("n-term")) {
            selector = ":is(.a, .b, .c, .d, .v)";
        } else if (t.classList.contains("c-term")) {
            selector = ":is(.w, .x, .y, .z)";
        } else if (t.classList.contains("ion")) {
            var cl = t.classList[1];
            selector = "." + cl;
        } else if (t.classList.contains("name")) {
            selector = ".p" + t.dataset.pos;
        } else {
            selector = ".p" + t.dataset.pos;
        }
        canvas.querySelectorAll(selector).forEach(element => toggle(element))
    })

    // Set the state for the whole canvas (is there anything highlighted still)
    if (state === true || container.querySelector(".permanent") != null) {
        container.querySelectorAll(".canvas").forEach(canvas => {
            canvas.classList.add("highlight");
        })
    } else {
        container.querySelectorAll(".canvas").forEach(canvas => {
            canvas.classList.remove("highlight");
        })
    }
}

/** The mouse moved in the spectrum graph.
 * @param {MouseEvent} event
*/
function SpectrumGraphMouseMove(event) {
    var data = event.target;
    if (data.classList.contains("point")) data = data.parentElement;
    var spectrum = data.parentElement;

    var baseline = data.getBoundingClientRect().top;
    var offset = event.clientY - baseline;

    var el = spectrum.querySelector(".ruler");
    el.style.top = offset + "px";
    var el = spectrum.querySelector("#ruler-value");

    let style = getComputedStyle(spectrum);
    let min = Number(style.getPropertyValue(`--y-min`));
    let max = Number(style.getPropertyValue(`--y-max`));
    el.innerText = fancyRound(max, min, offset / data.clientHeight * - 1 * (max - min) + max, 1);
}

/** Zoom the spectrum manually.
 * @param {Event} event
*/
function ManualZoom(event) {
    var spectrum = event.target.parentElement.parentElement;
    var min = Number(spectrum.querySelector(".mz-min").value);
    var max = Number(spectrum.querySelector(".mz-max").value);
    var maxI = Number(spectrum.querySelector(".intensity-max").value);

    spectrum.querySelectorAll(".canvas-wrapper").forEach(canvas => Zoom(canvas, min, max, maxI));
}

/** Zoom the spectrum manually.
 * @param {Event} event
*/
function ManualZoomSpectrumGraph(event) {
    var min_y = Number(document.getElementById("error-graph-y-min").value);
    var max_y = Number(document.getElementById("error-graph-y-max").value);

    document.getElementById("spectrum-wrapper")
        .querySelectorAll(".error-graph").forEach(canvas => ZoomSpectrumGraph(canvas, min_y, max_y));
}

/** Setup properties of the spectrum for publication
 * @param {Event} event
*/
function GraphicsSettings(event) {
    var cl = event.target.className;
    let spectrum_wrapper = document.getElementById("spectrum-wrapper");
    if (cl == "width") {
        spectrum_wrapper.style.setProperty("--width", event.target.value);
    } else if (cl == "height") {
        spectrum_wrapper.style.setProperty("--height", event.target.value);
    } else if (cl == "fs-peptide") {
        spectrum_wrapper.style.setProperty("--fs-peptide", event.target.value);
    } else if (cl == "stroke-peptide") {
        spectrum_wrapper.style.setProperty("--stroke-peptide", event.target.value);
    } else if (cl == "fs-spectrum") {
        spectrum_wrapper.style.setProperty("--fs-spectrum", event.target.value);
    } else if (cl == "stroke-spectrum") {
        spectrum_wrapper.style.setProperty("--stroke-spectrum", event.target.value);
    } else if (cl == "stroke-spectrum-unassigned") {
        spectrum_wrapper.style.setProperty("--stroke-spectrum-unassigned", event.target.value);
    }
}

/** Setup properties of the spectrum for publication
 * @param {Event} event
*/
function ErrorGraphSettings(event) {
    var id = event.target.id;
    let spectrum_wrapper = document.getElementById("spectrum-wrapper");
    if (id == "error-graph-y-type-absolute") {
        spectrum_wrapper.classList.toggle("error-graph-absolute");
        spectrum_wrapper.classList.remove("error-graph-relative");
        UpdateErrorGraph()
    } else if (id == "error-graph-y-type-relative") {
        spectrum_wrapper.classList.remove("error-graph-absolute");
        spectrum_wrapper.classList.toggle("error-graph-relative");
        UpdateErrorGraph()
    } else if (id == "error-graph-intensity") {
        spectrum_wrapper.classList.toggle("error-graph-intensity");
    } else {
        UpdateErrorGraph()
    }
}

/**
 * Update for the whole error graph the --y values for all points
 */
function UpdateErrorGraph() {
    let spectrum_wrapper = document.getElementById("spectrum-wrapper");
    let graph = document.getElementById("error-graph");
    let points = graph.querySelectorAll(".point");
    let relative = spectrum_wrapper.classList.contains("error-graph-relative");
    let assigned_mode = document.getElementById("error-graph-x-type-assigned").checked;
    let use_a = document.getElementById("error-graph-ion-a").checked;
    let use_b = document.getElementById("error-graph-ion-b").checked;
    let use_c = document.getElementById("error-graph-ion-c").checked;
    let use_x = document.getElementById("error-graph-ion-x").checked;
    let use_y = document.getElementById("error-graph-ion-y").checked;
    let use_z = document.getElementById("error-graph-ion-z").checked;
    let show_assigned = document.getElementById("error-graph-show-assigned").checked;
    document.getElementById("error-graph-y-title").innerText =
        "The error for all " +
        (assigned_mode ? "assigned " : (show_assigned ? "" : "unassigned ")) +
        "peaks as compared to " +
        (assigned_mode ? "the annotation" : "the closest theoretical ion") +
        " (" + (relative ? "ppm" : "Da") + ")";
    let points_y = Array();

    points.forEach((point) => {
        let y = Infinity;
        let label = "";
        if (assigned_mode) {
            if (relative) {
                y = Number(point.dataset.aRel);
            } else {
                y = Number(point.dataset.aAbs);
            }
        } else {
            if (!show_assigned && point.dataset.aRel != "undefined") {
                // hide
            } else if (relative) {
                if (use_a && Math.abs(y) > Math.abs(Number(point.dataset.uARelValue))) {
                    y = Number(point.dataset.uARelValue);
                    label = point.dataset.uARelFragment
                }
                if (use_b && Math.abs(y) > Math.abs(Number(point.dataset.uBRelValue))) {
                    y = Number(point.dataset.uBRelValue);
                    label = point.dataset.uBRelFragment
                }
                if (use_c && Math.abs(y) > Math.abs(Number(point.dataset.uCRelValue))) {
                    y = Number(point.dataset.uCRelValue);
                    label = point.dataset.uCRelFragment
                }
                if (use_x && Math.abs(y) > Math.abs(Number(point.dataset.uXRelValue))) {
                    y = Number(point.dataset.uXRelValue);
                    label = point.dataset.uXRelFragment
                }
                if (use_y && Math.abs(y) > Math.abs(Number(point.dataset.uYRelValue))) {
                    y = Number(point.dataset.uYRelValue);
                    label = point.dataset.uYRelFragment
                }
                if (use_z && Math.abs(y) > Math.abs(Number(point.dataset.uZRelValue))) {
                    y = Number(point.dataset.uZRelValue);
                    label = point.dataset.uZRelFragment
                }
            } else {
                if (use_a && Math.abs(y) > Math.abs(Number(point.dataset.uAAbsValue))) {
                    y = Number(point.dataset.uAAbsValue);
                    label = point.dataset.uAAbsFragment
                }
                if (use_b && Math.abs(y) > Math.abs(Number(point.dataset.uBAbsValue))) {
                    y = Number(point.dataset.uBAbsValue);
                    label = point.dataset.uBAbsFragment
                }
                if (use_c && Math.abs(y) > Math.abs(Number(point.dataset.uCAbsValue))) {
                    y = Number(point.dataset.uCAbsValue);
                    label = point.dataset.uCAbsFragment
                }
                if (use_x && Math.abs(y) > Math.abs(Number(point.dataset.uXAbsValue))) {
                    y = Number(point.dataset.uXAbsValue);
                    label = point.dataset.uXAbsFragment
                }
                if (use_y && Math.abs(y) > Math.abs(Number(point.dataset.uYAbsValue))) {
                    y = Number(point.dataset.uYAbsValue);
                    label = point.dataset.uYAbsFragment
                }
                if (use_z && Math.abs(y) > Math.abs(Number(point.dataset.uZAbsValue))) {
                    y = Number(point.dataset.uZAbsValue);
                    label = point.dataset.uZAbsFragment
                }
            }
        }
        if (Number.isNaN(y) || y == Infinity) {
            point.classList.toggle("hidden", true);
        } else {
            point.classList.toggle("hidden", false);
            point.style.setProperty("--y", y);
            point.dataset.label = label;
            points_y.push(y);
        }
    })

    // Call 
    invoke("density_graph", { data: points_y }).then((result) => {
        document.getElementById("error-graph-density").innerHTML = result;
    }).catch((error) => {
        console.error("density estimation wrong", error)
    })
}

/** Setup properties of the spectrum for publication
 * @param {Event} event
*/
let selected_colour = "default";
function PeptideSettings(event) {
    var t = event.target;
    let spectrum_wrapper = document.getElementById("spectrum-wrapper");
    if (t.id == "spectrum-compact") {
        spectrum_wrapper.classList.toggle("compact");
    } else if (t.name == "highlight") {
        selected_colour = t.value;
    } else if (t.id == "clear-colour") {
        spectrum_wrapper.querySelectorAll(".peptide>span").forEach(item => item.classList.remove("red", "green", "blue", "yellow", "purple", "default", "highlight"))
        spectrum_wrapper.querySelectorAll(".ion-series .permanent").forEach(item => item.classList.remove("permanent"))
        spectrum_wrapper.querySelectorAll(".peak").forEach(item => { item.classList.remove("red", "green", "blue", "yellow", "purple", "default", "highlight"); item.dataset.n = 0; })
        spectrum_wrapper.querySelectorAll(".point").forEach(item => { item.classList.remove("red", "green", "blue", "yellow", "purple", "default", "highlight"); item.dataset.n = 0; })
        spectrum_wrapper.querySelectorAll(".canvas").forEach(item => item.classList.remove("highlight"))
    }
}

/** Setup properties of the spectrum for publication
 * @param {Event} event
*/
function SpectrumSettings(event) {
    var t = event.target;
    var cl = t.className;
    let spectrum_wrapper = document.getElementById("spectrum-wrapper");
    let canvas = spectrum_wrapper.querySelector(".canvas-wrapper");
    if (cl == "theoretical") {
        spectrum_wrapper.classList.toggle("theoretical");
    } else if (t.id == "peak-colour-ion") {
        spectrum_wrapper.classList.add("legend-ion");
        spectrum_wrapper.classList.remove("legend-peptide", "legend-peptidoform", "legend-none");
    } else if (t.id == "peak-colour-peptide") {
        spectrum_wrapper.classList.remove("legend-ion", "legend-peptidoform", "legend-none");
        spectrum_wrapper.classList.add("legend-peptide");
    } else if (t.id == "peak-colour-peptidoform") {
        spectrum_wrapper.classList.remove("legend-ion", "legend-peptide", "legend-none");
        spectrum_wrapper.classList.add("legend-peptidoform");
    } else if (t.id == "peak-colour-none") {
        spectrum_wrapper.classList.remove("legend-ion", "legend-peptide", "legend-peptidoform");
        spectrum_wrapper.classList.add("legend-none");
    } else if (t.id == "peak-distance") {
        if (spectrum_wrapper.classList.contains("ruler")) {
            spectrum_wrapper.classList.remove("ruler");
            ruler_mode = false;
        } else {
            spectrum_wrapper.classList.add("ruler");
            ruler_mode = true;
        }
    } else if (t.id == "force-show-none") {
        spectrum_wrapper.classList.remove("force-show-label", "force-show-m-z", "force-show-hide");
    } else if (t.id == "force-show-label") {
        spectrum_wrapper.classList.remove("force-show-m-z", "force-show-hide");
        spectrum_wrapper.classList.add("force-show-label");
    } else if (t.id == "force-show-m-z") {
        spectrum_wrapper.classList.remove("force-show-label", "force-show-hide");
        spectrum_wrapper.classList.add("force-show-m-z");
    } else if (t.id == "force-show-hide") {
        spectrum_wrapper.classList.remove("force-show-label", "force-show-m-z");
        spectrum_wrapper.classList.add("force-show-hide");
    } else if (t.id == "force-show-clear") {
        spectrum_wrapper.querySelectorAll(".peak").forEach(element => { delete element.dataset.showLabelManually; delete element.dataset.showMZManually; });
    } else if (cl == "unassigned") {
        if (t.checked) { // Will be adding all background peaks
            spectrum_wrapper.classList.add('show-unassigned');
            if (Number(canvas.dataset.maxIntensity) == Number(canvas.dataset.initialMaxIntensityAssigned)) {
                canvas.dataset.maxIntensity = canvas.dataset.initialMaxIntensity;
                canvas.style.setProperty("--max-intensity", canvas.dataset.maxIntensity);
                UpdateSpectrumAxes(canvas)
            }
        } else { // Will be removing all background peaks
            spectrum_wrapper.classList.remove('show-unassigned');
            if (canvas.dataset.maxIntensity == canvas.dataset.initialMaxIntensity) {
                canvas.dataset.maxIntensity = canvas.dataset.initialMaxIntensityAssigned;
                canvas.style.setProperty("--max-intensity", canvas.dataset.maxIntensity);
                UpdateSpectrumAxes(canvas)
            }
        }
    } else if (t.id == "spectrum-label-charge") {
        spectrum_wrapper.classList.toggle("show-charge", t.checked);
    } else if (t.id == "spectrum-label-series") {
        spectrum_wrapper.classList.toggle("show-series", t.checked);
    } else if (t.id == "spectrum-label-glycan-id") {
        spectrum_wrapper.classList.toggle("show-glycan-id", t.checked);
    } else if (t.id == "spectrum-label-peptide-id") {
        spectrum_wrapper.classList.toggle("show-peptide-id", t.checked);
    } else if (t.id == "spectrum-label-neutral-losses") {
        spectrum_wrapper.classList.toggle("show-neutral-losses", t.checked);
    } else if (t.id == "spectrum-label-cross-links") {
        spectrum_wrapper.classList.toggle("show-cross-links", t.checked);
    } else if (t.id == "spectrum-label-ambiguous-amino-acids") {
        spectrum_wrapper.classList.toggle("show-ambiguous-amino-acids", t.checked);
    } else if (t.id == "spectrum-label-modifications") {
        spectrum_wrapper.classList.toggle("show-modifications", t.checked);
    } else if (t.id == "spectrum-label-charge-carriers") {
        spectrum_wrapper.classList.toggle("show-charge-carriers", t.checked);
    } else if (t.id == "rotate-label") {
        spectrum_wrapper.classList.toggle("rotate-label", t.checked);
    } else if (t.id == "spectrum-mz-min") {
        canvas.dataset.minMz = Number(t.value);
        canvas.style.setProperty("--min-mz", canvas.dataset.minMz);
        UpdateSpectrumAxes(canvas)
    } else if (t.id == "spectrum-mz-max") {
        canvas.dataset.maxMz = Number(t.value);
        canvas.style.setProperty("--max-mz", canvas.dataset.maxMz);
        UpdateSpectrumAxes(canvas)
    } else if (t.id == "spectrum-intensity-max") {
        canvas.dataset.maxIntensity = Number(t.value);
        canvas.style.setProperty("--max-intensity", canvas.dataset.maxIntensity);
        UpdateSpectrumAxes(canvas)
    } else if (t.id == "spectrum-ticks-y" || t.id == "spectrum-ticks-x") {
        UpdateSpectrumAxes(canvas)
    } else if (t.type == "range" || t.type == "number") { // Peak labels + masses
        var ele = document.getElementById(
            (t.type == "number") ?
                t.id.substring(0, t.id.length - 6)
                : t.id + "-value");
        ele.value = t.value;

        spectrum_wrapper.querySelectorAll(".canvas-spectrum").forEach(
            canvas => SpectrumUpdateLabels(canvas))
    } else if (t.id == "y-sqrt") {
        spectrum_wrapper.classList.toggle("y-sqrt", t.checked);
        UpdateSpectrumAxes(canvas)
    } else if (t.id == "y-percentage") {
        spectrum_wrapper.classList.toggle("y-percentage", t.checked);
        UpdateSpectrumAxes(canvas)
    }
}

/** Update the spectrum to only show the label for peaks within the given percentage group
 * @param {Element} canvas the canvas ('.canvas') to work on
*/
function SpectrumUpdateLabels(canvas) {
    const max = Number(canvas.parentElement.style.getPropertyValue("--max-intensity"));
    const label_value = Number(document.getElementById("spectrum-label-value").value);
    const mz_value = Number(document.getElementById("spectrum-m-z-value").value);
    const label_threshold = (100 - label_value) / 100 * max;
    const mz_threshold = (100 - mz_value) / 100 * max;

    setTimeout(() => {
        for (let index = 1; index < canvas.children.length; index++) {
            const peak = canvas.children[index];
            const intensity = Number(peak.style.getPropertyValue("--intensity"));
            // Update shown labels
            if (label_value != 0 && intensity >= label_threshold) {
                peak.dataset.showLabel = "true";
            } else {
                peak.dataset.showLabel = "";
            }
            if (mz_value != 0 && intensity >= mz_threshold) {
                peak.dataset.showMZ = "true";
            } else {
                peak.dataset.showMZ = "";
            }
            // Update cut peaks
            if (intensity > max) {
                peak.dataset.cut = "true";
            } else {
                peak.dataset.cut = "";
            }
        }
    }, 0)
}

var startPoint;
var selection;
var linked_selection;
var last_selection;
var ruler_mode = false;

function spectrumDragStart(event) {
    var wrapper = event.target;
    if (wrapper.classList.contains("peptide") || wrapper.classList.contains("canvas-wrapper")) wrapper = wrapper.parentElement;
    if (wrapper.classList.contains("canvas") || wrapper.classList.contains("y-axis") || wrapper.classList.contains("x-axis")) wrapper = wrapper.parentElement.parentElement;
    if (wrapper.classList.contains("wrapper")) {
        var canvas = wrapper.querySelector(".canvas");
        startPoint = event.pageX - wrapper.getBoundingClientRect().x;
        selection = canvas.querySelector(".selection");
        if (wrapper.classList.contains("first")) {
            linked_selection = wrapper.parentElement.querySelector(".selection.second")
        } else if (wrapper.classList.contains("second")) {
            linked_selection = wrapper.parentElement.querySelector(".selection.first")
        }
        wrapper.classList.add("dragging")
        selection.hidden = false;
        if (linked_selection != undefined) {
            linked_selection.hidden = false;
            linked_selection.classList.add("linked");
        }
        updateSelections(wrapper, canvas, startPoint, event.pageY - wrapper.getBoundingClientRect().y, startPoint);
    }
}

function updateSelections(wrapper, canvas, x, y, end_x) {
    var wrapper_box = wrapper.getBoundingClientRect();
    var canvas_box = canvas.getBoundingClientRect();
    var new_x = Math.min(Math.max(0, x - (canvas_box.left - wrapper_box.left)), canvas_box.width);
    var new_y = Math.min(Math.max(1, y - (canvas_box.top - wrapper_box.top)), canvas_box.height);
    var new_end_x = Math.min(Math.max(0, end_x - (canvas_box.left - wrapper_box.left)), canvas_box.width);
    var new_w = Math.min(Math.abs(new_end_x - new_x), canvas_box.width);
    var new_x = Math.min(new_x, new_end_x);
    if (selection != undefined) updateSelection(selection, new_x, new_y, new_w);
    if (linked_selection != undefined) updateSelection(linked_selection, new_x, new_y, new_w);
}


function updateSelection(select, offsetX, offsetY, width) {
    select.style.setProperty("left", offsetX + "px")
    select.style.setProperty("width", width + "px")
    select.style.setProperty("height", select.parentElement.getBoundingClientRect().height - offsetY + "px")

    if (select.classList.contains("first")) {
        select.style.setProperty("top", offsetY + "px")
    } else {
        select.style.setProperty("bottom", offsetY + "px")
    }
}

function spectrumDrag(event) {
    if (startPoint != undefined) {
        var canvas = event.target.querySelector(".canvas");
        var box = canvas.getBoundingClientRect();
        var offsetY = event.offsetY;
        if (selection.classList.contains("second")) offsetY = box.height - offsetY;
        updateSelections(event.target, canvas, startPoint, offsetY, event.offsetX)
    } else if (last_selection != undefined && event.target.classList.contains("canvas")) {
        if (last_selection.canvas == event.target) {
            selection = last_selection.selection;
            linked_selection = last_selection.linked_selection;
        } else {
            selection = last_selection.linked_selection;
            linked_selection = last_selection.selection;
        }
        startPoint = last_selection.startPoint;
        if (selection != undefined) selection.hidden = false;
        if (linked_selection != undefined) linked_selection.hidden = false;
        event.target.parentElement.parentElement.classList.add("dragging");
        last_selection = undefined;
    }
}

function spectrumDragOut(event) {
    if (event.target.classList.contains("wrapper")) {
        if (selection != undefined)
            last_selection = { selection: selection, linked_selection: linked_selection, startPoint: startPoint, canvas: event.target.querySelector(".canvas") };
        else last_selection = undefined;
        if (selection != undefined) selection.hidden = true;
        if (linked_selection != undefined) {
            linked_selection.hidden = true;
            linked_selection.classList.remove("linked");
        }
        selection = undefined;
        linked_selection = undefined;
        startPoint = undefined;
        first_anchor = undefined;

        document.querySelector("#spectrum-wrapper").classList.remove("distance");
        document.querySelectorAll(".wrapper.dragging").forEach(e => e.classList.remove("dragging"));
    }
}

function spectrumDragEnd(event) {
    var canvas = event.target.querySelector(".canvas-wrapper");
    if (startPoint != undefined) {
        var d = canvas.dataset;
        var box = canvas.querySelector(".canvas").getBoundingClientRect();
        var wrapper_box = event.target.getBoundingClientRect();
        var width = box.width;
        var height = box.height;
        var boxY = box.y + window.scrollY;
        var boxX = box.x + window.scrollX;
        startPoint = startPoint + wrapper_box.x - boxX;
        var min = Math.max(0, Math.min(startPoint, event.pageX - boxX)) / width;
        var max = Math.min(width, Math.max(startPoint, event.pageX - boxX)) / width;
        if (max - min < 0.005) {
            // Do not zoom if the user only clicked once and did not move considerably in the drag event
            spectrumDragOut(event);
            return;
        }
        var minMz = Number(d.minMz);
        var maxMz = Number(d.maxMz);
        var maxIntensity = Number(d.maxIntensity);
        var min = min * Math.max(1, maxMz - minMz) + minMz;
        var max = max * Math.max(1, maxMz - minMz) + minMz;
        var offsetY = Math.min(Math.max(0, event.pageY - boxY), height);
        if (selection.classList.contains("second")) offsetY = height - offsetY;
        var maxI = (1 - Math.max(0, offsetY / height)) * maxIntensity;

        Zoom(canvas, min, max, maxI);
        if (linked_selection != undefined) Zoom(linked_selection.parentElement, min, max, maxI);
        spectrumDragOut(event)
        last_selection = undefined;
    }
}

var first_anchor;

function spectrumDistanceStart(event) {
    startPoint = undefined;
    first_anchor = event.target;
    while (!first_anchor.classList.contains("peak")) {
        first_anchor = first_anchor.parentElement;
    }
    event.preventDefault();
    document.querySelector("#spectrum-wrapper").classList.add("distance");
}

function spectrumDistanceEnd(event) {
    if (first_anchor != undefined) {
        let second_anchor = event.target;
        while (!second_anchor.classList.contains("peak")) {
            second_anchor = second_anchor.parentElement;
        }
        if (first_anchor == second_anchor) {
            document.querySelector("#spectrum-wrapper").classList.remove("distance");
            first_anchor = undefined;
            return;
        }

        const canvasWrapper = document.querySelector("#spectrum-wrapper .canvas-wrapper");
        const intensity = Math.min(Number(first_anchor.style.getPropertyValue("--intensity")), Number(second_anchor.style.getPropertyValue("--intensity"))) / 2;
        const start_mz = Math.min(Number(first_anchor.style.getPropertyValue("--mz")), Number(second_anchor.style.getPropertyValue("--mz")));
        const end_mz = Math.max(Number(first_anchor.style.getPropertyValue("--mz")), Number(second_anchor.style.getPropertyValue("--mz")));
        const text = distanceLabel(Number(canvasWrapper.dataset.maxMz), Number(canvasWrapper.dataset.minMz), end_mz - start_mz);

        let distance = document.createElement("span");
        distance.className = "distance";
        distance.style.setProperty("--start-mz", start_mz);
        distance.style.setProperty("--end-mz", end_mz);
        distance.style.setProperty("--intensity", intensity);
        distance.dataset.distance = end_mz - start_mz;
        distance.dataset.label = text;
        distance.addEventListener("click", event => event.target.remove());

        document.querySelector("#spectrum-wrapper .canvas").appendChild(distance);
        first_anchor = undefined;
        document.querySelector("#spectrum-wrapper").classList.remove("distance");
    }
}

export function spectrumClearDistanceLabels() {
    document.querySelectorAll("#spectrum-wrapper .distance").forEach(element => element.remove());
}

function spectrumScroll(event) {
    event.preventDefault(); // Prevent page scrolling (also when not on the canvas)
    let target = event.target;
    while (!target.classList.contains("canvas")) {
        target = target.parentElement;
    }
    const wrapper = target.parentElement;
    const minMz = Number(wrapper.dataset.minMz);
    const maxMz = Number(wrapper.dataset.maxMz);

    if (Math.abs(event.deltaY) != 0 && event.ctrlKey) { // Vertical scroll + control -> zoom in on Y axis
        let factor = 1 + 5 * event.deltaY / 10000; // Prevent scrolling past 0
        Zoom(wrapper, minMz, maxMz, Number(wrapper.dataset.maxIntensity) * factor);
    } else if (Math.abs(event.deltaX) != 0 || (Math.abs(event.deltaY) != 0 && event.shiftKey)) { // Horizontal scroll -> pan X axis
        let delta = Math.abs(event.deltaX) != 0 ? event.deltaX : event.deltaY;
        delta = Math.max(-minMz, delta / 10000 * (maxMz - minMz)); // Prevent scrolling past 0
        Zoom(wrapper, minMz + delta, maxMz + delta, Number(wrapper.dataset.maxIntensity));
    } else if (Math.abs(event.deltaY) != 0) { // Vertical scroll -> zoom in on X axis
        // Smaller thing (95%) centred to location
        let box = target.getBoundingClientRect();
        let center = (event.pageX - box.x) / box.width;
        let delta = -event.deltaY / 10000 * 5 * (maxMz - minMz);
        Zoom(wrapper, Math.max(0, minMz + delta * center), maxMz - delta * (1 - center), Number(wrapper.dataset.maxIntensity));
    }
}

let last_max_intensity = undefined;
function Zoom(canvas_wrapper, min, max, maxI) {
    canvas_wrapper.classList.add("zoomed");
    canvas_wrapper.dataset.minMz = min;
    canvas_wrapper.dataset.maxMz = max;
    canvas_wrapper.dataset.maxIntensity = maxI;
    canvas_wrapper.style.setProperty("--min-mz", min);
    canvas_wrapper.style.setProperty("--max-mz", max);
    canvas_wrapper.style.setProperty("--max-intensity", canvas_wrapper.dataset.maxIntensity);

    if (last_max_intensity != maxI) {
        SpectrumUpdateLabels(canvas_wrapper.querySelector(".canvas"));
    }

    UpdateSpectrumAxes(canvas_wrapper)

    document.getElementById("spectrum-mz-min").value = fancyRound(max, min, min);
    document.getElementById("spectrum-mz-max").value = fancyRound(max, min, max);
    document.getElementById("spectrum-intensity-max").value = fancyRound(Number(canvas_wrapper.dataset.maxIntensity), 0, Number(canvas_wrapper.dataset.maxIntensity));
    last_max_intensity = maxI;
}

function ZoomSpectrumGraph(canvas, min_y, max_y) {
    canvas.classList.add("zoomed");

    var spectrum_graph = canvas.parentElement.parentElement.parentElement;
    canvas.parentElement.style.setProperty("--y-min", fancyRound(max_y, min_y, min_y));
    canvas.parentElement.style.setProperty("--y-max", fancyRound(max_y, min_y, max_y));
    spectrum_graph.querySelector(".y-min").value = fancyRound(max_y, min_y, min_y);
    spectrum_graph.querySelector(".y-max").value = fancyRound(max_y, min_y, max_y);

    spectrum_graph.parentElement.querySelector('.error-graph-y-axis .min').innerText = fancyRound(max_y, min_y, min_y);
    spectrum_graph.parentElement.querySelector('.error-graph-y-axis .max').innerText = fancyRound(max_y, min_y, max_y);

    if (max_y > 0 && min_y > 0) spectrum_graph.querySelector(".x-axis").classList.add("hug-bottom");
    else spectrum_graph.querySelector(".x-axis").classList.remove("hug-bottom");
}

function spectrumZoomOut(event) {
    var spectrum = document.querySelector("#spectrum-wrapper");
    var min, max, maxI = 0;
    spectrum.querySelectorAll(".canvas-wrapper").forEach(canvas => {
        var d = canvas.dataset;
        d.minMz = 0;
        d.maxMz = d.initialMaxMz;
        d.maxIntensity = d.initialMaxIntensity;
        canvas.style.setProperty("--min-mz", 0);
        canvas.style.setProperty("--max-mz", d.initialMaxMz);
        canvas.style.setProperty("--max-intensity", d.initialMaxIntensity);
        canvas.classList.remove("zoomed");

        UpdateSpectrumAxes(canvas)
        min = Number(d.minMz);
        max = Number(d.maxMz);
        maxI = Number(d.maxIntensity);

        SpectrumUpdateLabels(canvas.querySelector(".canvas-spectrum"));
    });
    spectrum.querySelector(".mz-min").value = fancyRound(max, min, min);
    spectrum.querySelector(".mz-max").value = fancyRound(max, min, max);
    spectrum.querySelector(".intensity-max").value = fancyRound(maxI, 0, maxI);
}

function fancyRound(max, min, value, additional = 0) {
    const diff = max - min;
    const digits = Math.min(Math.max(Math.ceil(-Math.log10(diff) + 2 + additional), 0), 100);
    return value.toFixed(digits + additional);
}

function distanceLabel(max, min, value) {
    const diff = max - min;
    const zoom_digits = Math.max(Math.ceil(-Math.log10(diff) + 2), 0);
    const number_digits = Math.max(Math.ceil(-Math.log10(value) + 3), 0);
    return value.toFixed(Math.min(Math.max(zoom_digits, number_digits), 100));
}

// Give the canvas element
function UpdateSpectrumAxes(canvas_wrapper) {
    const spectrum_wrapper = canvas_wrapper.parentElement.parentElement;

    // Update x-axis
    const x_axis = canvas_wrapper.children[2]; // x-axis
    const x_ticks = x_axis.children;
    const x_ticks_number = Number(document.getElementById("spectrum-ticks-x").value);
    const x_ticks_delta = x_ticks.length - x_ticks_number;
    if (x_ticks_delta < 0) {
        for (let i = 0; i < Math.abs(x_ticks_delta); i++) {
            let new_tick = document.createElement("span");
            new_tick.classList.add("tick");
            x_axis.appendChild(new_tick);
        }
    } else if (x_ticks_delta > 0) {
        for (let i = 0; i < Math.abs(x_ticks_delta); i++) {
            x_axis.children[0].remove();
        }
    }
    const min_mz = Number(canvas_wrapper.dataset.minMz);
    const max_mz = Number(canvas_wrapper.dataset.maxMz);
    const diff_mz = max_mz - min_mz;
    const digits_mz = Math.min(Math.max(Math.ceil(-Math.log10(diff_mz) + 2), 0), 100);
    for (let i = 0; i < x_ticks.length; i++) {
        x_ticks[i].innerText = (min_mz + i / (x_ticks.length - 1) * (diff_mz)).toFixed(digits_mz);
    }

    // Update spectrum graph axes
    canvas_wrapper.querySelector('.error-graph .x-axis .min').innerText = min_mz.toFixed(digits_mz);
    canvas_wrapper.querySelector('.error-graph .x-axis .max').innerText = max_mz.toFixed(digits_mz);

    // Update y-axis
    const y_axis = canvas_wrapper.children[0]; // y-axis
    const y_ticks = y_axis.children; // y-axis
    const y_ticks_number = Number(document.getElementById("spectrum-ticks-y").value);
    const y_ticks_delta = y_ticks.length - y_ticks_number;
    if (y_ticks_delta < 0) {
        for (let i = 0; i < Math.abs(y_ticks_delta); i++) {
            let new_tick = document.createElement("span");
            new_tick.classList.add("tick");
            y_axis.appendChild(new_tick);
        }
    } else if (y_ticks_delta > 0) {
        for (let i = 0; i < y_ticks_delta; i++) {
            y_axis.children[0].remove();
        }
    }
    const sqrt = spectrum_wrapper.classList.contains("y-sqrt");
    const percent = spectrum_wrapper.classList.contains("y-percentage");
    const max_y = Number(canvas_wrapper.dataset.maxIntensity);
    const max_y_initial = Number(canvas_wrapper.dataset.initialMaxIntensity);

    for (let i = 0; i < y_ticks.length; i++) {
        let v = i / (y_ticks.length - 1) * max_y;
        if (sqrt) v = Math.sqrt(v);
        if (percent) v = v / max_y_initial * 100;
        if (sqrt && percent) v = Math.sqrt(i / (y_ticks.length - 1) * max_y) / Math.sqrt(max_y_initial) * 100;

        if (i == 0)
            y_ticks[i].innerText = "0";
        else if (percent)
            y_ticks[i].innerText = Math.round(v);
        else
            y_ticks[i].innerText = Math.round(v).toExponential(2);
    }

    // Update distance labels 
    canvas_wrapper.querySelectorAll(".distance").forEach(element => {
        const distance = Number(element.dataset.distance);
        const digits_number = Math.max(Math.ceil(-Math.log10(distance) + 3), 0)
        element.dataset.label = distance.toFixed(Math.max(digits_mz, digits_number));
    });
}