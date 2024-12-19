"use strict";

const { invoke } = window.__TAURI__.core;

const { listen } = window.__TAURI__.event;

const { open } = window.__TAURI__.dialog;

const RAW_EXTENSIONS = ['mgf', 'mgf.gz', "mzml", "mzml.gz", "imzml", "imzml.gz", "mzmlb", "mzmlb.gz", "raw", "raw.gz"];
const IDENTIFIED_EXTENSIONS = ["csv", "csv.gz", "tsv", "tsv.gz", "txt", "txt.gz", "psmtsv", "psmtsv.gz", "fasta", "fasta.gz", "mztab", "mztab.gz", "deepnovo_denovo", "deepnovo_denovo.gz", "ssl", "ssl.gz"];

import { SetUpSpectrumInterface, spectrumClearDistanceLabels } from "./script.js";

listen('tauri://drag-drop', event => {
  document.querySelector("html").classList.remove("file-drop-hover");
  let opens = [];
  let peptides = false;
  let raw_files = false;

  for (let i = 0; i < event.payload.paths.length; i++) {
    let file = event.payload.paths[i];
    let extension = file.toLowerCase().split('.').pop();
    if (extension == "gz") {
      extension = file.toLowerCase().split('.').reverse()[1] + ".gz";
    }
    if (RAW_EXTENSIONS.includes(extension)) {
      document.querySelector("#load-raw-path").classList.add("loading")
      raw_files = true;
      opens.push(load_raw(file));
    } else if (IDENTIFIED_EXTENSIONS.includes(extension)) {
      document.querySelector("#load-identified-peptides").classList.add("loading")
      peptides = true;
      opens.push(load_identified_peptides(file));
    } else {
      console.error("Extension not recognised", extension);
    }
  }
  Promise.all(opens).then(() => {
    if (peptides) {
      update_identified_peptide_file_select();
      identified_peptide_details();
      document.querySelector("#load-identified-peptides").classList.remove("loading");
    }
    if (raw_files) {
      update_open_raw_files();
      document.querySelector("#load-raw-path").classList.remove("loading");
    }
  });
})

listen('tauri://drag-enter', event => {
  document.querySelector("html").classList.add("file-drop-hover");
})

listen('tauri://drag-leave', event => {
  document.querySelector("html").classList.remove("file-drop-hover");
})

document.addEventListener("dragend", () => document.querySelector("html").classList.remove("file-drop-hover"));

// function abort() {
//   console.log(controller);
//   controller.abort("User manual abort");
// }

/**
* @param e: Element
*/
async function select_raw_file(e) {
  let properties = {
    directory: false,
    multiple: true,
    filters: [{
      extensions: RAW_EXTENSIONS, name: "*"
    }]
  };
  open(properties).then((result) => {
    if (result != null) {
      document.querySelector("#load-raw-path").classList.add("loading")
      let opens = [];
      for (let file of result) {
        opens.push(load_raw(file));
      }
      Promise.all(opens).then(() => {
        update_open_raw_files();
        document.querySelector("#load-raw-path").classList.remove("loading");
      });
    }
  })
};

/**
* @param e: Element
*/
async function dialog_select_identified_peptides_file(e) {
  let properties = {
    directory: false,
    multiple: true,
    filters: [{
      extensions: IDENTIFIED_EXTENSIONS, name: "*"
    }]
  };
  open(properties).then((result) => {
    if (result != null) {
      document.querySelector("#load-identified-peptides").classList.add("loading")
      let opens = [];
      for (let file of result) {
        opens.push(load_identified_peptides(file));
      }
      Promise.all(opens).then(() => {
        update_identified_peptide_file_select();
        identified_peptide_details();
        document.querySelector("#load-identified-peptides").classList.remove("loading");
      });
    }
  })
};

async function load_raw(path) {
  return invoke("load_raw", { path: path }).then(() => {
    clearError("open-files-error");
  }).catch((error) => {
    showError("open-files-error", error);
  });
}

async function load_usi() {
  document.querySelector("#load-usi").classList.add("loading")
  return invoke("load_usi", { usi: document.getElementById("usi").value }).then((result) => {
    clearError("open-files-error");
    document.getElementById("usi").value = "";
    document.querySelector("#peptide").innerText = result.peptide;
    document.querySelector("#spectrum-charge").value = result.charge;
    if (result.mode != null) {
      document.querySelector("#spectrum-model").value = result.mode.toLowerCase();
    }
    if (result.warning != null) {
      showError("spectrum-error", result.warning);
    } else {
      clearError("spectrum-error");
    }
    document.querySelector("#load-usi").classList.remove("loading")
    update_open_raw_files();
  }).catch((error) => {
    showError("open-files-error", error);
    document.querySelector("#load-usi").classList.remove("loading")
  })
}

async function load_identified_peptides(path) {
  return invoke("load_identified_peptides_file", { path: path }).then((result) => {
    if (result == null) {
      clearError("open-files-error");
    } else {
      showError("open-files-error", result);
    }
    displayed_identified_peptide = undefined;
  }).catch((error) => {
    showError("open-files-error", error);
  });
}

async function select_identified_peptides_file(id) {
  const select = document.getElementById("details-identified-peptide-files");
  for (let child of select.children) {
    if (child.value == id) {
      document.getElementById("details-identified-peptide-index").value = 0;
      document.getElementById("number-of-identified-peptides").innerText = child.dataset.length;
      child.selected = true;
      identified_peptide_details();
      return;
    }
  }
}

async function update_identified_peptide_file_select() {
  invoke("get_identified_peptides_files").then((result) => {
    const select = document.getElementById("details-identified-peptide-files");
    select.innerText = "";
    let id = 0;
    for (let row in result) {
      let option = document.createElement("option");
      option.value = result[row][0];
      option.innerText = "F" + (result[row][0] + 1) + ":" + result[row][1];
      option.title = result[row][2];
      option.dataset.length = result[row][3];
      if (row == result.length - 1) {
        id = result[row][0];
        option.selected = true;
      }
      select.appendChild(option);
    }
    select_identified_peptides_file(id);
    if (result.length == 0) {
      document.querySelector("#peptides").style.display = "hidden";
    } else {
      document.querySelector("#peptides").style.display = "block";
    }
  }).catch((error) => {
    showError("open-files-error", error);
    document.querySelector("#load-identified-peptides").classList.remove("loading")
  })
}

async function load_clipboard() {
  document.querySelector("#load-clipboard").classList.add("loading");
  return navigator.clipboard
    .readText()
    .then(async (clipText) => {
      invoke("load_clipboard", { data: clipText }).then(() => {
        clearError("open-files-error");
        update_open_raw_files();
        document.querySelector("#load-clipboard").classList.remove("loading")
      }).catch((error) => {
        showError("open-files-error", error);
        document.querySelector("#load-clipboard").classList.remove("loading")
      })
    })
    .catch(() => {
      showError("open-files-error", "Could not load clipboard, did you give permission to read the clipboard?");
      document.querySelector("#load-clipboard").classList.remove("loading")
    });
}

async function update_open_raw_files() {
  invoke("get_open_raw_files", {}).then((result) => {
    let root = document.querySelector("#spectra");
    root.innerHTML = '';

    if (result.length == 0) {
      root.innerText = "No raw files opened, either drag some in or open files with the buttons above."
      root.classList.add("hint");
    } else {
      root.classList.remove("hint");

      for (let file of result) {
        let rawfile = document.createElement("div");
        rawfile.className = "rawfile";
        let header = document.createElement("div");
        rawfile.appendChild(header);
        header.className = "header";
        header.id = "rawfile-" + file.id;
        header.appendChild(createElement("span", { text: "R" + (file.id + 1) + ":" + file.path.split('\\').pop().split('/').pop(), title: file.path }));

        if (file.single) {
          let select = document.createElement("button");
          select.innerText = file.selected ? "Deselect" : "Select";
          select.addEventListener("click", e => {
            let r = e.target.parentElement.parentElement;
            if (r.dataset.selected == "true") {
              console.log("deselect")
              invoke("deselect_spectrum", { fileIndex: file.id, index: 0 }).then(() => {
                r.dataset.selected = false;
                select.innerText = "Select";
                update_selected_spectra()
              });
            } else {
              console.log("select")
              invoke("select_spectrum_index", { fileIndex: file.id, index: 0 }).then(() => {
                r.dataset.selected = true;
                select.innerText = "Deselect";
                clearError("open-files-error");
                update_selected_spectra();
              }).catch((error) => {
                showError("open-files-error", error);
              });
            }
          });
          header.appendChild(select);
          rawfile.dataset.single = true;
          rawfile.dataset.selected = file.selected;
        } else {
          let spectrum_selection = document.createElement("div");
          spectrum_selection.className = "spectrum-selection";
          spectrum_selection.appendChild(createElement("label", {
            text: "Spectrum index",
            for: "spectrum-" + file.id + "-index",
            title: "The index of the spectrum is the 0 based index in the file, so the first spectrum has index 0, the second 1, etc"
          }))
          let group_index = document.createElement("div");
          group_index.className = "combined-input";
          let input_index = document.createElement("input");
          input_index.id = "spectrum-" + file.id + "-index";
          input_index.setAttribute("min", "0");
          input_index.setAttribute("type", "number");
          let select_on_index = () => {
            if (input_index.value.trim() != "") {
              invoke("select_spectrum_index", { fileIndex: file.id, index: Number(input_index.value) }).then(
                () => {
                  clearError("open-files-error");
                  input_index.value = "";
                  update_selected_spectra();
                }
              ).catch((error) => {
                showError("open-files-error", error);
              })
            }
          };
          input_index.addEventListener("focusout", select_on_index);
          input_index.addEventListener("keydown", event => { if (event.keyCode == 13) { select_on_index() } else { } });
          group_index.appendChild(input_index);
          group_index.appendChild(createElement("span", { text: "/" }))
          group_index.appendChild(createElement("span", { text: file.spectra, id: "spectrum-" + file.id + "-scans" }))
          spectrum_selection.appendChild(group_index);
          spectrum_selection.appendChild(createElement("label", {
            text: "Native ID",
            for: "spectrum-" + file.id + "-native-id",
            title: "The native ID is the identifier given to a spectrum in the original raw file, the format for this is different between different manufacturers. This can be used to track the same spectrum through multiple subsequent programs."
          }))
          let input_native_id = document.createElement("input");
          input_native_id.id = "spectrum-" + file.id + "-native-id";
          let select_on_scan = () => {
            if (input_native_id.value.trim() != "") {
              invoke("select_spectrum_native_id", { fileIndex: file.id, nativeId: input_native_id.value }).then(
                () => {
                  clearError("open-files-error");
                  input_native_id.value = "";
                  update_selected_spectra();
                }
              ).catch((error) => {
                showError("open-files-error", error);
              })
            }
          };
          input_native_id.addEventListener("focusout", select_on_scan);
          input_native_id.addEventListener("keydown", event => { if (event.keyCode == 13) { select_on_scan() } else { } });
          spectrum_selection.appendChild(input_native_id);
          header.appendChild(spectrum_selection);
          rawfile.dataset.single = false;
        }
        let close = document.createElement("button");
        close.innerText = "Close file";
        close.addEventListener("click", () => invoke("close_raw_file", { fileIndex: file.id }).then(update_open_raw_files()));
        header.appendChild(close);

        let spectra_list = document.createElement("ul");
        spectra_list.id = "rawfile-" + file.id + "-spectra";
        rawfile.appendChild(spectra_list);

        root.appendChild(rawfile);
      }

      // By definition always refresh the selected as well
      update_selected_spectra();
    }
  });
}

/// Refresh the selected spectra
async function update_selected_spectra() {
  invoke("get_selected_spectra", {}).then(
    (result) => {
      for (let i in result) {
        let index = Number(result[i][0]);
        let rawfile = document.getElementById("rawfile-" + index + "-spectra");
        rawfile.innerHTML = '';

        for (let element of result[i][1]) {
          let li = document.createElement("li");
          let single = rawfile.parentElement.dataset.single === "single";

          if (!single) {
            let close = document.createElement("button");
            close.innerText = "Deselect";
            close.addEventListener("click", () => {
              invoke("deselect_spectrum", { fileIndex: index, index: single ? 0 : element.id }).then(update_selected_spectra())
            });
            li.appendChild(close);
          }
          li.appendChild(createElement("span", { text: element.short }));
          let tooltip = document.createElement("div");
          tooltip.className = "tooltip";
          tooltip.innerHTML = element.description;
          li.appendChild(tooltip);
          li.setAttribute("data-index", String(element.id));

          rawfile.appendChild(li);
        }
      }
    }
  )
}

let displayed_identified_peptide = undefined;
async function identified_peptide_details() {
  let select = document.querySelector("#details-identified-peptide-files");
  if (select.children.length == 0) return;
  let file_id = Number(select.options[select.selectedIndex].value);
  let index = Number(document.querySelector("#details-identified-peptide-index").value);
  if (displayed_identified_peptide != [file_id, index]) {
    invoke("identified_peptide_details", { file: file_id, index: index }).then((result) => {
      document.querySelector("#identified-peptide-details").innerHTML = result;
      displayed_identified_peptide = [file_id, index];
      clearError("spectrum-error");
    }).catch((error) => {
      showError("spectrum-error", error);
      document.querySelector("#identified-peptide-details").innerText = "ERROR";
    })
  }
}

async function search_peptide() {
  if (document.querySelector("#search-peptide-input").value != "") {
    document.querySelector("#search-peptide").classList.add("loading");
    invoke("search_peptide", { text: document.querySelector("#search-peptide-input").value, minimalMatchScore: Number(document.querySelector("#search-peptide-minimal-match").value), minimalPeptideScore: Number(document.querySelector("#search-peptide-minimal-peptide").value) }).then((result) => {
      document.querySelector("#resulting-peptides").innerHTML = result;
      document.querySelector("#search-peptide").classList.remove("loading");
      clearError("spectrum-error");
    }).catch((error) => {
      showError("spectrum-error", error);
      document.querySelector("#search-peptide").classList.remove("loading");
      document.querySelector("#resulting-peptides").innerText = "ERROR";
    })
  } else {
    document.querySelector("#resulting-peptides").innerText = "";
  }
}

async function close_identified_peptide_file() {
  let select = document.querySelector("#details-identified-peptide-files");
  invoke("close_identified_peptides_file", { file: Number(select.options[select.selectedIndex].value) }).then(() => {
    clearError("spectrum-error");
    if (select.children.length <= 1) {
      document.querySelector("#peptides").style.display = "none";
    } else {
      update_identified_peptide_file_select();
      identified_peptide_details();
    }
  }).catch((error) => {
    showError("spectrum-error", error);
    document.querySelector("#identified-peptide-details").innerText = "ERROR";
  })
}

export async function load_peptide(file, peptide) {
  select_identified_peptides_file(file);
  document.querySelector("#details-identified-peptide-index").value = peptide;
  identified_peptide_details();
}
window.load_peptide = load_peptide;

async function search_modification() {
  if (document.querySelector("#search-modification").value != "") {
    document.querySelector("#search-modification-button").classList.add("loading");
    invoke("search_modification", { text: document.querySelector("#search-modification").value, tolerance: Number(document.querySelector("#search-modification-tolerance").value) }).then((result) => {
      document.querySelector("#search-modification-result").innerHTML = result;
      document.querySelector("#search-modification-button").classList.remove("loading");
      clearError("search-modification-error");
    }).catch((error) => {
      showError("search-modification-error", error);
      document.querySelector("#search-modification-button").classList.remove("loading");
    })
  } else {
    document.querySelector("#resulting-modification-result").innerText = "";
  }
}

async function details_formula(event) {
  const formula = document.getElementById("details-formula");
  invoke("details_formula", { text: formula.innerText }).then((result) => {
    document.querySelector("#details-formula-result").innerHTML = result;
    if (formula.querySelector(".error") != null) {
      formula.innerText = formula.innerText;
      moveCursorToEnd(formula);
    }
    clearError("details-formula-error");
  }).catch((error) => {
    showError("details-formula-error", error, false);
    formula.innerHTML = showContext(error, formula.innerText);
    moveCursorToEnd(formula);
  })
}

async function load_identified_peptide() {
  let select = document.querySelector("#details-identified-peptide-files");
  if (select.children.length == 0) return;
  let file_id = Number(select.options[select.selectedIndex].value);
  let index = Number(document.querySelector("#details-identified-peptide-index").value);
  document.getElementById("load-identified-peptide").classList.add("loading");
  invoke("load_identified_peptide", { index: index, file: file_id }).then((result) => {
    document.querySelector("#peptide").innerText = result.peptide;
    document.querySelector("#spectrum-charge").value = result.charge;
    if (result.mode != null) {
      document.querySelector("#spectrum-model").value = result.mode.toLowerCase();
    }
    if (result.warning != null) {
      showError("spectrum-error", result.warning);
    } else {
      clearError("spectrum-error");
    }
    update_selected_spectra();
    document.getElementById("load-identified-peptide").classList.remove("loading");
  }).catch(() => {
    document.querySelector("#identified-peptide-details").innerText = "ERROR";
    document.getElementById("load-identified-peptide").classList.remove("loading");
  })
}


function get_location(id) {
  let loc = document.querySelector(id);
  let t = loc.children[0].options[Number(loc.children[0].value)].dataset.value;
  if (["SkipN", "SkipC", "TakeC"].includes(t)) {
    return { [t]: Number(loc.children[1].value) };
  } else if (t == "TakeN1") {
    return { ["TakeN"]: { "skip": 0, "take": Number(loc.children[1].value) } };
  } else if (t == "SkipNC") {
    return { [t]: [Number(loc.children[1].value), Number(loc.children[2].value)] };
  } else if (t == "TakeN") {
    return { [t]: { "skip": Number(loc.children[1].value), "take": Number(loc.children[2].value) } };
  } else {
    return t;
  }
}

function get_losses(ion) {
  let loss = loadSeparatedInput("model-" + ion + "-loss");
  document.getElementsByName("model-" + ion + "-loss-selection").forEach(element => {
    if (element.checked) {
      loss.push(element.value);
    }
  });
  return loss;
}

function get_charge_range(ion) {
  let start_type = document.getElementById("model-" + ion + "-charge-start-type").value;
  let start_value = Number(document.getElementById("model-" + ion + "-charge-start-value").value);
  let end_type = document.getElementById("model-" + ion + "-charge-end-type").value;
  let end_value = Number(document.getElementById("model-" + ion + "-charge-end-value").value);
  return { start: { [start_type]: start_value }, end: { [end_type]: end_value } };
}

async function annotate_spectrum() {
  document.querySelector("#annotate-button").classList.add("loading");
  document.querySelector("#peptide").innerText = document.querySelector("#peptide").innerText.trim();
  var charge = document.querySelector("#spectrum-charge").value == "" ? null : Number(document.querySelector("#spectrum-charge").value);
  var noise_threshold = Number(document.querySelector("#noise-filter").value);
  var model = {
    a: [get_location("#model-a-location"), get_losses("a"), get_charge_range("a")],
    b: [get_location("#model-b-location"), get_losses("b"), get_charge_range("b")],
    c: [get_location("#model-c-location"), get_losses("c"), get_charge_range("c")],
    d: [get_location("#model-d-location"), get_losses("d"), get_charge_range("d")],
    v: [get_location("#model-v-location"), get_losses("v"), get_charge_range("v")],
    w: [get_location("#model-w-location"), get_losses("w"), get_charge_range("w")],
    x: [get_location("#model-x-location"), get_losses("x"), get_charge_range("x")],
    y: [get_location("#model-y-location"), get_losses("y"), get_charge_range("y")],
    z: [get_location("#model-z-location"), get_losses("z"), get_charge_range("z")],
    precursor: [get_losses("precursor"), get_charge_range("precursor")],
    immonium: [document.querySelector("#model-immonium-enabled").checked, get_charge_range("immonium")],
    m: document.querySelector("#model-m-enabled").checked,
    modification_neutral: document.querySelector("#model-modification-neutral-enabled").checked,
    modification_diagnostic: [document.querySelector("#model-modification-diagnostic-enabled").checked, get_charge_range("diagnostic")],
    cleave_cross_links: document.querySelector("#model-cleave-cross-links-enabled").checked,
    glycan: [
      document.querySelector("#model-glycan-enabled").checked,
      [Number(document.querySelector("#model-glycan-composition-min").value), Number(document.querySelector("#model-glycan-composition-max").value)],
      get_losses("glycan"),
      get_charge_range("glycan-oxonium"),
      get_charge_range("glycan-other")],
  };
  invoke("annotate_spectrum", {
    tolerance: [Number(document.querySelector("#spectrum-tolerance").value), document.querySelector("#spectrum-tolerance-unit").value],
    charge: charge,
    filter: noise_threshold,
    model: document.querySelector("#spectrum-model").value,
    peptide: document.querySelector("#peptide").innerText, customModel: model,
    massMode: document.querySelector("#spectrum-mass-mode").value,
    mzRange: [optional_number(document.querySelector("#model-mz-range-min").value), optional_number(document.querySelector("#model-mz-range-max").value)]
  }).then((result) => {
    document.querySelector("#spectrum-results-wrapper").innerHTML = result.spectrum;
    document.querySelector("#spectrum-fragment-table").innerHTML = result.fragment_table;
    document.querySelector("#spectrum-error").innerText = "";
    document.querySelector("#spectrum-wrapper").classList.remove("hidden"); // Remove hidden class if this is the first run
    document.querySelector("#spectrum-mz-max").value = result.mz_max;
    document.querySelector("#spectrum-intensity-max").value = result.intensity_max;
    document.querySelector("#spectrum-label").value = 90;
    document.querySelector("#spectrum-label-value").value = 90;
    document.querySelector("#spectrum-m-z").value = 0;
    document.querySelector("#spectrum-m-z-value").value = 0;
    SetUpSpectrumInterface();
    document.querySelector("#annotate-button").classList.remove("loading");
    clearError("spectrum-error");
  }).catch((error) => {
    showError("spectrum-error", error, false);
    document.querySelector("#peptide").innerHTML = showContext(error, document.querySelector("#peptide").innerText);
    document.querySelector("#annotate-button").classList.remove("loading");
  })
}

function formatError(error, showContext = true) {
  console.error(error);
  if (typeof error == "string") {
    return "<div class='raw'>" + error + "</div>";
  } else {
    let msg = "<p class='title'>" + error.content.short_description + "</p><p class='description'>" + error.content.long_description + "</p>";
    if (error.content.version != "") {
      msg += "<p class='version'>Version: " + error.content.version + "</p>";
    }
    if (showContext) {
      if (error.content.context.hasOwnProperty('Line')) {
        let Line = error.content.context.Line;
        msg += "<div class='context'>" + (Line.line_index != null ? ("<span class='line-number'>" + (Line.line_index + 1) + "</span>") : "") + "<pre>" + Line.line + "\n" + " ".repeat(Line.offset) + "^".repeat(Line.length) + "</pre></div>";
      } else if (error.content.context.hasOwnProperty('Show')) {
        msg += "<div class='error'>" + error.content.context.Show.line + "</div>";
      } else if (error.content.context.hasOwnProperty('FullLine')) {
        let FullLine = error.content.context.FullLine;
        msg += "<div class='error'>" + FullLine.line + "</div>";
      } else if (error.content.context == "None") {
        // Empty
      } else {
        msg += "<pre>" + error.content.context + "</pre>";
      }
    }
    if (error.content.suggestions.length > 0) {
      msg += "<p>Did you mean any of the following?</p><ul>";
      for (let suggestion in error.content.suggestions) {
        msg += "<li>" + error.content.suggestions[suggestion] + "</li>";
      }
      msg += "</ul>";
    }
    if (error.content.underlying_errors.length > 0) {
      msg += "<label><input type='checkbox'></input>Show " + String(error.content.underlying_errors.length) + " underlying errors</label><ul>";
      for (let underlying_error in error.content.underlying_errors) {
        msg += "<li class='underlying-error'>" + formatError(error.content.underlying_errors[underlying_error], showContext) + "</li>";
      }
      msg += "</ul>";
    }
    return msg;
  }
}

/**
 * @param {Object} error - The error object
 * @param {String} fallback - The fallback full original text if the error is unsupported or None
 * @returns {String} String representation of HTML for use in `element.innerHTML = result;`
 */
function showContext(error, fallback) {
  if (error.content.context.hasOwnProperty('Line')) {
    let Line = error.content.context.Line;
    return Line.line.slice(0, Line.offset) + "<span class='error'>" + Line.line.slice(Line.offset, Line.offset + Line.length) + "</span>" + Line.line.slice(Line.offset + Line.length, Line.line.length);
  } else if (error.content.context.hasOwnProperty('FullLine')) {
    let FullLine = error.content.context.FullLine;
    return "<span class='error'>" + FullLine.line + "</span>";
  } else if (error.content.context.hasOwnProperty('Show')) {
    return "<span class='error'>" + error.content.context.Show.line + "</span>";
  } else if (error.content.context = "None") {
    return fallback;
  } else {
    console.error("Error type not handled", error);
    return fallback;
  }
}

/** @param e {MouseEvent}  */
function resizeDown(e) {
  document.querySelector(".resize-wrapper").classList.add("active");
  document.querySelector(".resize-wrapper").dataset.start_x = e.clientX;
  document.querySelector(".resize-wrapper").dataset.left_width = document.querySelector(".resize-wrapper").firstElementChild.getBoundingClientRect().width - 16;
  document.addEventListener("mousemove", resizeMove);
  document.addEventListener("mouseup", resizeUp);
}

/** @param e {MouseEvent}  */
function resizeMove(e) {
  let wrapper = document.querySelector(".resize-wrapper");
  let first = wrapper.firstElementChild;
  first.style.width =
    Math.max(10, Math.min(90, (Number(wrapper.dataset.left_width) + (e.clientX - Number(wrapper.dataset.start_x))) /
      wrapper.getBoundingClientRect().width * 100)) + "%";
  e.preventDefault();
}

function resizeUp() {
  document.querySelector(".resize-wrapper").classList.remove("active");
  document.removeEventListener("mousemove", resizeMove);
  document.removeEventListener("mouseup", resizeUp);
}

// Setup
window.addEventListener("DOMContentLoaded", () => {
  document.querySelector(".resize").addEventListener("mousedown", resizeDown);
  document.querySelectorAll(".collapsible>legend").forEach(element => element.addEventListener("click", (e) => document.getElementById(e.target.parentElement.dataset.linkedItem).toggleAttribute("checked")));
  document
    .querySelector("#load-raw-path")
    .addEventListener("click", (event) => select_raw_file(event.target));
  document
    .querySelector("#load-usi")
    .addEventListener("click", () => load_usi());
  enter_event("#usi", load_usi)
  document
    .querySelector("#load-identified-peptides")
    .addEventListener("click", (event) => dialog_select_identified_peptides_file(event.target));
  document
    .querySelector("#details-identified-peptide-files")
    .addEventListener("change", (event) => {
      select_identified_peptides_file(Number(event.target.options[event.target.selectedIndex].value))
    });
  document
    .querySelector("#close-identified-peptide-file")
    .addEventListener("click", (event) => close_identified_peptide_file());
  document
    .querySelector("#load-identified-peptide")
    .addEventListener("click", (event) => load_identified_peptide());
  document
    .querySelector("#load-clipboard")
    .addEventListener("click", () => load_clipboard());
  document
    .querySelector("#search-peptide", search_peptide)
    .addEventListener("click", () => search_peptide());
  document
    .querySelector("#search-modification-button")
    .addEventListener("click", () => search_modification());
  document
    .querySelector("#details-formula")
    .addEventListener("input", details_formula);
  enter_event("#search-peptide-input", search_peptide)
  enter_event("#search-modification", search_modification)
  add_event("#details-identified-peptide-index", ["change", "focus"], identified_peptide_details)
  document
    .querySelector("#annotate-button")
    .addEventListener("click", () => annotate_spectrum());
  document
    .querySelector("#theme")
    .addEventListener("change", (e) => {
      document.body.classList.remove("theme-light", "theme-dark", "theme-auto");
      document.body.classList.add(e.target.value);
    });
  document
    .querySelector("#distance-labels-clear")
    .addEventListener("click", () => spectrumClearDistanceLabels());
  document
    .querySelector("#peptide")
    .addEventListener("focus", (event) => {
      event.target.innerHTML = event.target.innerText;
    });
  enter_event("#peptide", annotate_spectrum)

  // Set up all separated inputs
  document.querySelectorAll(".separated-input").forEach(t => {
    t.addEventListener("click", event => {
      if (event.target.classList.contains("separated-input")) {
        event.target.querySelector(".input").focus({ focusVisible: true })
        event.preventDefault();
      }
    });
  });
  document.querySelectorAll(".separated-input .input").forEach(t => {
    t.addEventListener("keydown", async event => {
      let input = event.target;
      let values_container = input.parentElement;
      let outer = input.parentElement.parentElement;
      outer.classList.toggle("error", false);
      if ((event.keyCode == 13 || event.keyCode == 9) && input.innerText.trim() != "") { // Enter or Tab
        addValueSeparatedElement(values_container, input.innerText);
        event.preventDefault();
      } else if (event.keyCode == 8 && !input.hasChildNodes()) { // Backspace
        if (values_container.children.length >= 3) {
          let target = values_container.children[values_container.children.length - 3];
          input.innerText = target.innerText.slice(0, -1);
          target.remove();
          moveCursorToEnd(input);
          event.preventDefault();
        }
      }
    });
    t.addEventListener("focusout", async event => {
      let input = event.target;
      if (input.innerText.trim() != "") {
        input.parentElement.parentElement.classList.toggle("error", false);
        addValueSeparatedElement(input.parentElement, input.innerText);
      }
    });
  });
  document.querySelectorAll(".separated-input .clear").forEach(t => {
    t.addEventListener("click", e => {
      clearSeparatedInput(e.target.parentElement.parentElement);
    })
  });

  // Custom mods
  document.getElementById("custom-mod-create").addEventListener("click", e => {
    loadCustomModification();
    document.getElementById("custom-mod-id").value = e.target.dataset.newId;
    document.getElementById("custom-mod-example-id").innerText = "CUSTOM:" + e.target.dataset.newId;
    document.getElementById("custom-mod-dialog").showModal()
  });
  document.getElementById("custom-mod-formula").addEventListener("focusin", e => e.target.innerText = e.target.innerText);
  document.getElementById("custom-mod-formula").addEventListener("focusout", async e => {
    e.target.parentElement.classList.remove("error");
    if (e.target.innerText.trim() != "") {
      let text = e.target.innerText;
      invoke("validate_molecular_formula", { text: text })
        .then(value => e.target.innerHTML = value)
        .catch(error => {
          e.target.parentElement.parentElement.classList.add("error");
          e.target.parentElement.parentElement.querySelector("output.error").innerHTML = formatError(error, false);
          e.target.innerHTML = showContext(error, text);
        });
    }
  });
  document.getElementById("custom-mod-name").addEventListener("input", e => document.getElementById("custom-mod-example-name").innerText = "C:" + e.target.value);
  document.querySelectorAll(".list-input").forEach(t => {
    t.querySelector(".create").addEventListener("click", e => {
      e.target.parentElement.classList.add("creating");
    })
    t.querySelector(".save").addEventListener("click", async e => {
      let listInput = e.target.parentElement.parentElement;
      listInput.classList.remove("creating");
      let new_element = document.createElement("li");
      new_element.classList.add("element");
      if (listInput.classList.contains("single")) {
        new_element.innerHTML = await invoke("validate_custom_single_specificity", {
          placementRules: loadSeparatedInput("custom-mod-single-placement-rules"),
          neutralLosses: loadSeparatedInput("custom-mod-single-neutral-losses"),
          diagnosticIons: loadSeparatedInput("custom-mod-single-diagnostic-ions")
        }).catch(error => {
          console.error(error)
        });
      } else if (listInput.classList.contains("linker")) {
        new_element.innerHTML = await invoke("validate_custom_linker_specificity", {
          asymmetric: document.getElementById("custom-mod-linker-asymmetric").checked,
          placementRules: loadSeparatedInput("custom-mod-linker-placement-rules"),
          secondaryPlacementRules: loadSeparatedInput("custom-mod-linker-secondary-placement-rules"),
          stubs: loadSeparatedInput("custom-mod-linker-stubs"),
          diagnosticIons: loadSeparatedInput("custom-mod-linker-diagnostic-ions")
        }).catch(error => {
          console.error(error)
        });
      }
      new_element.children[0].title = "Edit";
      new_element.children[0].addEventListener("click", e => {
        let data = JSON.parse(e.target.dataset.value);
        if (listInput.classList.contains("single")) {
          populateSeparatedInput("custom-mod-single-placement-rules", data.placement_rules);
          populateSeparatedInput("custom-mod-single-neutral-losses", data.neutral_losses);
          populateSeparatedInput("custom-mod-single-diagnostic-ions", data.diagnostic_ions);
        } else if (listInput.classList.contains("linker")) {
          document.getElementById("custom-mod-linker-asymmetric").checked = data.asymmetric;
          populateSeparatedInput("custom-mod-linker-placement-rules", data.placement_rules);
          populateSeparatedInput("custom-mod-linker-secondary-placement-rules", data.secondary_placement_rules);
          populateSeparatedInput("custom-mod-linker-stubs", data.stubs);
          populateSeparatedInput("custom-mod-linker-diagnostic-ions", data.diagnostic_ions);
        }
        listInput.classList.add("creating");
        e.target.parentElement.classList.add("hidden");
      });
      let delete_button = document.createElement("button");
      delete_button.classList.add("delete");
      delete_button.appendChild(document.createTextNode("x"));
      delete_button.addEventListener("click", e => e.target.parentElement.remove());
      delete_button.title = "Delete";
      new_element.appendChild(delete_button);
      listInput.querySelector(".values").appendChild(new_element);
      listInput.querySelectorAll(".modal .separated-input").forEach(s => clearSeparatedInput(s));
    })
    t.querySelector(".cancel").addEventListener("click", e => {
      e.target.parentElement.parentElement.classList.remove("creating");
      e.target.parentElement.parentElement.querySelectorAll(".hidden").forEach(e => e.classList.remove("hidden"));
      e.target.parentElement.parentElement.querySelectorAll(".modal .separated-input").forEach(s => clearSeparatedInput(s));
    })
  });
  document.getElementById("custom-mod-save").addEventListener("click", e => {
    e.target.parentElement.parentElement.querySelectorAll(".hidden").forEach(e => e.remove());
    document.getElementById("custom-mod-dialog").close();
    let mod = {
      id: Number(document.getElementById("custom-mod-id").value),
      name: document.getElementById("custom-mod-name").value,
      formula: document.getElementById("custom-mod-formula").innerText,
      description: document.getElementById("custom-mod-description").value,
      synonyms: loadSeparatedInput("custom-mod-synonyms"),
      cross_ids: loadSeparatedInput("custom-mod-cross-ids"),
      linker: document.getElementById("custom-mod-type-linker").checked,
      single_specificities: loadListInput("custom-mod-single-specificities"),
      linker_specificities: loadListInput("custom-mod-linker-specificities"),
      linker_length: Number(document.getElementById("custom-mod-linker-length").value),
    };
    invoke("update_modification", {
      customModification: mod
    })
      .then(() => updateCustomModifications())
      .catch(error => console.error(error))
  });
  document.getElementById("custom-mod-cancel").addEventListener("click", () => document.getElementById("custom-mod-dialog").close());
  invoke("get_custom_mods_path").then(path =>
    document.getElementById("custom-modifications-path").innerText = path
  );
  updateCustomModifications();

  // Refresh interface for hot reload
  invoke("refresh").then((result) => {
    if (result[0] > 0) {
      update_selected_spectra();
    }
    document.querySelector("#number-of-identified-peptides").innerText = result[1];
    if (result[1] > 0) {
      document.querySelector("#peptides").style.display = "block";
      identified_peptide_details();
    }
    update_identified_peptide_file_select();
    update_open_raw_files();
  })
});

function updateCustomModifications() {
  invoke("get_custom_modifications")
    .then(modifications => {
      let container = document.getElementById("custom-mods");
      container.innerText = "";
      let highest_id = -1;
      for (let modification of modifications) {
        let new_element = document.createElement("li");
        new_element.dataset.id = modification[0];
        new_element.innerHTML = modification[1];
        let edit_button = document.createElement("button");
        edit_button.classList.add("edit");
        edit_button.appendChild(document.createTextNode("Edit"));
        edit_button.addEventListener("click", () =>
          invoke("get_custom_modification", { id: modification[0] })
            .then(result => {
              loadCustomModification(result);
              document.getElementById("custom-mod-dialog").showModal();
            })
            .catch(error => console.error(error)));
        new_element.appendChild(edit_button);
        let duplicate_button = document.createElement("button");
        duplicate_button.classList.add("duplicate");
        duplicate_button.appendChild(document.createTextNode("Duplicate"));
        duplicate_button.addEventListener("click", () =>
          invoke("duplicate_custom_modification", { id: modification[0], newId: Number(document.getElementById("custom-mod-create").dataset.newId) })
            .then(result => {
              loadCustomModification(result);
              document.getElementById("custom-mod-dialog").showModal();
            })
            .catch(error => console.error(error)));
        new_element.appendChild(duplicate_button);
        let delete_button = document.createElement("button");
        delete_button.classList.add("delete");
        delete_button.classList.add("secondary");
        delete_button.appendChild(document.createTextNode("Delete"));
        delete_button.addEventListener("click", () =>
          invoke("delete_custom_modification", { id: modification[0] })
            .then(() => {
              updateCustomModifications();
            })
            .catch(error => console.error(error)));
        new_element.appendChild(delete_button);
        container.appendChild(new_element);
        highest_id = Math.max(highest_id, modification[0]);
      }
      document.getElementById("custom-mod-create").dataset.newId = highest_id + 1;
    })
    .catch(error => console.error(error))
}

/**
 * @param {Object?} modification - If null clear
 */
function loadCustomModification(modification = null) {
  if (modification == null) {
    document.getElementById("custom-mod-id").value = 0;
    document.getElementById("custom-mod-example-id").innerText = "CUSTOM:0";
    document.getElementById("custom-mod-name").value = "";
    document.getElementById("custom-mod-example-name").innerText = "C:NAME";
    document.getElementById("custom-mod-formula").innerText = "";
    document.getElementById("custom-mod-description").value = "";
    clearSeparatedInput(document.getElementById("custom-mod-synonyms").parentElement);
    clearSeparatedInput(document.getElementById("custom-mod-cross-ids").parentElement);
    document.getElementById("custom-mod-type-single").checked = true;
    document.getElementById("custom-mod-type-linker").checked = false;
    document.getElementById("custom-mod-single-specificities").innerText = "";
    document.getElementById("custom-mod-linker-length").value = "";
    document.getElementById("custom-mod-linker-specificities").innerText = "";
  } else {
    document.getElementById("custom-mod-id").value = modification.id;
    document.getElementById("custom-mod-example-id").innerText = "CUSTOM:" + modification.id;
    document.getElementById("custom-mod-name").value = modification.name;
    document.getElementById("custom-mod-example-name").innerText = "C:" + modification.name;
    document.getElementById("custom-mod-formula").innerText = modification.formula;
    document.getElementById("custom-mod-description").value = modification.description;
    populateSeparatedInput("custom-mod-synonyms", modification.synonyms);
    populateSeparatedInput("custom-mod-cross-ids", modification.cross_ids);
    document.getElementById("custom-mod-type-single").checked = !modification.linker;
    document.getElementById("custom-mod-type-linker").checked = modification.linker;
    document.getElementById("custom-mod-single-specificities").innerText = "";
    for (let specificity of modification.single_specificities) {
      addValueListInput(
        document.getElementById("custom-mod-single-specificities").parentElement,
        specificity[0],
        specificity[1],
        specificity[2],
        false, [], [], [], []
      )
    }
    document.getElementById("custom-mod-linker-length").value = modification.linker_length;
    document.getElementById("custom-mod-linker-specificities").innerText = "";
    for (let specificity of modification.linker_specificities) {
      addValueListInput(
        document.getElementById("custom-mod-linker-specificities").parentElement,
        [], [], [],
        specificity[0],
        specificity[1],
        specificity[2],
        specificity[3],
        specificity[4],
      )
    }
  }
}

function loadListInput(id) {
  let listInput = document.getElementById(id).parentElement;
  let values = [...document.getElementById(id).querySelectorAll(".element>span")].map(e => JSON.parse(e.dataset.value));
  if (listInput.classList.contains("single")) {
    return values.map(v => [v.placement_rules, v.neutral_losses, v.diagnostic_ions]);
  } else if (listInput.classList.contains("linker")) {
    return values.map(v => [v.asymmetric, v.placement_rules, v.secondary_placement_rules, v.stubs, v.diagnostic_ions]);
  }
}

async function addValueListInput(listInput, singlePlacementRules, singleNeutralLosses, singleDiagnosticIons, linkerAsymmetric, linkerPlacementRules, linkerSecondaryPlacementRules, linkerStubs, linkerDiagnosticIons) {
  listInput.classList.remove("creating");
  let new_element = document.createElement("li");
  new_element.classList.add("element");
  if (listInput.classList.contains("single")) {
    new_element.innerHTML = await invoke("validate_custom_single_specificity", {
      placementRules: singlePlacementRules,
      neutralLosses: singleNeutralLosses,
      diagnosticIons: singleDiagnosticIons,
    }).catch(error => {
      console.error(error)
    });
  } else if (listInput.classList.contains("linker")) {
    new_element.innerHTML = await invoke("validate_custom_linker_specificity", {
      asymmetric: linkerAsymmetric,
      placementRules: linkerPlacementRules,
      secondaryPlacementRules: linkerSecondaryPlacementRules,
      stubs: linkerStubs,
      diagnosticIons: linkerDiagnosticIons,
    }).catch(error => {
      console.error(error)
    });
  }
  new_element.children[0].title = "Edit";
  new_element.children[0].addEventListener("click", e => {
    let data = JSON.parse(e.target.dataset.value);
    if (listInput.classList.contains("single")) {
      populateSeparatedInput("custom-mod-single-placement-rules", data.placement_rules);
      populateSeparatedInput("custom-mod-single-neutral-losses", data.neutral_losses);
      populateSeparatedInput("custom-mod-single-diagnostic-ions", data.diagnostic_ions);
    } else if (listInput.classList.contains("linker")) {
      document.getElementById("custom-mod-linker-asymmetric").checked = data.asymmetric;
      populateSeparatedInput("custom-mod-linker-placement-rules", data.placement_rules);
      populateSeparatedInput("custom-mod-linker-secondary-placement-rules", data.secondary_placement_rules);
      populateSeparatedInput("custom-mod-linker-stubs", data.stubs);
      populateSeparatedInput("custom-mod-linker-diagnostic-ions", data.diagnostic_ions);
    }
    listInput.classList.add("creating");
    e.target.parentElement.classList.add("hidden");
  });
  let delete_button = document.createElement("button");
  delete_button.classList.add("delete");
  delete_button.appendChild(document.createTextNode("x"));
  delete_button.addEventListener("click", e => e.target.parentElement.remove());
  delete_button.title = "Delete";
  new_element.appendChild(delete_button);
  listInput.querySelector(".values").appendChild(new_element);
}

function moveCursorToEnd(contentEle) {
  const range = document.createRange();
  const selection = window.getSelection();
  range.setStart(contentEle, contentEle.childNodes.length);
  range.collapse(true);
  selection.removeAllRanges();
  selection.addRange(range);
};

/** 
 * @param {String} id - The ID of the `.separated-input .values` element. 
 * @returns {String[]} - List of all elements, each of those as a string.
*/
function loadSeparatedInput(id) {
  return [...document.getElementById(id).querySelectorAll(".element")].map(c => { return c.dataset.value; });
}

/**
 * @param {Element} element - The outer `.separate-input` element.
 */
function clearSeparatedInput(element) {
  element.querySelectorAll(".element").forEach(c => c.remove());
  element.querySelector(".input").innerText = "";
  element.classList.remove("error");
  element.querySelector(".error").innerText = "";
}

/**
 * @param {String} id - The ID of the `.separated-input .values` element.
 * @param {String[]} values - The separate values to populate the field with.
 */
function populateSeparatedInput(id, values) {
  let element = document.getElementById(id);
  clearSeparatedInput(element.parentElement);
  for (let value of values) {
    addValueSeparatedElement(element, value);
  }
}

/**
 * @param {Element} element - The `.separated-input .values` element.
 * @param {String} value - The value to add
 */
async function addValueSeparatedElement(element, value) {
  let input = element.querySelector(".input");
  let outer = element.parentElement;
  let verified_value = undefined;
  switch (input.dataset.type) {
    case "molecular_formula":
      verified_value = await invoke("validate_molecular_formula", { text: value })
        .catch(error => {
          input.innerHTML = showContext(error, value);
          outer.querySelector("output.error").innerHTML = formatError(error, false);
        });
      break;
    case "neutral_loss":
      verified_value = await invoke("validate_neutral_loss", { text: value })
        .catch(error => {
          input.innerHTML = showContext(error, value);
          outer.querySelector("output.error").innerHTML = formatError(error, false);
        });
      break;
    case "placement_rule":
      verified_value = await invoke("validate_placement_rule", { text: value })
        .catch(error => {
          input.innerHTML = showContext(error, value);
          outer.querySelector("output.error").innerHTML = formatError(error, false);
        });
      break;
    case "stub":
      verified_value = await invoke("validate_stub", { text: value })
        .catch(error => {
          input.innerHTML = showContext(error, value);
          outer.querySelector("output.error").innerHTML = formatError(error, false);
        });
      break;
    case "cross_id":
      verified_value = value.trim();
      if (!value.includes(':')) {
        outer.querySelector("output.error").innerText = "Cross ID should contain a colon ':'";
        value = undefined;
      }
      break;
    case "text":
      verified_value = value.trim();
      break;
    default: console.error("Invalid separated input type");
  }
  if (verified_value !== undefined) {
    let new_element = document.createElement("span");
    new_element.classList.add("element");
    new_element.innerHTML = verified_value;
    new_element.dataset.value = new_element.innerText;
    new_element.addEventListener("click", e => {
      let input = e.target.parentElement.querySelector(".input");
      input.innerText = e.target.innerText.slice(0, -1);
      moveCursorToEnd(input);
      e.target.remove();
    });
    let delete_button = document.createElement("button");
    delete_button.classList.add("delete");
    delete_button.appendChild(document.createTextNode("x"));
    delete_button.addEventListener("click", e => e.target.parentElement.remove());
    delete_button.title = "Delete";
    new_element.appendChild(delete_button);

    element.insertBefore(new_element, input);
    input.innerText = "";
  } else {
    outer.classList.toggle("error", true);
    moveCursorToEnd(input);
  }
}

function add_event(selector, events, callback) {
  for (let i = 0; i < events.length; i++) {
    document.querySelector(selector).addEventListener(events[i], callback);
  }
}

function enter_event(selector, callback) {
  document
    .querySelector(selector)
    .addEventListener("keydown", event => { if (event.keyCode == 13) { callback(event); event.preventDefault(); } else { } });
}

function optional_number(text) {
  return text === "" ? null : Number(text);
}

function createElement(element, settings = undefined) {
  let node = document.createElement(element);
  if (settings != undefined) {
    if (settings.text != undefined) {
      node.innerText = settings.text;
    }
    if (settings.html != undefined) {
      node.innerHTML = settings.html;
    }
    if (settings.id != undefined) {
      node.id = settings.id;
    }
    if (settings.title != undefined) {
      node.title = settings.title;
    }
    if (settings.for != undefined) {
      node.setAttribute("for", settings.for)
    }
  }
  return node;
}

function showError(id, error, showContext = true) {
  let node = document.getElementById(id);
  if (error.warning) {
    node.classList.add("warning");
  } else {
    node.classList.add("error");
  }
  node.classList.remove("hidden");
  node.innerHTML = formatError(error, showContext);
}

function clearError(id) {
  let node = document.getElementById(id);
  node.classList.remove("error");
  node.classList.remove("warning");
  node.classList.add("hidden");
  node.innerHTML = "";
}