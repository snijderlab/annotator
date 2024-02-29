const { invoke } = window.__TAURI__.tauri;

const { listen } = window.__TAURI__.event
// const controller = new AbortController();

listen('tauri://file-drop', event => {
  document.querySelector("html").classList.remove("file-drop-hover");
  for (let i = 0; i < event.payload.length; i++) {
    let file = event.payload[i];
    if (file.toLowerCase().endsWith(".mgf") || file.toLowerCase().endsWith(".mgf.gz")) {
      load_mgf(file);
    } else {
      load_identified_peptides(file);
    }
  }
})

listen('tauri://file-drop-hover', event => {
  document.querySelector("html").classList.add("file-drop-hover");
})

listen('tauri://file-drop-cancelled', event => {
  document.querySelector("html").classList.remove("file-drop-hover");
})

// function abort() {
//   console.log(controller);
//   controller.abort("User manual abort");
// }

/**
* @param e: Element
*/
async function select_mgf_file(e) {
  let properties = {
    //defaultPath: 'C:\\',
    directory: false,
    filters: [{
      extensions: ['mgf', 'mgf.gz'], name: "*"
    }]
  };
  window.__TAURI__.dialog.open(properties).then((result) => {
    e.dataset.filepath = result;
    load_mgf(e.dataset.filepath);
  })
};

/**
* @param e: Element
*/
async function select_identified_peptides_file(e) {
  let properties = {
    //defaultPath: 'C:\\',
    directory: false,
    filters: [{
      extensions: ["csv", "csv.gz", "psmtsv", "psmtsv.gz", "fasta", "fasta.gz"], name: "*"
    }]
  };
  window.__TAURI__.dialog.open(properties).then((result) => {
    e.dataset.filepath = result;
    load_identified_peptides(e.dataset.filepath);
  })
};

async function load_mgf(path) {
  document.querySelector("#load-mgf-path").classList.add("loading")
  invoke("load_mgf", { path: path }).then((result) => {
    document.querySelector("#spectrum-log").innerText = "Loaded " + result + " spectra";
    document.querySelector("#loaded-path").classList.remove("error");
    document.querySelector("#loaded-path").innerText = path.split('\\').pop().split('/').pop();
    document.querySelector("#number-of-scans").innerText = result;
    spectrum_details();
    document.querySelector("#load-mgf-path").classList.remove("loading")
  }).catch((error) => {
    console.log(error);
    document.querySelector("#loaded-path").classList.add("error");
    document.querySelector("#loaded-path").innerText = error;
    document.querySelector("#load-mgf-path").classList.remove("loading")
  })
}

async function load_identified_peptides(path) {
  document.querySelector("#load-identified-peptides").classList.add("loading")
  invoke("load_identified_peptides", { path: path }).then((result) => {
    document.querySelector("#identified-peptides-log").innerText = "Loaded " + result + " peptides";
    document.querySelector("#loaded-identified-peptides-path").classList.remove("error");
    document.querySelector("#loaded-identified-peptides-path").innerText = path.split('\\').pop().split('/').pop();
    document.querySelector("#number-of-identified-peptides").innerText = result;
    document.querySelector("#peptides").style.display = "block";
    displayed_identified_peptide = undefined;
    identified_peptide_details();
    document.querySelector("#load-identified-peptides").classList.remove("loading")
  }).catch((error) => {
    console.log(error);
    document.querySelector("#loaded-identified-peptides-path").classList.add("error");
    document.querySelector("#loaded-identified-peptides-path").innerText = error;
    document.querySelector("#load-identified-peptides").classList.remove("loading")
  })
}

async function load_clipboard() {
  document.querySelector("#load-clipboard").classList.add("loading");
  navigator.clipboard
    .readText()
    .then(async (clipText) => {
      invoke("load_clipboard", { data: clipText }).then((result) => {
        document.querySelector("#spectrum-log").innerText = "Loaded " + result + " spectra";
        document.querySelector("#loaded-path").classList.remove("error");
        document.querySelector("#loaded-path").innerText = "Clipboard";
        document.querySelector("#number-of-scans").innerText = result;
        spectrum_details();
      }).catch((error) => {
        document.querySelector("#load-clipboard").classList.remove("loading")
        console.log(error);
        document.querySelector("#loaded-path").classList.add("error");
        document.querySelector("#loaded-path").innerText = error;
        document.querySelector("#load-clipboard").classList.remove("loading")
      })
    });
}

async function find_scan_number() {
  invoke("find_scan_number", { scanNumber: Number(document.querySelector("#scan-number").value) }).then(
    (result) => {
      document.querySelector("#details-spectrum-index").value = result;
      spectrum_details();
    }
  ).catch((error) => {
    console.log(error);
    document.querySelector("#spectrum-error").classList.remove("hidden");
    document.querySelector("#spectrum-error").innerText = error;
  })
}

async function spectrum_details() {
  invoke("spectrum_details", { index: Number(document.querySelector("#details-spectrum-index").value) }).then(
    (result) => document.querySelector("#spectrum-details").innerText = result
  ).catch((error) => {
    console.log(error);
    document.querySelector("#spectrum-error").classList.remove("hidden");
    document.querySelector("#spectrum-error").innerText = error;
    document.querySelector("#spectrum-details").innerText = "ERROR";
  })
}

let displayed_identified_peptide = undefined;
async function identified_peptide_details() {
  let index = Number(document.querySelector("#details-identified-peptide-index").value);
  if (displayed_identified_peptide != index) {
    invoke("identified_peptide_details", { index: index }).then((result) => {
      document.querySelector("#identified-peptide-details").innerHTML = result;
      displayed_identified_peptide = index;
    }).catch((error) => {
      console.log(error);
      document.querySelector("#spectrum-error").classList.remove("hidden");
      document.querySelector("#spectrum-error").innerText = error;
      document.querySelector("#identified-peptide-details").innerText = "ERROR";
    })
  }
}

async function search_peptide() {
  if (document.querySelector("#search-peptide-input").value != "") {
    document.querySelector("#search-peptide").classList.add("loading");
    invoke("search_peptide", { text: document.querySelector("#search-peptide-input").value }).then((result) => {
      document.querySelector("#resulting-peptides").innerHTML = result;
      document.querySelector("#search-peptide").classList.remove("loading");
    }).catch((error) => {
      console.log(error);
      document.querySelector("#spectrum-error").classList.remove("hidden");
      document.querySelector("#spectrum-error").innerText = error;
      document.querySelector("#identified-peptide-details").innerText = "ERROR";
      document.querySelector("#search-peptide").classList.remove("loading");
    })
  } else {
    document.querySelector("#resulting-peptides").innerText = "";
  }
}

async function load_identified_peptide() {
  let index = Number(document.querySelector("#details-identified-peptide-index").value);
  invoke("load_identified_peptide", { index: index }).then((result) => {
    document.querySelector("#peptide").innerText = result.peptide;
    document.querySelector("#spectrum-charge").value = result.charge;
    document.querySelector("#spectrum-model").value = result.mode.toLowerCase();
    document.querySelector("#details-spectrum-index").value = result.scan_index;
    console.log(result)
  }).catch((error) => {
    console.log(error);
    document.querySelector("#spectrum-error").classList.remove("hidden");
    document.querySelector("#spectrum-error").innerText = error;
    document.querySelector("#identified-peptide-details").innerText = "ERROR";
  })
}


function get_location(id) {
  let loc = document.querySelector(id);
  let t = loc.children[0].options[Number(loc.children[0].value)].dataset.value;
  // [[{\"SkipN\":1},[\"Water\"]],[\"All\",[\"Water\"]],[{\"SkipNC\":[1,2]},[\"Water\"]],[{\"TakeN\":{\"skip\":2,\"take\":1}},[\"Water\"]]]
  if (["SkipN", "SkipC", "TakeC"].includes(t)) {
    let obj = {};
    obj[t] = Number(loc.children[1].value);
    return obj;
  } else if (t == "SkipNC") {
    let obj = {};
    obj[t] = [Number(loc.children[1].value), Number(loc.children[2].value)];
    return obj;
  } else if (t == "TakeN") {
    return { t: { "skip": Number(loc.children[1].value), "take": Number(loc.children[2].value) } };
  } else {
    return t;
  }
}

function get_noise_filter(id) {
  let loc = document.querySelector(id);
  let t = loc.children[0].options[Number(loc.children[0].value)].dataset.value;
  if (["Relative", "Absolute"].includes(t)) {
    let obj = {};
    obj[t] = Number(loc.children[1].value);
    return obj;
  } else if (t == "TopX") {
    let obj = {};
    obj[t] = [Number(loc.children[1].value), Number(loc.children[2].value)];
    return obj;
  } else {
    return t;
  }
}

//import { SpectrumSetUp } from "./stitch-assets/script.js";
async function annotate_spectrum() {
  document.querySelector("#annotate-button").classList.add("loading");
  document.querySelector("#peptide").innerText = document.querySelector("#peptide").innerText.trim();
  var charge = document.querySelector("#spectrum-charge").value == "" ? null : Number(document.querySelector("#spectrum-charge").value);
  var noise_threshold = get_noise_filter("#noise-filter");
  console.log(noise_threshold, charge);
  var model = [
    [get_location("#model-a-location"), document.querySelector("#model-a-loss").value],
    [get_location("#model-b-location"), document.querySelector("#model-b-loss").value],
    [get_location("#model-c-location"), document.querySelector("#model-c-loss").value],
    [get_location("#model-d-location"), document.querySelector("#model-d-loss").value],
    [get_location("#model-v-location"), document.querySelector("#model-v-loss").value],
    [get_location("#model-w-location"), document.querySelector("#model-w-loss").value],
    [get_location("#model-x-location"), document.querySelector("#model-x-loss").value],
    [get_location("#model-y-location"), document.querySelector("#model-y-loss").value],
    [get_location("#model-z-location"), document.querySelector("#model-z-loss").value],
    [get_location("#model-z-location"), document.querySelector("#model-precursor-loss").value], // First element is discarded
  ];
  invoke("annotate_spectrum", { index: Number(document.querySelector("#details-spectrum-index").value), ppm: Number(document.querySelector("#spectrum-ppm").value), charge: charge, filter: noise_threshold, model: document.querySelector("#spectrum-model").value, peptide: document.querySelector("#peptide").innerText, cmodel: model }).then((result) => {
    document.querySelector("#spectrum-results-wrapper").innerHTML = result[0];
    document.querySelector("#spectrum-fragments").innerHTML = result[1];
    document.querySelector("#spectrum-log").innerText = result[2];
    document.querySelector("#spectrum-error").innerText = "";
    document.querySelector("#spectrum-wrapper").classList.remove("hidden"); // Remove hidden class if this is the first run
    document.querySelector("#spectrum-error").classList.add("hidden");
    SetUpSpectrumInterface();
    document.querySelector("#annotate-button").classList.remove("loading");
  }).catch((error) => {
    console.log(error);
    document.querySelector("#spectrum-error").classList.remove("hidden");
    document.querySelector("#spectrum-error").innerHTML = "<p class='title'>" + error.short_description + "</p><p class='description'>" + error.long_description + "</p>";
    if (error.context.hasOwnProperty('Line')) {
      let Line = error.context.Line;
      document.querySelector("#peptide").innerHTML = Line.line.slice(0, Line.offset) + "<span class='error'>" + Line.line.slice(Line.offset, Line.offset + Line.length) + "</span>" + Line.line.slice(Line.offset + Line.length, Line.line.length);
    }
    document.querySelector("#annotate-button").classList.remove("loading");
  })
}

/** @param e {MouseEvent}  */
function resizeDown(e) {
  document.querySelector(".resize-wrapper").classList.add("active");
  document.querySelector(".resize-wrapper").dataset.start_x = e.clientX;
  document.querySelector(".resize-wrapper").dataset.left_width = document.querySelector(".resize-wrapper").firstElementChild.getBoundingClientRect().width - 16;
  document.addEventListener("mousemove", resizeMove);
  document.addEventListener("mouseup", resizeUp);
  document.querySelector(".resize-wrapper").style.userSelect = 'none';
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

window.addEventListener("DOMContentLoaded", () => {
  document.querySelector(".resize").addEventListener("mousedown", resizeDown);
  document
    .querySelector("#load-mgf-path")
    .addEventListener("click", (event) => select_mgf_file(event.target));
  document
    .querySelector("#load-identified-peptides")
    .addEventListener("click", (event) => select_identified_peptides_file(event.target));
  document
    .querySelector("#load-identified-peptide")
    .addEventListener("click", (event) => load_identified_peptide());
  document
    .querySelector("#load-clipboard")
    .addEventListener("click", () => load_clipboard());
  document
    .querySelector("#search-peptide", search_peptide)
    .addEventListener("click", () => search_peptide());
  enter_event("#search-peptide-input", search_peptide)
  enter_event("#scan-number", find_scan_number)
  add_event("#details-spectrum-index", ["change", "focus"], spectrum_details)
  add_event("#details-identified-peptide-index", ["change", "focus"], identified_peptide_details)
  document
    .querySelector("#annotate-button")
    .addEventListener("click", () => annotate_spectrum());
  document
    .querySelector("#peptide")
    .addEventListener("focus", (event) => {
      event.target.innerHTML = event.target.innerText;
    });
  enter_event("#peptide", annotate_spectrum)

  // Refresh interface for hot reload
  invoke("refresh").then((result) => {
    document.querySelector("#number-of-scans").innerText = result[0];
    if (result[0] > 0) {
      spectrum_details();
    }
    document.querySelector("#number-of-identified-peptides").innerText = result[1];
    if (result[1] > 0) {
      document.querySelector("#peptides").style.display = "block";
      identified_peptide_details();
    }
  })
});

function add_event(selector, events, callback) {
  for (let i = 0; i < events.length; i++) {
    document.querySelector(selector).addEventListener(events[i], callback);
  }
}

function enter_event(selector, callback) {
  document
    .querySelector(selector)
    .addEventListener("keypress", event => { if (event.keyCode == 13) { callback(event); event.preventDefault(); } else { } });
}