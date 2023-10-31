const { invoke } = window.__TAURI__.tauri;

/**
* @param e: Element
*/
async function select_mgf_file(e) {
  e.classList.add("loading")
  let properties = {
    //defaultPath: 'C:\\',
    directory: false,
    filters: [{
      extensions: ['mgf', 'mgf.gz'], name: "*"
    }]
  };
  window.__TAURI__.dialog.open(properties).then((result) => {
    e.dataset.filepath = result;
    load_mgf(e);
  })
};

/**
* @param e: Element
*/
async function select_identified_peptides_file(e) {
  e.classList.add("loading")
  let properties = {
    //defaultPath: 'C:\\',
    directory: false,
    filters: [{
      extensions: ["csv", "csv.gz", "psmtsv", "psmtsv.gz", "fasta", "fasta.gz"], name: "*"
    }]
  };
  window.__TAURI__.dialog.open(properties).then((result) => {
    e.dataset.filepath = result;
    load_identified_peptides(e);
  })
};

/**
* @param e: Element
*/
async function load_mgf(e) {
  invoke("load_mgf", { path: e.dataset.filepath }).then((result) => {
    document.querySelector("#spectrum-log").innerText = "Loaded " + result + " spectra";
    document.querySelector("#loaded-path").classList.remove("error");
    document.querySelector("#loaded-path").innerText = e.dataset.filepath.split('\\').pop().split('/').pop();
    document.querySelector("#number-of-scans").innerText = result;
    spectrum_details();
    e.classList.remove("loading")
  }).catch((error) => {
    console.log(error);
    document.querySelector("#loaded-path").classList.add("error");
    document.querySelector("#loaded-path").innerText = error;
    e.classList.remove("loading")
  })
}

/**
* @param e: Element
*/
async function load_identified_peptides(e) {
  invoke("load_identified_peptides", { path: e.dataset.filepath }).then((result) => {
    document.querySelector("#identified-peptides-log").innerText = "Loaded " + result + " peptides";
    document.querySelector("#loaded-identified-peptides-path").classList.remove("error");
    document.querySelector("#loaded-identified-peptides-path").innerText = e.dataset.filepath.split('\\').pop().split('/').pop();
    document.querySelector("#number-of-identified-peptides").innerText = result;
    displayed_identified_peptide = undefined;
    identified_peptide_details();
    e.classList.remove("loading")
  }).catch((error) => {
    console.log(error);
    document.querySelector("#loaded-identified-peptides-path").classList.add("error");
    document.querySelector("#loaded-identified-peptides-path").innerText = error;
    e.classList.remove("loading")
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

//import { SpectrumSetUp } from "./stitch-assets/script.js";
async function annotate_spectrum() {
  document.querySelector("#annotate-button").classList.add("loading");
  document.querySelector("#peptide").innerText = document.querySelector("#peptide").innerText.trim();
  var charge = document.querySelector("#spectrum-charge").value == "" ? null : Number(document.querySelector("#spectrum-charge").value);
  var noise_threshold = document.querySelector("#noise-threshold").value == "" ? null : Number(document.querySelector("#noise-threshold").value);
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
  invoke("annotate_spectrum", { index: Number(document.querySelector("#details-spectrum-index").value), ppm: Number(document.querySelector("#spectrum-ppm").value), charge: charge, noise_threshold: noise_threshold, model: document.querySelector("#spectrum-model").value, peptide: document.querySelector("#peptide").innerText, cmodel: model }).then((result) => {
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

window.addEventListener("DOMContentLoaded", () => {
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
    .querySelector("#search-peptide")
    .addEventListener("click", () => search_peptide());
  document
    .querySelector("#search-peptide-input")
    .addEventListener("keypress", event => event.keyCode == 13 ? search_peptide() : {});
  let dsi = document
    .querySelector("#details-spectrum-index");
  dsi.addEventListener("change", () => spectrum_details())
  dsi.addEventListener("focus", () => spectrum_details());
  let dipi = document
    .querySelector("#details-identified-peptide-index");
  dipi.addEventListener("change", () => identified_peptide_details())
  dipi.addEventListener("focus", () => identified_peptide_details());
  document
    .querySelector("#annotate-button")
    .addEventListener("click", () => annotate_spectrum());
  document
    .querySelector("#peptide")
    .addEventListener("focus", (event) => {
      event.target.innerHTML = event.target.innerText;
    });

  // Refresh interface for hot reload
  invoke("refresh").then((result) => {
    document.querySelector("#number-of-scans").innerText = result[0];
    if (result[0] > 0) {
      spectrum_details();
    }
    document.querySelector("#number-of-identified-peptides").innerText = result[1];
    if (result[1] > 0) {
      identified_peptide_details();
    }
  })
});
