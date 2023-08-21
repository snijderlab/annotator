const { invoke } = window.__TAURI__.tauri;

let sequenceInputA;
let sequenceInputB;
let sequenceType;
let alignmentScore;

async function align() {
  alignmentScore.innerHTML = await invoke("align_sequences", { template: sequenceInputA.value, reads: sequenceInputB.value, alignmentType: sequenceType.value });
}

async function load_cif() {
  try {
    var result = await invoke("load_cif", { path: document.querySelector("#load-path").value, minLength: Number(document.querySelector("#load-min-length").value), warn: true });
    sequenceInputB.value = result[0];
    document.querySelector("#error-log").innerText = result[1];
  } catch (error) {
    console.log(error);
    document.querySelector("#error-log").innerText = error;
  }
}

async function load_mgf() {
  try {
    let result = await invoke("load_mgf", { path: document.querySelector("#load-mgf-path").dataset.filepath });
    document.querySelector("#spectrum-log").innerText = result;
    document.querySelector("#spectrum-error").innerText = "";
  } catch (error) {
    console.log(error);
    document.querySelector("#spectrum-error").innerText = error;
  }
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
  try {
    var charge = document.querySelector("#spectrum-charge").value == "" ? null : Number(document.querySelector("#spectrum-charge").value);
    var noise_threshold = document.querySelector("#noise-threshold").value == "" ? null : Number(document.querySelector("#noise-threshold").value);
    var model = [
      [get_location("#model-a-location"), document.querySelector("#model-a-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-b-location"), document.querySelector("#model-b-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-c-location"), document.querySelector("#model-c-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-d-location"), document.querySelector("#model-d-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-v-location"), document.querySelector("#model-v-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-w-location"), document.querySelector("#model-w-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-x-location"), document.querySelector("#model-x-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-y-location"), document.querySelector("#model-y-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-z-location"), document.querySelector("#model-z-loss").value.split(' ').filter(a => a != "")],
      [get_location("#model-z-location"), document.querySelector("#model-precursor-loss").value.split(' ').filter(a => a != "")], // First element is discarded
    ];
    var result = await invoke("annotate_spectrum", { index: Number(document.querySelector("#spectrum-index").value), ppm: Number(document.querySelector("#spectrum-ppm").value), mass: document.querySelector("#mass-system").value, charge: charge, noise_threshold: noise_threshold, model: document.querySelector("#spectrum-model").value, peptide: document.querySelector("#peptide").value, cmodel: model });
    document.querySelector("#spectrum-results-wrapper").innerHTML = result[0];
    document.querySelector("#spectrum-fragments").innerHTML = result[1];
    document.querySelector("#spectrum-log").innerText = result[2];
    document.querySelector("#spectrum-error").innerText = "";
    SpectrumSetUp();
  } catch (error) {
    console.log(error);
    document.querySelector("#spectrum-error").innerText = error;
  }
}

window.addEventListener("DOMContentLoaded", () => {
  // sequenceInputA = document.querySelector("#greet-input-a");
  // sequenceInputB = document.querySelector("#greet-input-b");
  // sequenceType = document.querySelector("#input-alignment-type");
  // alignmentScore = document.querySelector("#reads-alignment");
  // document
  //   .querySelector("#greet-button")
  //   .addEventListener("click", () => align());
  // document
  //   .querySelector("#load-button")
  //   .addEventListener("click", () => load_cif());
  document
    .querySelector("#load-mgf-button")
    .addEventListener("click", () => load_mgf());
  document
    .querySelector("#annotate-button")
    .addEventListener("click", () => annotate_spectrum());
});
