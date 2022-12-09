const { invoke } = window.__TAURI__.tauri;

let sequenceInputA;
let sequenceInputB;
let sequenceType;
let alignmentScore;

async function align() {
  // Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
  alignmentScore.innerHTML = await invoke("align_sequences", { template: sequenceInputA.value, reads: sequenceInputB.value, alignmentType: sequenceType.value });
}

async function load_cif() {
  // Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
  try {
    var result = await invoke("load_cif", { path: document.querySelector("#load-path").value, minLength: Number(document.querySelector("#load-min-length").value), warn: true });
    console.log(result);
    sequenceInputB.value = result[0];
    document.querySelector("#error-log").innerText = result[1];
  } catch (error) {
    console.log(error);
    document.querySelector("#error-log").innerText = error;
  }
}

window.addEventListener("DOMContentLoaded", () => {
  sequenceInputA = document.querySelector("#greet-input-a");
  sequenceInputB = document.querySelector("#greet-input-b");
  sequenceType = document.querySelector("#input-alignment-type");
  alignmentScore = document.querySelector("#reads-alignment");
  document
    .querySelector("#greet-button")
    .addEventListener("click", () => align());
  document
    .querySelector("#load-button")
    .addEventListener("click", () => load_cif());
});
