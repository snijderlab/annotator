use std::io::Write;
use std::{fs::File, io::BufWriter};

fn main() {
    let file = File::create("../src/index.html").unwrap();
    println!("{:?}", file.metadata());
    let mut writer = BufWriter::new(file);
    write!(
        writer,
r#"<!DOCTYPE html>
<html lang="en">

<head>
  <meta charset="UTF-8" />
  <link rel="stylesheet" href="stitch-assets/styles.css" />
  <link rel="stylesheet" href="style.css" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Annotator</title>
  <script type="module" src="/main.js" defer></script>
  <script src="stitch-assets/script.js"></script>
</head>

<body>
  <button class="print" onclick="window.print()">Export</button>
  <div class="input-settings">
    <h2>Load spectra</h2>
    <label for="load-mgf-path">File</label>
    <button type="button" id="load-mgf-path">Select file</button>
    <label for="load-clipboard">Clipboard</label>
    <button id="load-clipboard" type="button">Load Clipboard</button>
    </div>
  <pre id="loaded-path"></pre>
  <div class="input-settings">
    <h2>Spectrum details</h2>
    <label for="spectrum-index">Spectrum index</label>
    <div class="combined-input">
      <input type="number" id="details-spectrum-index" value="0" min="0" />
      <span>/</span>
      <span id="number-of-scans">0</span>
    </div>
  </div>
  <pre id="spectrum-details"></pre>
  <pre id="loaded-path"></pre>
  <div class="input-settings">
    <h2>Annotate</h2>
    <label for="spectrum-ppm">PPM </label>
    <input type="number" id="spectrum-ppm" value="20" />
    <label for="spectrum-charge">Max charge </label>
    <input type="number" id="spectrum-charge" value="" placeholder="Empty takes peptide charge from raw data" />
    <label for="noise-threshold">Intensity noise filter</label>
    <input type="number" id="noise-threshold" value="" min="0" max="1" placeholder="Empty does not filter" />
    <label for="spectrum-model">Model </label>
    <select id="spectrum-model">
    <option value="all">All</option>
    <option value="ethcd">Ethcd</option>
    <option value="etcid">Etcid</option>
    <option value="cidhcd">CidHcd</option>
    <option value="etd">Etd</option>
    <option value="custom">Custom</option>
    </select>
    <fieldset class="custom-model">
    <legend>Custom model</legend>
    <p>Ion</p>
    <p>Location</p>
    <p>Loss</p>"#).unwrap();
    for ion in ["a", "b", "c", "d", "v", "w", "x", "y", "z"] {
        write!(
            writer,
            r#"<label>{0}</label>
          <div id="model-{0}-location" class="location">
        <select onchange="this.className=this.options[Number(this.value)].dataset.cls;">
          <option value="0" data-cls="arg-0" data-value="All">All</option>
          <option value="1" data-cls="arg-0" data-value="None" selected>None</option>
          <option value="2" data-cls="arg-1" data-value="SkipN">SkipN</option>
          <option value="3" data-cls="arg-1" data-value="SkipC">SkipC</option>
          <option value="4" data-cls="arg-1" data-value="TakeC">TakeC</option>
          <option value="5" data-cls="arg-2" data-value="SkipNC">SkipNC</option>
          <option value="6" data-cls="arg-2" data-value="TakeN">TakeN</option>
          </select>
          <input type="number" value="1" min="1">
          <input type="number" value="1" min="1">
          </div>
          <input type="text" id="model-{0}-loss" value=""/>"#,
            ion
        )
        .unwrap();
    }
    write!(
        writer,
        r#"<label>precursor</label>
      <input type="text" id="model-precursor-loss" value="" class="col-2"/>
      </fieldset>
      <label class="wide" for="peptide">Peptide sequence </label>
      <div class="peptide-input wide" id="peptide" contentEditable="true"></div>
      <button id="annotate-button" type="button" class="col-2">Annotate</button>
    </div>
  <div id="spectrum-error" class="hidden"></div>
  <div id='spectrum-wrapper' class="spectrum hidden" onload='SpectrumSetUp()'>
    <div class='settings render-setup'><p>Render setup</p>
      <label for='spectrum-width'>Width</label>
      <input id='spectrum-width' class='width' type='text' value='100%'/>
      <label for='spectrum-height'>Height</label>
      <input id='spectrum-height' class='height' type='text' value='250px'/>
      <label for='spectrum-fs-peptide'>Peptide font size</label>
      <input id='spectrum-fs-peptide' class='fs-peptide' type='text' value='1.25rem'/>
      <label for='spectrum-peptide-stroke'>Peptide stroke width</label>
      <input id='spectrum-peptide-stroke' class='stroke-peptide' type='text' value='2px'/>
      <label for='spectrum-fs-spectrum'>Spectrum font size</label>
      <input id='spectrum-fs-spectrum' class='fs-spectrum' type='text' value='1rem'/>
      <label for='spectrum-spectrum-stroke'>Spectrum stroke width</label>
      <input id='spectrum-spectrum-stroke' class='stroke-spectrum' type='text' value='2px'/>
      <input id='spectrum-compact' class='compact' type='checkbox'/>
      <label for='spectrum-compact'>Compact peptide</label>
    </div>
    <div class='settings spectrum-graph-setup'><p>Spectrum graph</p>
      <input type='radio' name='y-axis' id='absolute' value='absolute' checked/><label for='absolute'>Absolute</label>
      <input type='radio' name='y-axis' id='relative' value='relative'/><label for='relative'>Relative</label>
      <input type='checkbox' name='intensity' id='intensity' value='intensity'/><label for='intensity'>Intensity</label>
      <div class='manual-zoom'>
        <label for='y-min'>Y Min</label>
        <input id='y-min' class='y-min' type='number' value='-20'/>
        <label for='y-max'>Y Max</label>
        <input id='y-max' class='y-max' type='number' value='20'/>
      </div>
    </div>
    <div class='wrapper show-unassigned'  id="spectrum-results-wrapper">
    </div>
  </div>
  <details>
  <summary>Logs</summary>
  <div id="spectrum-fragments"></div>
  <pre id="spectrum-log"></pre>
  </details>
</body>

</html>"#
    )
    .unwrap();
    tauri_build::build()
}
