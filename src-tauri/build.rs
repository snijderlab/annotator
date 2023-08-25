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
  <script src="/basic.js" defer></script>
</head>

<body>
  <button class="print" onclick="window.print()">Export 🖶</button>
  <div class="input-settings">
    <h2>Load spectra</h2>
    <label for="load-mgf-path">Path</label>
    <button type="button" onclick="select_file(this)" id="load-mgf-path">Select file</button>
    <button id="load-mgf-button" type="button">Load</button>
  </div>
  <div class="input-settings">
    <h2>Annotate</h2>
    <label for="spectrum-index">Spectrum index</label>
    <input type="number" id="spectrum-index" value="0" />
    <label for="spectrum-ppm">PPM </label>
    <input type="number" id="spectrum-ppm" value="20" />
    <label for="mass-system">Mass type </label>
    <select id="mass-system">
    <option value="monoisotopic">MonoIsotopic</option>
    <option value="averageweight">AverageWeight</option>
    </select>
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
      <textarea class="wide" id="peptide"></textarea>
      <button id="annotate-button" type="button">Annotate</button>
    </div>
  <pre id="spectrum-error"></pre>
  <div id="spectrum-results-wrapper"></div>
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
