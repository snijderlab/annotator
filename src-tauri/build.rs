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
    <button id="load-clipboard" type="button">Load Clipboard</button>
    </div>
  <pre id="loaded-path"></pre>
  <div class="input-settings">
    <h2>Annotate</h2>
    <label for="spectrum-index">Spectrum index</label>
    <input type="number" id="spectrum-index" value="0" />
    <label for="spectrum-ppm">PPM </label>
    <input type="number" id="spectrum-ppm" value="20" />
    <label for="mass-system">Mass type </label>
    <select id="mass-system" disabled>
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
    <!--
    <div class='legend'>
      <span class='title'>Ion legend</span>
      <div class='ion-series'>
          <div class='top'>
              <span class='ion w' tabindex='0'>w</span>
              <span class='ion x' tabindex='0'>x</span>
              <span class='ion y' tabindex='0'>y</span>
              <span class='ion z' tabindex='0'>z</span>
              <span class='ion c-term' tabindex='0'>C-term</span>
          </div><div class='bottom'>
              <span class='ion n-term' tabindex='0'>N-term</span>
              <span class='ion a' tabindex='0'>a</span>
              <span class='ion b' tabindex='0'>b</span>
              <span class='ion c' tabindex='0'>c</span>
              <span class='ion d' tabindex='0'>d</span>
              <span class='ion v' tabindex='0'>v</span>
          </div>
      </div>
      <span class='other'>Other</span>
      <input id='spectrum_unassigned' type='checkbox' checked class='unassigned'/>
      <label for='spectrum_unassigned' class='unassigned' tabindex='0'>Unassigned</label>
      <label class='label'>
      Ion
      <sup>Charge</sup>
      <sub style='margin-left:-6ch;margin-right:.5rem;'>Position</sub>
      Show for top:
      <input id='spectrum_label' type='range' min='0' max='100' value='100'/>
      <input id='spectrum_label_value' type='number' min='0' max='100' value='100'/>
      %
      </label>
      <input id='spectrum_mass_label' type='checkbox' class='mass-label'/>
      <label for='spectrum_mass_label' class='mass-label' tabindex='0'>Show top masses</label>
    </div>
    -->
    <div class='wrapper unassigned'  id="spectrum-results-wrapper">
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
