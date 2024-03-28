use std::io::Write;
use std::{fs::File, io::BufWriter};

fn create_loss_modal(id: &str) -> String {
    format!(
        r#"<div><button onclick='document.getElementById("model-{id}-loss-selection-dialog").showModal();'>Select</button><output id='model-{id}-loss-selection-output'>0 selected</output>
  <dialog id="model-{id}-loss-selection-dialog" onclose='var num = 0; document.getElementsByName("model-{id}-loss-selection").forEach(e=>num += e.checked);document.getElementById("model-{id}-loss-selection-output").innerText = num + " selected";'>
    <p>Losses</p>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H2O"/>Water</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H4O2"/>Double water</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H6O3"/>Triple water</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H"/>Hydrogen</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H2"/>Double hydrogen</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H3"/>Triple hydrogen</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-NH3"/>Ammonia</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-CO"/>Carbon monoxide</label>
    <p>Gains</p>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H2O"/>Water</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H4O2"/>Double water</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H6O3"/>Triple water</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H"/>Hydrogen</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H2"/>Double hydrogen</label>
    <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H3"/>Triple hydrogen</label>
    <button autofocus onclick='this.parentElement.close()'>Close</button>
  </dialog></div>"#
    )
}

fn main() {
    let version = option_env!("CARGO_PKG_VERSION").unwrap_or("undefined");
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
      <button class="secondary print" onclick="window.print()">Export</button>
      <!-- <button class="print" id="abort">Abort</button> -->
      <div class="input-flex">
        <button type="button" id="load-mgf-path">Load raw data file</button>
        <button type="button" id="load-clipboard">Load Clipboard</button>
        <button type="button" id="load-identified-peptides">Load identified peptides file</button>
      </div>
      <pre id="loaded-path"></pre>
      <pre id="loaded-identified-peptides-path"></pre>
      <div class="input-settings">
        <h2>Spectrum details</h2>
        <label for="spectrum-index">Spectrum index</label>
        <div class="combined-input">
          <input type="number" id="details-spectrum-index" value="0" min="0" />
          <span>/</span>
          <span id="number-of-scans">0</span>
        </div>
        <label for="scan-number">Scan number</label>
        <input type="number" id="scan-number" value="0" min="0" />
      </div>
      <pre id="spectrum-details"></pre>
      
      <div id="peptides" style="display:none">
        <h2>Peptide details</h2>
        <div class="resize-wrapper">
          <div>
            <div class="input-flex">
              <label for="search-peptide-input">Search peptide</label>
              <input id="search-peptide-input" type="text"></input>
              <label for="search-peptide-minimal-match">Minimal match score</label>
              <input id="search-peptide-minimal-match" type="number" min="0" max="1" value="0"></input>
              <label for="search-peptide-minimal-peptide">Minimal peptide score</label>
              <input id="search-peptide-minimal-peptide" type="number" min="0" max="1" value="0"></input>
              <button id="search-peptide">Search</button>
            </div>
            <div id="resulting-peptides">Go and search!</div>
          </div>
          <div class="resize"></div>
          <div>
            <div class="input-flex">
              <label for="spectrum-index">Identified peptide index</label>
              <div class="combined-input">
                <input type="number" id="details-identified-peptide-index" value="0" min="0" />
                <span>/</span>
                <span id="number-of-identified-peptides">0</span>
              </div>
              <button id="load-identified-peptide" type="button">Load identified peptide</button>
            </div>
            <div id="identified-peptide-details"></div>
          </div>
        </div>
      </div>
      <pre id="loaded-path"></pre>
      <div class="input-settings">
        <h2>Annotate</h2>
        <label for="spectrum-ppm">PPM </label>
        <input type="number" id="spectrum-ppm" value="20" />
        <label for="spectrum-charge">Max charge </label>
        <input type="number" id="spectrum-charge" value="" placeholder="Empty takes peptide charge from raw data" />
        
        
        <label for="noise-filter">Noise filter</label>
        <div id="noise-filter" class="noise-filter select-input">
          <select onchange="this.className=this.options[Number(this.value)].dataset.cls;">
            <option value="0" data-cls="arg-0" data-value="None" selected>None</option>
            <option value="1" data-cls="arg-1" data-value="Relative">Relative intensity</option>
            <option value="2" data-cls="arg-1" data-value="Absolute">Absolute intensity</option>
            <option value="3" data-cls="arg-2" data-value="TopX">TopX</option>
          </select>
          <input type="number" value="1" min="0">
          <input type="number" value="1" min="0">
        </div>
    
        <label for="spectrum-model">Model </label>
        <select id="spectrum-model">
          <option value="all" title="All possible ions with single water loss from all and additionally double water loss from glycans">All</option>
          <option value="ethcd" title="b+c+w+y+z+glycan with single water loss from all and additionally double water loss from glycans">EThcD/ETcaD</option>
          <option value="cidhcd" title="a+b+d+y+precursor with single water loss">CID/HCD</option>
          <option value="etd" title="c+y+z with single water loss and precursor with single water and ammonia loss">ETD</option>
          <option value="none" title="Only the base precursor peak, no losses">None</option>
          <option value="custom">Custom</option>
        </select>
        <fieldset class="custom-model">
          <legend>Custom model</legend>
          <p>Ion</p>
          <p>Location</p>
          <p>Loss</p>
          <p>Custom loss</p>"#).unwrap();
    for ion in ["a", "b", "c", "d", "v", "w", "x", "y", "z"] {
        write!(
                writer,
                r#"<label>{0}</label>
              <div id="model-{0}-location" class="location select-input">
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
              {1}
              <input type="text" id="model-{0}-loss" value="" title="Supply all losses as +/- followed by the chemical formula, supply multiple by separating them by commas. Example: '+H2O,-H2O'."/>"#,
                ion, create_loss_modal(ion)
            )
            .unwrap();
    }
    write!(
            writer,
            r#"<label>precursor</label>
            <div class="empty"></div>
            {}
            <input type="text" id="model-precursor-loss" value="" title="Supply all losses as +/- followed by the chemical formula, supply multiple by separating them by commas. Example: '+H2O,-H2O'."/>
            <label>immonium</label>
            <label><input id='model-immonium-enabled' type='checkbox' switch/>Enable</label>
            <label>glycan</label>
            <label><input id='model-glycan-enabled' type='checkbox' switch/>Enable</label>
            {}
            <input type="text" id="model-glycan-loss" value="" title="Supply all losses as +/- followed by the chemical formula, supply multiple by separating them by commas. Example: '+H2O,-H2O'."/>
          </fieldset>
          <label class="wide" for="peptide">Peptide sequence </label>
          <div class="peptide-input wide" id="peptide" contentEditable="plaintext-only"></div>
          <button id="annotate-button" type="button" class="col-2">Annotate</button>
        </div>
      <output id="spectrum-error" class="hidden error"></output>
      <div id='spectrum-wrapper' class="spectrum show-unassigned hidden legend-ion" onload='SpectrumSetUp()'>
        <div class='all-settings'>
          <fieldset class='settings graphics-settings'>
            <legend>Graphics settings</legend>
            <label class='row align' for='spectrum-width'><span>Width</span><input id='spectrum-width' class='width' type='text' value='100%'/></label>
            <label class='row align' for='spectrum-height'><span>Height</span><input id='spectrum-height' class='height' type='text' value='250px'/></label>
            <label class='row align' for='spectrum-fs-peptide'><span>Peptide font size</span><input id='spectrum-fs-peptide' class='fs-peptide' type='text' value='1.25rem'/></label>
            <label class='row align' for='spectrum-peptide-stroke'><span>Peptide stroke width</span><input id='spectrum-peptide-stroke' class='stroke-peptide' type='text' value='2px'/></label>
            <label class='row align' for='spectrum-fs-spectrum'><span>Spectrum font size</span><input id='spectrum-fs-spectrum' class='fs-spectrum' type='text' value='1rem'/></label>
            <label class='row align' for='spectrum-spectrum-stroke'><span>Spectrum stroke width</span><input id='spectrum-spectrum-stroke' class='stroke-spectrum' type='text' value='2px'/></label>
          </fieldset>
    
          <fieldset class='settings spectrum-graph-settings'>
            <legend>Spectrum graph settings</legend>
    
            <div class='row'>
              <span class='title'>Y axis mode</span>
              <div class='select-box' id='spectrum-graph-type'>
                <label for='spectrum-graph-type-absolute' tabindex='0'><input type='radio' name='spectrum-graph-type' value='absolute' id='spectrum-graph-type-absolute' checked>Absolute</label>
                <label for='spectrum-graph-type-relative' tabindex='0'><input type='radio' name='spectrum-graph-type' value='relative' id='spectrum-graph-type-relative'>Relative</label>
              </div>
            </div>
    
            <label for='intensity'><input type='checkbox' name='intensity' id='intensity' value='intensity' switch/>Intensity</label>
            
            <div class='manual-zoom'>
              <span>Y</span>
              <label for='y-min'>min</label>
              <input id='spectrum-graph-y-min' class='y-min' type='number' value='-20'/>
              <label for='y-max'>max</label>
              <input id='spectrum-graph-y-max' class='y-max' type='number' value='20'/>
            </div>
          </fieldset>
    
          <fieldset class='settings peptide-settings'>
            <legend>Peptide settings</legend>
            <label for='spectrum-compact' title='Display the peptide ion support in a more compact way'><input id='spectrum-compact' class='compact' type='checkbox' switch/>Compact peptide</label>
    
            <div class='row'>
              <span class='title'>Highlight</span>
              <div class='select-box' id='highlight'>
                <label title='Highlight a region in a peptide, while not overriding the ion colours' tabindex='0'><input type='radio' name='highlight' value='default' id='highlight-default' checked>Default</label>
                <label title='Annotate regions in a peptide in red' class='colour' tabindex='0'><input type='radio' name='highlight' value='red' id='highlight-red'></label>
                <label title='Annotate regions in a peptide in green' class='colour' tabindex='0'><input type='radio' name='highlight' value='green' id='highlight-green'></label>
                <label title='Annotate regions in a peptide in blue' class='colour' tabindex='0'><input type='radio' name='highlight' value='blue' id='highlight-blue'></label>
                <label title='Annotate regions in a peptide in yellow' class='colour' tabindex='0'><input type='radio' name='highlight' value='yellow' id='highlight-yellow'></label>
                <label title='Annotate regions in a peptide in purple' class='colour' tabindex='0'><input type='radio' name='highlight' value='purple' id='highlight-purple'></label>
                <label title='Remove the annotation for regions in a peptide' tabindex='0'><input type='radio' name='highlight' value='remove' id='highlight-remove'>X</label>
              </div>
              <button id='clear-colour' class='clear-colour' title='Remove all annotations on all peptides' tabindex='0'>Clear</button>
            </div>
    
          </fieldset>
    
          <fieldset class='settings spectrum-settings'>
            <legend>Spectrum settings</legend>
            
            <label for='theoretical' title='Show the theoretical peptide spectrum on the x axis'><input id='theoretical' class='theoretical' type='checkbox' switch/>Theoretical spectrum</label>
            
            <label for='unassigned' title='Show the unassigned peaks in the spectrum'><input id='unassigned' class='unassigned' type='checkbox' switch checked/>Unassigned</label>
    
            <div class='row'>
              <span class='title'>Ion colouration mode</span>
              <div class='select-box' id='peak-colour'>
                <label for='peak-colour-ion' tabindex='0'><input type='radio' name='peak-colour' value='ion' id='peak-colour-ion' checked>Ion</label>
                <label for='peak-colour-peptide' tabindex='0'><input type='radio' name='peak-colour' value='peptide' id='peak-colour-peptide'>Peptide</label>
                <label for='peak-colour-none' tabindex='0'><input type='radio' name='peak-colour' value='none' id='peak-colour-none'>None</label>
              </div>
            </div>
    
            <label class='has-slider label row align'>
              <span>Show labels for top:</span>
              <input id='spectrum-label' type='range' min='0' max='100' value='100'/>
              <input id='spectrum-label-value' type='number' min='0' max='100' value='100'/>
              %
            </label>
            
            <label class='has-slider masses row align'>
              <span>Show masses for top:</span>
              <input id='spectrum-masses' type='range' min='0' max='100' value='0'/>
              <input id='spectrum-masses-value' type='number' min='0' max='100' value='0'/>
              %
            </label>
    
            <div class='row'>
              <span class='title'>Manually force show</span>
              <div class='select-box' id='force-show'>
                <label tabindex='0'><input type='radio' name='force-show' value='none' id='force-show-none' checked>None</label>
                <label tabindex='0'><input type='radio' name='force-show' value='label' id='force-show-label'>Label</label>
                <label tabindex='0'><input type='radio' name='force-show' value='mass' id='force-show-mass'>Mass</label>
              </div>
              <button tabindex='0' id='force-show-clear'>Clear</button>
            </div>
            
            <div>
              <span>Mz</span>
              <label for='spectrum-mz-min'>min</label>
              <input id='spectrum-mz-min' class='mz-min' type='number' value='0' minimum='0'/>
              <label for='spectrum-mz-max'>max</label>
              <input id='spectrum-mz-max' class='mz-max' type='number' value='42' minimum='0'/>
            </div>
    
            <div>
              <span>Intensity</span>
              <label for='spectrum-intensity-max'>max</label>
              <input id='spectrum-intensity-max' class='intensity-max' type='number' value='42' minimum='0'/>
            </div>
    
            <button id='reset-zoom' class='reset-zoom' title='Reset the zoom to the default' tabindex='0'>Reset zoom</button>
          </fieldset>
        </div>
        <div class='legend'>
          <span class='title'>Ion legend</span>
          <div class='ion-series'>
              <div class='top'>
                  <span class='ion w' tabindex='0'>w</span>
                  <span class='ion x' tabindex='0'>x</span>
                  <span class='ion y' tabindex='0'>y</span>
                  <span class='ion z' tabindex='0'>z</span>
                  <span class='ion c-term' tabindex='0'>C-term</span>
              </div>
              <div class='bottom'>
                  <span class='ion n-term' tabindex='0'>N-term</span>
                  <span class='ion a' tabindex='0'>a</span>
                  <span class='ion b' tabindex='0'>b</span>
                  <span class='ion c' tabindex='0'>c</span>
                  <span class='ion d' tabindex='0'>d</span>
                  <span class='ion v' tabindex='0'>v</span>
              </div>
              <div class='side'>
                  <span class='ion precursor' tabindex='0'>Precursor</span>
                  <span class='ion multi mp' tabindex='0'>Multi</span>
                  <span class='ion glycan' tabindex='0'>Glycan</span>
                  <span class='other'>Other</span>
              </div>
          </div>
        </div>
        <output class='wrapper show-unassigned' id="spectrum-results-wrapper"></output>
    
        <input type="checkbox" id="collapsible-fragment-table">
        <fieldset class="collapsible" data-linked-item="collapsible-fragment-table">
          <legend>Fragment Table</legend>
          <output class="collapsible-content" id="spectrum-fragment-table"></output>
        </fieldset>
    
      </div>
    
      <input type="checkbox" id="collapsible-tools">
      <fieldset class="collapsible" data-linked-item="collapsible-tools">
        <legend>Tools</legend>
        <div class="flex-input collapsible-content">
          <label for="search-modification">Modification</label>
          <input id="search-modification" type="text"></input>
          <label for="search-modification-tolerance" title="If searching for a numeric modification, the tolerance on the mass">Tolerance (in Da)</label>
          <input id="search-modification-tolerance" type="number" value="0.1" min="0"></input>
          <button id="search-modification-button">Search</button>
        </div>
        <output class="collapsible-content" id="search-modification-result"></output>
        <output class="error collapsible-content hidden" id="search-modification-error"></output>
      </fieldset>
      <input type="checkbox" id="collapsible-logs">
      <fieldset class="collapsible" data-linked-item="collapsible-logs">
      <legend>Logs</legend>
      <output class="collapsible-content" id="spectrum-log"></output>
      <output class="collapsible-content" id="identified-peptides-log"></output>
      </fieldset>
      <div class="grow"></div>
      <footer>
        <p>Annotator by Snijderlab</p>
        <a href='https://github.com/snijderlab/annotator'>Open source at github (bugs can be reported here)</a>
        <p>Licensed under MIT OR Apache-2.0</p>
        <p>Version {}</p>
      </footer>
    </body>
    
    </html>"#,
    create_loss_modal("precursor"),
    create_loss_modal("glycan"),
    version
        )
        .unwrap();
    tauri_build::build()
}
