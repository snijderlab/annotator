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
      <!-- <button class="secondary capture" id="capture">Export bitmap</button> -->
      <button class="secondary cancel-drop" onclick='document.querySelector("html").classList.remove("file-drop-hover")'>Cancel drop</button>
      <!-- <button class="print" id="abort">Abort</button> -->
      <div class="input-flex">
        <button type="button" id="load-mgf-path">Load raw data file</button>
        <button type="button" id="load-clipboard">Load Clipboard</button>
        <button type="button" id="load-identified-peptides">Load identified peptides file</button>
      </div>
      <output id="loaded-path"></output>
      <output id="loaded-identified-peptides-path"></output>
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
      <output id="spectrum-details"></output>
      
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
      <output id="loaded-path"></output>
      <div class="input-settings">
        <h2>Annotate</h2>
        <label for="spectrum-tolerance">Tolerance</label>
        <div class="grouped-input">
          <input type="number" id="spectrum-tolerance" value="20" />
          <select id="spectrum-tolerance-unit">
            <option value="ppm">ppm</option>
            <option value="th">Thompson (mz)</option>
          </select>
        </div>
        
        <label for="spectrum-mass-mode">Match mode</label>
        <select id="spectrum-mass-mode">
          <option value="monoisotopic" title="Match peaks based on the monoisotopic mass, meaning using only the most abundant isotope for each element (eg H1, C12 etc).">Monoisotopic</option>
          <option value="average_weight" title="Match the peaks based on average weight, meaning to use the average weight of each element based on its natural distribution.">Average weight</option>
          <option value="most_abundant" title="Match the peaks based on the most abundant isotope, meaning that for each fragment the isotopic distribution is generated and the most abundant species is selected. This mode should be good for very high mass matching, but it will take more time to be calculated.">Most abundant</option>
        </select>

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
              <div class="separated-input">
                <div class="values" id="model-{0}-loss">
                  <div class="input context" placeholder="Add molecular formula or mass preceded by '+' or '-'" data-type="neutral_loss" contentEditable="plaintext-only"></div>
                  <button class="clear">Clear</button>
                </div>
                <output class="error"></output>
              </div>"#,
                ion, create_loss_modal(ion)
            )
            .unwrap();
    }
    write!(
            writer,
            r#"<label>precursor</label>
            <div class="empty"></div>
            {}
            <div class="separated-input">
              <div class="values" id="model-precursor-loss">
                <div class="input context" placeholder="Add molecular formula or mass preceded by '+' or '-'" data-type="neutral_loss" contentEditable="plaintext-only"></div>
                <button class="clear">Clear</button>
              </div>
              <output class="error"></output>
            </div>
            <div class="grid-row">
              <label>immonium</label>
              <label><input id='model-immonium-enabled' type='checkbox' switch/>Enable</label>
            </div>
            <div class="grid-row">
              <label title="side chain loss from precursor">m</label>
              <label><input id='model-m-enabled' type='checkbox' switch/>Enable</label>
            </div>
            <div class="grid-row">
              <label title="Allow modification specific neutral losses, as defined by the database">modification neutral losses</label>
              <label><input id='model-modification-neutral-enabled' type='checkbox' switch/>Enable</label>
            </div>
            <div class="grid-row">
              <label title="Allow modification specific diagnostic ions, as defined by the database">modification diagnostic ions</label>
              <label><input id='model-modification-diagnostic-enabled' type='checkbox' switch/>Enable</label>
            </div>
            <label>glycan</label>
            <label><input id='model-glycan-enabled' type='checkbox' switch/>Enable</label>
            {}
            <div class="separated-input">
              <div class="values" id="model-glycan-loss">
                <div class="input context" placeholder="Add molecular formula or mass preceded by '+' or '-'" data-type="neutral_loss" contentEditable="plaintext-only"></div>
                <button class="clear">Clear</button>
              </div>
              <output class="error"></output>
            </div>
          </fieldset>
          <label class="wide" for="peptide">Peptide sequence </label>
          <div class="peptide-input wide context" id="peptide" contentEditable="plaintext-only"></div>
          <button id="annotate-button" type="button" class="col-2">Annotate</button>
        </div>
      <output id="spectrum-error" class="hidden error"></output>
      <div id='spectrum-wrapper' class="spectrum show-unassigned legend-ion hidden show-charge show-series show-glycan-id show-peptide-id show-neutral-losses" onload='SpectrumSetUp()'>
        <div class='all-settings'>
          <fieldset class='settings graphics-settings'>
            <legend>Graphics settings</legend>
            <label class='row align' for='spectrum-width'><span>Width</span><input id='spectrum-width' class='width' type='text' value='100%'/></label>
            <label class='row align' for='spectrum-height'><span>Height</span><input id='spectrum-height' class='height' type='text' value='250px'/></label>
            <label class='row align' for='spectrum-fs-peptide'><span>Peptide font size</span><input id='spectrum-fs-peptide' class='fs-peptide' type='text' value='1.25rem'/></label>
            <label class='row align' for='spectrum-peptide-stroke'><span>Peptide stroke width</span><input id='spectrum-peptide-stroke' class='stroke-peptide' type='text' value='2px'/></label>
            <label class='row align' for='spectrum-fs-spectrum'><span>Spectrum font size</span><input id='spectrum-fs-spectrum' class='fs-spectrum' type='text' value='1rem'/></label>
            <label class='row align' for='spectrum-spectrum-stroke'><span>Spectrum stroke width</span><input id='spectrum-spectrum-stroke' class='stroke-spectrum' type='text' value='2px'/></label>
            <label class='row align' for='spectrum-spectrum-stroke-unassigned'><span>Spectrum unassigned stroke width</span><input id='spectrum-spectrum-stroke-unassigned' class='stroke-spectrum-unassigned' type='text' value='1px'/></label>
          </fieldset>
    
          <fieldset class='settings error-graph-settings'>
            <legend>Error graph settings</legend>
    
            <div class='row'>
              <span class='title'>Y axis mode</span>
              <div class='select-box' id='error-graph-y-type'>
              <label tabindex='0'><input type='radio' name='error-graph-y-type' value='relative' id='error-graph-y-type-relative' checked>Relative</label>
              <label tabindex='0'><input type='radio' name='error-graph-y-type' value='absolute' id='error-graph-y-type-absolute'>Absolute</label>
              </div>
            </div>
    
            <label><input type='checkbox' name='error-graph-intensity' id='error-graph-intensity' value='error-graph-intensity' switch/>Intensity</label>
            
            <div class='manual-zoom'>
              <span>Y</span>
              <label for='y-min'>min</label>
              <input id='error-graph-y-min' class='y-min' type='number' value='-20'/>
              <label for='y-max'>max</label>
              <input id='error-graph-y-max' class='y-max' type='number' value='20'/>
            </div>

            <div class='row'>
              <span class='title'>X axis mode</span>
              <div class='select-box' id='error-graph-x-type'>
                <label tabindex='0'><input type='radio' name='error-graph-x-type' value='assigned' id='error-graph-x-type-assigned' checked>Assigned</label>
                <label tabindex='0'><input type='radio' name='error-graph-x-type' value='unassigned' id='error-graph-x-type-unassigned'>Unassigned</label>
              </div>
            </div>

            <fieldset>
              <legend>Unassigned mode settings</legend>
              
              <div class='row'>
                <span class='title'>Reference ions</span>
                <div class='error-graph-ion-selection multi-checkbox'>
                  <label><input id='error-graph-ion-a' type='checkbox'/>a</label>
                  <label><input id='error-graph-ion-b' type='checkbox'/>b</label>
                  <label><input id='error-graph-ion-c' type='checkbox'/>c</label>
                  <label><input id='error-graph-ion-x' type='checkbox'/>x</label>
                  <label><input id='error-graph-ion-y' type='checkbox'/>y</label>
                  <label><input id='error-graph-ion-z' type='checkbox'/>z</label>
                </div>
              </div>

              <label><input type='checkbox' name='error-graph-sh ow-assigned' id='error-graph-show-assigned' value='error-graph-show-assigned' switch/>Show assigned peaks</label>

            </fieldset>
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
            
            <label class='has-slider mz row align'>
              <span>Show m/z for top:</span>
              <input id='spectrum-m-z' type='range' min='0' max='100' value='0'/>
              <input id='spectrum-m-z-value' type='number' min='0' max='100' value='0'/>
              %
            </label>
    
            <div class='row'>
              <span class='title'>Manually force show</span>
              <div class='select-box' id='force-show'>
                <label tabindex='0'><input type='radio' name='force-show' value='none' id='force-show-none' checked>None</label>
                <label tabindex='0'><input type='radio' name='force-show' value='label' id='force-show-label'>Label</label>
                <label tabindex='0'><input type='radio' name='force-show' value='m-z' id='force-show-m-z'>m/z</label>
              </div>
              <button tabindex='0' id='force-show-clear'>Clear</button>
            </div>

            <div>
              <span class='title'>Show in label</span>
                <div class='spectrum-label-parts-selection multi-checkbox'>
                  <label><input checked id='spectrum-label-charge' type='checkbox'/>Charge</label>
                  <label><input checked id='spectrum-label-series' type='checkbox'/>Series number</label>
                  <label><input checked id='spectrum-label-glycan-id' type='checkbox'/>Glycan position</label>
                  <label><input checked id='spectrum-label-peptide-id' type='checkbox'/>Peptide number</label>
                  <label><input checked id='spectrum-label-neutral-losses' type='checkbox'/>Neutral losses</label>
                </div>
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
                  <span class='ion oxonium' tabindex='0'>Oxonium</span>
                  <span class='ion Y' tabindex='0'>Y</span>
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
        <output class="error collapsible-content hidden" id="search-modification-error"></output>
        <output class="collapsible-content" id="search-modification-result"></output>
        <div class="flex-input collapsible-content">
          <label for="details-formula">Formula</label>
          <input id="details-formula" type="text"></input>
        </div>
        <output class="error collapsible-content hidden" id="details-formula-error"></output>
        <output class="collapsible-content" id="details-formula-result"></output>
      </fieldset>
      <input type="checkbox" id="collapsible-custom-mods">
      <fieldset class="collapsible" data-linked-item="collapsible-custom-mods">
        <legend>Custom mods</legend>
        <dialog id="custom-mod-dialog">
          <h1>Custom modification</h1>

          <h2>Basic features</h2>
          <div class="basic">
            <label>Formula<span class="input context" id="custom-mod-formula" type="text" placeholder="Molecular formula or mass" contentEditable="plaintext-only"></span></label>
            <label>Id<input id="custom-mod-id" type="number" disabled value="0"></input></label>
            <label>Name<input id="custom-mod-name" type="text" placeholder="Identifying name"></input></label>
            <p class="justify-end">Use as follows:</p>
            <span>Id<span class="example" id="custom-mod-example-id">CUSTOM:0</span></span>
            <span>Name<span class="example" id="custom-mod-example-name">C:NAME</span></span>
            <output class="error"></output>
          </div>

          <h2>Metadata</h2>
          <label for="custom-mod-description">Description</label>
          <textarea id="custom-mod-description" class="wide" placeholder="Describe the modification"></textarea>
          <label for="custom-mod-synonyms">Synonyms</label>
          <div class="separated-input">
            <div class="values" id="custom-mod-synonyms">
              <div class="input context" placeholder="Add synonyms" data-type="text" contentEditable="plaintext-only"></div>
              <button class="clear">Clear</button>
            </div>
            <output class="error"></output>
          </div>
          <label for="custom-mod-cross-ids">Cross IDs</label>
          <div class="separated-input">
            <div class="values" id="custom-mod-cross-ids">
              <div class="input context" placeholder="Specify like: 'System:ID'" data-type="cross_id" contentEditable="plaintext-only"></div>
              <button class="clear">Clear</button>
            </div>
            <output class="error"></output>
          </div>

          <div class='row'>
            <span class='title'>Modification type</span>
            <div class='select-box' id='custom-mod-type'>
              <label for='custom-mod-type-single' tabindex='0'><input type='radio' name='custom-mod-type' value='single' id='custom-mod-type-single' checked>Modification</label>
              <label for='custom-mod-type-linker' tabindex='0'><input type='radio' name='custom-mod-type' value='linker' id='custom-mod-type-linker'>Cross Linker</label>
            </div>
          </div>

          <div class="single">
            <label for="custom-mod-single-specificities">Placement rules</label>
            <div class="list-input single">
              <ul class="values" id="custom-mod-single-specificities"></ul>
              <button class="create" id="custom-mod-single-specificities-create">Create</button>
              <div class="modal" id="custom-mod-single-specificities-create">
                <label for="custom-mod-single-placement-rules">Placement rules</label>
                <div class="separated-input">
                  <div class="values" id="custom-mod-single-placement-rules">
                    <div class="input context" title="Add placement rules with 'AMINOACIDS@Position' or 'Position'. Position can be any of the following: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm." placeholder="Add placement rules with 'AMINOACIDS@Position' or 'Position'" data-type="placement_rule" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <label for="custom-mod-single-neutral-losses">Neutral losses</label>
                <div class="separated-input">
                  <div class="values" id="custom-mod-single-neutral-losses">
                    <div class="input context" placeholder="Add molecular formula or mass preceded by '+' or '-'" data-type="neutral_loss" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <label for="custom-mod-single-diagnostic-ions">Diagnostic ions</label>
                <div class="separated-input">
                  <div class="values" id="custom-mod-single-diagnostic-ions">
                    <div class="input context" placeholder="Add molecular formula or mass" data-type="molecular_formula" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <button class="save" id="custom-mod-single-save">Save</button>
                <button class="cancel secondary" id="custom-mod-single-cancel">Cancel</button>
              </div>
            </div>
          </div>
          <div class="linker">
            <label class="row">Length<input type="number" id="custom-mod-linker-length"></input></label>
            <label for="custom-mod-linker-specificities">Placement rules</label>
            <div class="list-input linker">
              <ul class="values" id="custom-mod-linker-specificities"></ul>
              <button class="create" id="custom-mod-linker-specificities-create">Create</button>
              <div class="modal" id="custom-mod-linker-specificities-create">
                <label class="block">Asymmetric linker<input type="checkbox" switch id="custom-mod-linker-asymmetric"></input></label>
                <label for="custom-mod-linker-placement-rules">Placement rules</label>
                <div class="separated-input">
                  <div class="values" id="custom-mod-linker-placement-rules">
                    <div class="input context" title="Add placement rules with 'AMINOACIDS@Position' or 'Position'. Position can be any of the following: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm." placeholder="Add placement rules with 'AMINOACIDS@Position' or 'Position'" data-type="placement_rule" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <label class="asymmetric" for="custom-mod-linker-placement-rules">Placement rules second position</label>
                <div class="separated-input asymmetric">
                  <div class="values" id="custom-mod-linker-secondary-placement-rules">
                    <div class="input context" title="Add placement rules with 'AMINOACIDS@Position' or 'Position'. Position can be any of the following: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm." placeholder="Add placement rules with 'AMINOACIDS@Position' or 'Position'" data-type="placement_rule" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <label for="custom-mod-linker-stubs">Breakage</label>
                <div class="separated-input">
                  <div class="values" id="custom-mod-linker-stubs">
                    <div class="input context" placeholder="Add breakage with a mass or molecular formula as 'formula1:formula2'." data-type="stub" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <label for="custom-mod-linker-diagnostic-ions">Diagnostic ions</label>
                <div class="separated-input">
                  <div class="values" id="custom-mod-linker-diagnostic-ions">
                    <div class="input context" placeholder="Add molecular formula or mass" data-type="molecular_formula" contentEditable="plaintext-only"></div>
                    <button class="clear">Clear</button>
                  </div>
                  <output class="error"></output>
                </div>
                <button class="save" id="custom-mod-linker-save">Save</button>
                <button class="cancel secondary" id="custom-mod-linker-cancel">Cancel</button>
                </div>
            </div>
          </div>
          <div class="row">
            <button id="custom-mod-save">Save</button>
            <button id="custom-mod-cancel" class="secondary">Cancel</button>
          </div>
        </dialog>
        <button id="custom-mod-create" data-new-id="0">Create new</button>
        <ol id="custom-mods"></ol>
      </fieldset>
      <input type="checkbox" id="collapsible-logs">
      <fieldset class="collapsible" data-linked-item="collapsible-logs">
      <legend>Logs</legend>
      <output class="collapsible-content" id="spectrum-log"></output>
      <output class="collapsible-content" id="identified-peptides-log"></output>
      </fieldset>
      <div class="grow"></div>
      <footer>
        <p>Annotator by Douwe Schulte at the Snijderlab</p>
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
