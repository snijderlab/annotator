use std::io::Write;
use std::{fs::File, io::BufWriter};

fn create_loss_modal(id: &str) -> String {
    format!(
        r##"<div>
            <button onclick='document.getElementById("model-{id}-loss-selection-dialog").showModal();'>Select</button>
            <output id='model-{id}-loss-selection-output' class='selected-neutral-loss'>0 selected</output>
            <dialog class="neutral-loss" id="model-{id}-loss-selection-dialog" onclose='var num = 0; document.getElementsByName("model-{id}-loss-selection").forEach(e=>num += e.checked); num += document.querySelectorAll("#model-{id}-loss-selection-dialog .element").length;document.getElementById("model-{id}-loss-selection-output").innerText = num + " selected";'>
              <p>Losses</p>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H1O1"/>OH</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H2O1"/>Water</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H4O2"/>Double water</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H6O3"/>Triple water</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H1"/>Hydrogen</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H2"/>Double hydrogen</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H3"/>Triple hydrogen</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-H3N1"/>Ammonia</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-C1O1"/>Carbon monoxide</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-C1H1O2"/>COOH (seen with ETD on Asp)</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="-C2H3O2"/>C<sub>2</sub>H<sub>3</sub>O<sub>2</sub> (seen with ETD on Glu)</label>
              <p>Gains</p>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H2O1"/>Water</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H4O2"/>Double water</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H6O3"/>Triple water</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H1"/>Hydrogen</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H2"/>Double hydrogen</label>
              <label class='block'><input type="checkbox" name="model-{id}-loss-selection" value="+H3"/>Triple hydrogen</label>
              <p>Custom</p>
              <div class="separated-input">
                <div class="values" id="model-{id}-loss">
                  <div class="input context" placeholder="Add molecular formula or mass preceded by '+' or '-'" data-type="neutral_loss" contentEditable="plaintext-only"></div>
                  <button class="clear">Clear</button>
                </div>
                <output class="error"></output>
              </div>
              <button autofocus onclick='this.parentElement.close()'>Close</button>
            </dialog>
          </div>"##
    )
}

#[derive(Copy, Clone)]
enum ChargeRange {
    One,
    OneToPrecursor,
    Precursor,
}

impl ChargeRange {
    fn settings(
        self,
    ) -> (
        &'static str,
        &'static str,
        usize,
        &'static str,
        &'static str,
        usize,
    ) {
        match self {
            Self::One => (" selected", "", 1, " selected", "", 1),
            Self::OneToPrecursor => (" selected", "", 1, "", " selected", 0),
            Self::Precursor => ("", " selected", 0, "", " selected", 0),
        }
    }
}

fn create_charge_range_fields(id: &str, default: ChargeRange, comment: &str) -> String {
    let settings = default.settings();
    format!(
        r##"<div class="charge-range">
          <select id="model-{id}-charge-start-type">
            <option value="Absolute" title="An absolute charge"{}>Absolute</option>
            <option value="Relative" title="Relative to the precursor charge"{}>Precursor</option>
          </select>
          <input id="model-{id}-charge-start-value" type="number" value="{}">
          <span>—</span>
          <select id="model-{id}-charge-end-type">
            <option value="Absolute" title="An absolute charge"{}>Absolute</option>
            <option value="Relative" title="Relative to the precursor charge"{}>Precursor</option>
          </select>
          <input id="model-{id}-charge-end-value" type="number" value="{}">
          {}
        </div>"##,
        settings.0,
        settings.1,
        settings.2,
        settings.3,
        settings.4,
        settings.5,
        if comment.is_empty() {
            String::new()
        } else {
            format!("<span> ({comment})</span>")
        },
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
        <button type="button" id="load-raw-path">Load raw data file</button>
        <button type="button" id="load-clipboard">Load Clipboard</button>
        <button type="button" id="load-identified-peptides">Load identified peptides file</button>
      </div>
      <output class="wrap" id="open-files-error"></output>
      
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
              <label for="details-identified-peptide-files">File</label>
              <select id="details-identified-peptide-files"></select>
              <button id="close-identified-peptide-file" type="button">Close file</button>
              <label for="details-identified-peptide-index">Peptide index</label>
              <div class="combined-input">
                <input type="number" id="details-identified-peptide-index" value="0" min="0" />
                <span>/</span>
                <span id="number-of-identified-peptides">0</span>
              </div>
              <button id="load-identified-peptide" title="Find this peptide in the raw data file and populate all setting fields with this peptide's data" type="button">Load</button>
            </div>
            <div id="identified-peptide-details"></div>
          </div>
        </div>
      </div>

      <div>
        <h2>Spectrum details</h2>
        <ol id="spectra"></ol>
      </div>
      <p class="wrap" id="spectrum-details"></p>

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
        
        <label for="model-mz-range-min">m/z range</label>
        <div class="row">
          <input style='flex-grow:1' type="number" id="model-mz-range-min" min="0" value="" placeholder="Empty imposes no bounds" />
          <span style='padding: 0 .5em'>—</span>
          <input style='flex-grow:1' type="number" id="model-mz-range-max" min="0" value="" placeholder="Empty imposes no bounds" />
        </div>
        
        <label for="noise-filter" title="Determine the noise level from the spectrum and remove everything below this factor times the noise level">Noise filter</label>
        <input id="noise-filter" type="number" value="1.0" min="0.0">
    
        <label for="spectrum-model">Model </label>
        <select id="spectrum-model">
          <option value="all" title="All possible ions with single water loss from all and additionally double water loss from glycans">All</option>
          <option value="ethcd" title="b+c+w+y+z+glycan with single water loss from all and additionally double water loss from glycans">EThcD/ETcaD</option>
          <option value="cidhcd" title="a+b+d+y+precursor with single water loss">CID/HCD</option>
          <option value="etd" title="c+y+z with single water loss and precursor with single water, ammonia, and ETD specific losses">ETD</option>
          <option value="none" title="Only the base precursor peak, no losses">None</option>
          <option value="custom">Custom</option>
        </select>
        <fieldset class="custom-model">
          <legend>Custom model</legend>
          <p>Ion</p>
          <p>Location</p>
          <p>Neutral Loss/Gain</p>
          <p>Charge range</p>"#).unwrap();
    for ion in ["a", "b", "c", "d", "v", "w", "x", "y", "z"] {
        write!(
                writer,
                r#"<label>{0}</label>
              <div id="model-{0}-location" class="location select-input">
              <select onchange="this.className=this.options[Number(this.value)].dataset.cls;">
                <option value="0" data-cls="arg-0" data-value="All" title="All backbone bonds produce fragments ions">All</option>
                <option value="1" data-cls="arg-0" data-value="None" selected title="No fragments are generated">None</option>
                <option value="2" data-cls="arg-1" data-value="SkipN" title="Select a number of amino acids from the N terminal that do not produce a fragment, the rest does produce fragments.">Disallow x from N terminal</option>
                <option value="3" data-cls="arg-1" data-value="SkipC" title="Select a number of amino acids from the C terminal that do not produce a fragment, the rest does produce fragments.">Disallow x from C terminal</option>
                <option value="4" data-cls="arg-1" data-value="TakeN" title="Select a number of amino acids from the N terminal that do produce a fragment, the rest does not produce fragments.">Allow x from N terminal</option>
                <option value="5" data-cls="arg-1" data-value="TakeC" title="Select a number of amino acids from the C terminal that do produce a fragment, the rest does not produce fragments.">Allow x from C terminal</option>
                <option value="6" data-cls="arg-2" data-value="SkipNC" title="Select an offset from the N terminal that do not produce fragments, then select a number of amino acids that do.">Disallow x from N and allow y</option>
              </select>
              <input type="number" value="1" min="1">
              <input type="number" value="1" min="1">
              </div>
              {1}{2}"#,
                ion, create_loss_modal(ion), create_charge_range_fields(ion, ChargeRange::OneToPrecursor, "")
            )
            .unwrap();
    }
    write!(
            writer,
            r#"<label>precursor</label>
            <div class="empty"></div>
            {}{}
            <div class="grid-row">
              <label>glycan</label>
              <label><input id='model-glycan-enabled' type='checkbox' switch/>Enable fragments from structure (GNO)</label>
              {}{}
            </div>
            <div class='grid-row'>
              <label>glycan</label>
              <span style="grid-column: span 2">Enable fragments from compositions between <input id='model-glycan-composition-min' type='number' min='0' value='0'/> — <input id='model-glycan-composition-max' type='number' min='0' value='0'/> monosaccharides</span>
              {}
            </div>
            <div class="grid-row">
              <label title="Allow modification specific diagnostic ions, as defined by the database">modification diagnostic ions</label>
              <label><input id='model-modification-diagnostic-enabled' type='checkbox' switch/>Enable</label>
              <span class='empty'></span>
              {}
            </div>
            <div class="grid-row">
              <label>immonium</label>
              <label><input id='model-immonium-enabled' type='checkbox' switch/>Enable</label>
              <span class='empty'></span>
              {}
            </div>
            <div class="grid-row">
              <label title="side chain loss from precursor as seen in electron based fragmentation">precursor side chain loss</label>
              <label><input id='model-m-enabled' type='checkbox' switch/>Enable</label>
            </div>
            <div class="grid-row">
              <label title="Allow modification specific neutral losses, as defined by the database">modification neutral losses</label>
              <label><input id='model-modification-neutral-enabled' type='checkbox' switch/>Enable</label>
            </div>
            <div class="grid-row">
              <label title="Allow MS cleavable cross-links to be cleaved">MS cleavable cross-links</label>
              <label><input id='model-cleave-cross-links-enabled' type='checkbox' switch/>Enable</label>
            </div>
          </fieldset>
          <label class="wide" for="peptide">Peptide sequence </label>
          <div class="peptide-input wide context" id="peptide" contentEditable="plaintext-only"></div>
          <button id="annotate-button" type="button" class="col-2 center">Annotate</button>
          <svg class='orbitrap' width='100' height='50' version='1.1' viewBox='0 0 26.458 13.229' xmlns='http://www.w3.org/2000/svg'>
            <path class='base' d='m26.458 7.4586c-3.299-0.010534-4.2389 0.22021-6.6141 1.4114-2.5075 1.2575-3.9311 1.4114-6.6141 1.4114s-4.1066-0.1538-6.6141-1.4114c-2.3751-1.1911-3.3151-1.4219-6.6141-1.4114m26.456-1.6879c-3.299 0.010534-4.2389-0.2202-6.6141-1.4114-2.5075-1.2575-3.9311-1.4114-6.6141-1.4114s-4.1066 0.1538-6.6141 1.4114c-2.3751 1.1911-3.3151 1.4219-6.6141 1.4114' />
            <path class='trace' d='m17.074 9.9441c0.52615 1.7417 1.1464 3.0166 2.0544 3.0166 1.1421-3e-6 1.8118-2.1616 2.0104-4.7007 0.08406-1.0746 0.08373-2.2168-8.64e-4 -3.2913-0.19988-2.5389-0.87019-4.7001-2.0096-4.7001-0.90806-3.63e-6 -1.5283 1.2749-2.0544 3.0166m-1.8388 6.9381c-0.55028 1.5971-1.226 2.7374-2.1908 2.7374-0.97053-2e-6 -1.6482-1.1555-2.2002-2.7691-0.68781-2.0107-1.1804-4.733-1.7569-6.6688-0.57644-1.9358-1.2028-3.2543-2.1269-3.2543-0.86681 3.63e-6 -1.3781 2.1348-1.534 4.6535m0 3.3852c0.15581 2.5186 0.66712 4.6535 1.534 4.6535 0.92431-2e-6 1.5507-1.3188 2.0829-3.1062 0.62072-2.0847 1.1133-4.807 1.8012-6.8176 0.55194-1.6132 1.2296-2.7684 2.1999-2.7684 0.96475 3.63e-6 1.6403 1.1401 2.1906 2.7369' />
          </svg>
        </div>
      <output id="spectrum-error" class="hidden error"></output>
      <div id='spectrum-wrapper' class="spectrum show-unassigned legend-ion hidden show-charge show-series show-glycan-id show-peptide-id show-neutral-losses show-cross-links show-ambiguous-amino-acids show-modifications" onload='SpectrumSetUp()'>
        <div class='legend'>
          <span class='title'>Ion legend</span>
          <div class='ion-series'>
              <div class='top'>
                  <span class='ion w' tabindex='0'>w</span>
                  <span class='ion v' tabindex='0'>v</span>
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
              </div>
              <div class='side'>
                  <span class='ion p' tabindex='0'>Precursor</span>
                  <span class='ion multi mp mpp' tabindex='0'>Multi</span>
                  <span class='ion oxonium' tabindex='0'>Oxonium</span>
                  <span class='ion Y' tabindex='0'>Y</span>
                  <span class='ion neutral-loss' tabindex='0'>Neutral loss</span>
                  <span class='ion diagnostic' tabindex='0'>Diagnostic ions</span>
                  <span class='other'>Other</span>
              </div>
          </div>
        </div>
        <output class='wrapper show-unassigned' id="spectrum-results-wrapper"></output>

        <input type="checkbox" id="collapsible-settings">
        <fieldset class="collapsible all-settings" data-linked-item="collapsible-settings" id="settings">
          <legend>Settings</legend>

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

          <fieldset class='settings peaks-settings'>
            <legend>Peaks settings</legend>
            
            <label for='theoretical' title='Show the theoretical peptide spectrum on the x axis'><input id='theoretical' class='theoretical' type='checkbox' switch/>Show theoretical spectrum</label>
            
            <label for='unassigned' title='Show the unassigned peaks in the spectrum'><input id='unassigned' class='unassigned' type='checkbox' switch checked/>Show unassigned peaks</label>
    
            <div class='row'>
              <span class='title'>Ion colour mode</span>
              <div class='select-box' id='peak-colour'>
                <label for='peak-colour-ion' tabindex='0'><input type='radio' name='peak-colour' value='ion' id='peak-colour-ion' checked>Ion</label>
                <label for='peak-colour-peptide' tabindex='0'><input type='radio' name='peak-colour' value='peptide' id='peak-colour-peptide'>Peptide</label>
                <label for='peak-colour-peptidoform' tabindex='0'><input type='radio' name='peak-colour' value='peptidoform' id='peak-colour-peptidoform'>Peptidoform</label>
                <label for='peak-colour-none' tabindex='0'><input type='radio' name='peak-colour' value='none' id='peak-colour-none'>None</label>
              </div>
            </div>

            <label>Measure distance <button id='peak-distance'>Use ruler</button></label>
          </fieldset>
    
          <fieldset class='settings spectrum-settings'>
            <legend>Label settings</legend>
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
              <div class='select-box' id='force-show' style='min-width: fit-content;'>
                <label tabindex='0'><input type='radio' name='force-show' value='none' id='force-show-none' checked>None</label>
                <label tabindex='0'><input type='radio' name='force-show' value='label' id='force-show-label'>Label</label>
                <label tabindex='0'><input type='radio' name='force-show' value='m-z' id='force-show-m-z'>m/z</label>
                <label tabindex='0'><input type='radio' name='force-show' value='hide' id='force-show-hide' title='Select a peak to not show label & m/z regardless of the show range selected'>Hide</label>
                <!--<label tabindex='0'><input type='radio' name='force-show' value='distance' id='force-show-distance'>Distance</label>-->
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
                  <label><input checked id='spectrum-label-cross-links' type='checkbox'/>Cross-links</label>
                  <label><input checked id='spectrum-label-ambiguous-amino-acids' type='checkbox'/>Ambiguous amino acids</label>
                  <label><input checked id='spectrum-label-modifications' type='checkbox'/>Modifications of unknown position</label>
                  <label><input id='spectrum-label-charge-carriers' type='checkbox'/>Charge carriers</label>
                </div>
            </div>

            <label><input id='rotate-label' class='rotate-label' type='checkbox' switch/>Rotate labels</label>
          </fieldset>
          
          <fieldset class='settings spectrum-settings'>
            <legend>Spectrum settings</legend>
            
            <div>
              <span>Mz</span>
              <label for='spectrum-mz-min'>min</label>
              <input id='spectrum-mz-min' class='mz-min' type='number' value='0' min='0'/>
              <label for='spectrum-mz-max'>max</label>
              <input id='spectrum-mz-max' class='mz-max' type='number' value='42' min='0'/>
            </div>
    
            <div>
              <span>Intensity</span>
              <label for='spectrum-intensity-max'>max</label>
              <input id='spectrum-intensity-max' class='intensity-max' type='number' value='42' min='0'/>
            </div>

            <div>
              <span>Tick marks</span>
              <label for='spectrum-ticks-x'>X</label>
              <input id='spectrum-ticks-x' type='number' value='5' min='2'/>
              <label for='spectrum-ticks-y'>Y</label>
              <input id='spectrum-ticks-y' type='number' value='5' min='2'/>
            </div>

            <label title='Show the intensity in square root, this emphasizes the lower intensity peaks'><input id='y-sqrt' class='y-sqrt' type='checkbox' switch/>Square root intensity</label>
            <label><input id='y-percentage' class='y-percentage' type='checkbox' switch/>Intensity percent</label>
    
            <button id='reset-zoom' class='reset-zoom' title='Reset the zoom to the default' tabindex='0'>Reset zoom</button>
          </fieldset>
        </fieldset>

        <input type="checkbox" id="collapsible-fragment-table">
        <fieldset class="collapsible" data-linked-item="collapsible-fragment-table">
          <legend>Fragment Table</legend>
          <output class="collapsible-content" id="spectrum-fragment-table"></output>
        </fieldset>
      </div>
    
      <input type="checkbox" id="collapsible-tools">
      <fieldset class="collapsible" data-linked-item="collapsible-tools" id="tools">
        <legend>Tools</legend>
        <div class="flex-input collapsible-content">
          <label for="search-modification">Search modification</label>
          <input id="search-modification" type="text" title="Search for a ProForma modification, it will display its details. If you give a mass, formula, or glycan composition modification it instead will return all modifications with that composition."></input>
          <label for="search-modification-tolerance" title="If searching for a mass modification, the tolerance on the mass">Tolerance (in Da)</label>
          <input id="search-modification-tolerance" type="number" value="0.1" min="0"></input>
          <button id="search-modification-button">Search</button>
        </div>
        <output class="error collapsible-content hidden" id="search-modification-error"></output>
        <output class="collapsible-content" id="search-modification-result"></output>
        <div class="flex-input collapsible-content">
          <label for="details-formula">Isotopic distribution for formula</label>
          <span class="input context" id="details-formula" type="text" placeholder="Molecular formula" contentEditable="plaintext-only"></span>
        </div>
        <output class="error collapsible-content hidden" id="details-formula-error"></output>
        <output class="collapsible-content" id="details-formula-result"></output>
      </fieldset>
      <input type="checkbox" id="collapsible-custom-mods">
      <fieldset class="collapsible" data-linked-item="collapsible-custom-mods" id="custom-modifications">
        <legend>Custom modifications</legend>
        <p>Path to configuration file: <span style='-webkit-user-select:all;user-select:all;' id='custom-modifications-path'>Not loaded</span></p>
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
                    <div class="input context" placeholder="Add molecular formula or mass (M not MH+)" data-type="molecular_formula" contentEditable="plaintext-only"></div>
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
                    <div class="input context" placeholder="Add molecular formula or mass (M not MH+)" data-type="molecular_formula" contentEditable="plaintext-only"></div>
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
      <div class="grow"></div>
      <footer>
        <p>Annotator by Douwe Schulte at the Snijderlab</p>
        <a target='_blank' href='https://github.com/snijderlab/annotator'>Open source at github (bugs can be reported here)</a>
        <p>Licensed under MIT OR Apache-2.0</p>
        <p>Version {}</p>
      </footer>
    </body>
    
    </html>"#,
    create_loss_modal("precursor"),
    create_charge_range_fields("precursor", ChargeRange::Precursor, ""),
    create_loss_modal("glycan"),
    create_charge_range_fields("glycan-other", ChargeRange::OneToPrecursor, "Y"),
    create_charge_range_fields("glycan-oxonium", ChargeRange::One, "Oxonium"),
    create_charge_range_fields("diagnostic", ChargeRange::One, ""),
    create_charge_range_fields("immonium", ChargeRange::One, ""),
    version
        )
        .unwrap();
    tauri_build::build()
}
