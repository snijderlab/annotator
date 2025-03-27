use std::io::Write;
use std::{fs::File, io::BufWriter};

#[path = "src/html_builder.rs"]
mod html_builder;

use html_builder::*;

fn create_loss_modal(id: &str) -> HtmlElement {
    let mut dialog = vec![
        HtmlTag::h2.new().content("Normal losses and gains").clone(),
        HtmlTag::div
            .new()
            .class("col-2")
            .children([
                HtmlTag::div
                    .new()
                    .children(
                        std::iter::once(HtmlTag::h3.new().content("Losses").clone()).chain(
                            HtmlElement::input_list(
                                format!("model-{id}-loss-selection"),
                                "checkbox",
                                "block",
                                [
                                    ("-H1O1", "OH"),
                                    ("-H2O1", "Water"),
                                    ("-H4O2", "Double water"),
                                    ("-H6O3", "Triple water"),
                                    ("-H1", "Hydrogen"),
                                    ("-H2", "Double hydrogen"),
                                    ("-H3", "Triple hydrogen"),
                                    ("-H3N1", "Ammonia"),
                                    ("-C1O1", "Carbon monoxide"),
                                ],
                            ),
                        ),
                    )
                    .clone(),
                HtmlTag::div
                    .new()
                    .children(
                        std::iter::once(HtmlTag::h3.new().content("Gains").clone()).chain(
                            HtmlElement::input_list(
                                format!("model-{id}-gain-selection"),
                                "checkbox",
                                "block",
                                [
                                    ("+H2O1", "Water"),
                                    ("+H4O2", "Double water"),
                                    ("+H6O3", "Triple water"),
                                    ("+H1", "Hydrogen"),
                                    ("+H2", "Double hydrogen"),
                                    ("+H3", "Triple hydrogen"),
                                ],
                            ),
                        ),
                    )
                    .clone(),
            ])
            .clone(),
        HtmlTag::h3.new().content("Custom").clone(),
        HtmlElement::separated_input(
            format!("model-{id}-loss"),
            "Add molecular formula or mass preceded by '+' or '-'",
            "neutral_loss",
        ),
    ];
    if id == "glycan" {
        dialog.extend(
        [HtmlTag::h2.new().content("Glycan attachment to fragments").clone(),
            HtmlTag::p.new().content("Glycans can be attached to amino acid in different ways and each of these ways can result in different fragmentation behaviour. Here rules on the behaviour can be defined that determine which part of a glycan can be attached to any fragment with a glycan. Free indicates that the entire glycan is removed. Core indicates that the core, the first HexNAc possibly with some fucoses, remains attached to the fragment. Full indicates that the full glycan is present on all glycan containing fragments. If multiple options are checked all options will be generated.").clone()
          ].into_iter().chain([
            HtmlTag::div.new().class("list-input glycan-fragments").content(
              HtmlTag::ul.new().class("values").id("model-glycan-fragments").content(r#"<li class="element"><span data-value='{"fallback": true, "selection": [], "free": false, "core": false, "full": true }' title="Edit">All undefined, Fragments: full</span></li>"#))
              .content(HtmlTag::button.new().class("create").id("model-glycan-fragments-create").content("New rule"))
              .content(HtmlTag::div.new().class("modal").id("model-glycan-fragments-create").children([
                HtmlTag::label.new().id("model-glycan-fragments-selection-label").header("for", "model-glycan-fragments-selection").content("Amino acid location of a glycan").clone(),
                HtmlElement::separated_input("model-glycan-fragments-selection", "Amino acid code", "amino_acid"),
                HtmlTag::label.new().id("model-glycan-fragments-other").content("All undefined attachment locations").clone(),
                HtmlTag::label.new().header("for", "model-glycan-fragments-form").content("Allowed glycan forms on fragments").clone(),
                HtmlTag::label.new().content(HtmlTag::input.new().header("name", "model-glycan-fragments-form").header("type", "checkbox").id("model-glycan-fragments-free")).content("Free").clone(),
                HtmlTag::label.new().content(HtmlTag::input.new().header("name", "model-glycan-fragments-form").header("type", "checkbox").id("model-glycan-fragments-core")).content("Core").clone(),
                HtmlTag::label.new().content(HtmlTag::input.new().header("name", "model-glycan-fragments-form").header("type", "checkbox").id("model-glycan-fragments-full")).content("Full").clone(),
                HtmlTag::button.new().class("save").id("model-glycan-fragments-save").content("Save").clone(),
                HtmlTag::button.new().class("cancel secondary").id("model-glycan-fragments-cancel").content("Cancel").clone(),
              ])).content(HtmlTag::output.new().class("error").header2("hidden")).clone()
        ]))
    } else {
        dialog.extend(
        [HtmlTag::h2.new().content("Amino acid specific losses and gains").clone(),
          HtmlTag::p.new().content("Some losses are only seen from certain amino acids, so all losses specified here will only be generated if the indicated amino acid is present in the fragment. When defining rules multiple amino acids can be combined with multiple losses, for example <code>N,D:-C1H4,-C1H4O1</code> indicates that N and D can both lose either C<sub>1</sub>H<sub>4</sub> or C<sub>1</sub>H<sub>4</sub>O<sub>1</sub>.").clone()
        ].into_iter().chain(
          HtmlElement::input_list(format!("model-{id}-aa-loss-selection"), "checkbox", "block", [
            ("N:-C1H1O2", "N: COOH"),
            ("Q:-C2H3O2", "Q: C<sub>2</sub>H<sub>3</sub>O<sub>2</sub>"),
          ])).chain([
            HtmlElement::separated_input(format!("model-{id}-aa-loss"), "Amino acid codes followed by a colon and the losses separated by commas", "aa_neutral_loss"),
            HtmlTag::h2.new().content("Amino acid side chain losses").clone(),
            HtmlTag::p.new().content("For some fragmentation techniques, notably ETD, the loss of side chains can be seen. Sometimes even multiple side chains can be lost from a single fragment. The amino acids that can lose a side chain can be selected below, if no selection is given all amino acids will be allowed.").clone(),
            HtmlTag::label.new().content(HtmlTag::input.new().id(format!("model-{id}-aa-side-chain-loss-number")).header("type", "number").header("min", "0").header("value", "0").header("max", "255").clone()).content("Maximal number of lost side chains").clone(),
            HtmlElement::separated_input(format!("model-{id}-aa-side-chain-loss-selection"), "Add amino acid code", "amino_acid"),])
          )
    };
    dialog.push(
        HtmlTag::button
            .new()
            .header2("autofocus")
            .header("onclick", "this.parentElement.close()")
            .content("Close")
            .clone(),
    );
    HtmlTag::div.new().children([
    HtmlTag::button.new().header("onclick", format!("document.getElementById('model-{id}-loss-selection-dialog').showModal();")).content("Select"),
    HtmlTag::output.new().id(format!("model-{id}-loss-selection-output")).class("selected-neutral-loss").content("0 selected"),
    HtmlTag::dialog.new().class("neutral-loss").id(format!("model-{id}-loss-selection-dialog")).header("onclose", format!(r##"var num = 0; document.getElementsByName("model-{id}-loss-selection").forEach(e=>num += e.checked); num += document.querySelectorAll("#model-{id}-loss-selection-dialog .element").length;document.getElementById("model-{id}-loss-selection-output").innerText = num + " selected";"##)).children(dialog)
  ]).clone()
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
            Self::One => ("selected", "", 1, "selected", "", 1),
            Self::OneToPrecursor => ("selected", "", 1, "", "selected", 0),
            Self::Precursor => ("", "selected", 0, "", "selected", 0),
        }
    }
}

fn create_charge_range_fields(id: &str, default: ChargeRange) -> HtmlElement {
    let settings = default.settings();
    HtmlTag::div
        .new()
        .class("charge-range")
        .children([
            HtmlTag::span.new().children([
                HtmlTag::select
                    .new()
                    .id(format!("model-{id}-charge-start-type"))
                    .children([
                        HtmlTag::option
                            .new()
                            .value("Absolute")
                            .title("An absolute charge")
                            .header2(settings.0)
                            .content("Absolute"),
                        HtmlTag::option
                            .new()
                            .value("Relative")
                            .title("Relative to the precursor charge")
                            .header2(settings.1)
                            .content("Precursor"),
                    ]),
                HtmlTag::input
                    .new()
                    .id(format!("model-{id}-charge-start-value"))
                    .header("type", "number")
                    .value(settings.2.to_string()),
            ]),
            HtmlTag::span.new().content("-"),
            HtmlTag::span.new().children([
                HtmlTag::select
                    .new()
                    .id(format!("model-{id}-charge-end-type"))
                    .children([
                        HtmlTag::option
                            .new()
                            .value("Absolute")
                            .title("An absolute charge")
                            .header2(settings.3)
                            .content("Absolute"),
                        HtmlTag::option
                            .new()
                            .value("Relative")
                            .title("Relative to the precursor charge")
                            .header2(settings.4)
                            .content("Precursor"),
                    ]),
                HtmlTag::input
                    .new()
                    .id(format!("model-{id}-charge-end-value"))
                    .header("type", "number")
                    .value(settings.5.to_string()),
            ]),
        ])
        .clone()
}

fn main() {
    let version = option_env!("CARGO_PKG_VERSION").unwrap_or("undefined");
    let file = File::create("../frontend/index.html").unwrap();
    println!("{:?}", file.metadata());
    let mut writer = BufWriter::new(file);
    write!(
            writer,
    r#"<!DOCTYPE html>
    <html lang="en">
    
    <head>
      <meta charset="UTF-8" />
      <link rel="stylesheet" href="styles.css" />
      <meta name="viewport" content="width=device-width, initial-scale=1.0" />
      <title>Annotator</title>
      <script type="module" src="/main.js" defer></script>
      <script src="script.js"></script>
    </head>
    
    <body class="theme-auto" id="body">
      <div class='header'>
        <button print" onclick="window.print()">Export</button>
        <button class="cancel-drop" onclick='document.querySelector("html").classList.remove("file-drop-hover")'>Cancel drop</button>
        <div class="theme">
          <span>Theme</span>
          <div class='select-box' id='theme'>
            <label tabindex='0'><input type='radio' name='theme' value='theme-light' id='theme-light'>Light</label>
            <label tabindex='0'><input type='radio' name='theme' value='theme-auto' id='theme-auto' checked>Auto</label>
            <label tabindex='0'><input type='radio' name='theme' value='theme-dark' id='theme-dark'>Dark</label>
          </div>
        </div>
      </div>
      <div class="input-flex">
        <div class="joined-button" id="load-raw-path"><button type="button" id="load-raw-file">Load raw data file</button><button type="button" id="load-raw-folder" title="Open a Bruker TDF .d directory">folder</button></div>
        <button type="button" id="load-clipboard">Load Clipboard</button>
        <button type="button" id="load-identified-peptides">Load identified peptides file</button>
      </div>
      <div class="input-flex">
        <div class="usi">
          <a target='_blank' href="https://www.psidev.info/usi" title="Universal Spectrum Identifier, a standard to reference any publicly available spectrum">USI:</a>
          <input id="usi" type="text" placeholder="mzspec:ID:FILE:scan:SCAN"></input>
          <button type="button" id="load-usi">Load</button>
        </div>
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
          <option value="all" title="All possible ions with single water loss from all">All</option>
          <option value="ethcd" title="b+c+w+y+z+glycan with single water loss from all and precursor decharging">EThcD/ETcaD</option>
          <option value="hot_eacid" title="a+b+c+d+x+w+y+z+glycan with single water loss from all and precursor decharging">Hot EACID</option>
          <option value="ead" title="a+b+c+d+x+w+y+z+glycan with single water loss from all and precursor decharging">EAD</option>
          <option value="cidhcd" title="a+b+d+y+precursor with single water loss">CID/HCD</option>
          <option value="etd" title="c+y+z with single water loss and precursor with single water, ammonia, and ETD specific losses">ETD</option>
          <option value="td_etd" title="c+z with single water, and ammonia loss and one, two, and three hydrogen gain and precursor with single water, ammonia, and ETD specific losses and one, two, and three hydrogen gain">Top-down ETD</option>
          <option value="none" title="Only the base precursor peak, no losses">None</option>
          <option value="custom">Custom</option>
        </select>"#).unwrap();

    write!(
            writer,
            r#"<label class="wide" for="peptide">Peptide sequence </label>
          <div class="peptide-input wide context" id="peptide" contentEditable="plaintext-only"></div>
          <button id="annotate-button" type="button" class="col-2 center">Annotate</button>
          <button id="save-spectrum" type="button" class="secondary center" title="Save the (merged) selected spectrum, with the noise filter applied.">Save selected spectrum</button>
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
            <label for='peptide-intensities' title='Display the intensities in the peptide support'><input id='peptide-intensities' class='compact' type='checkbox' switch/>Display intensities</label>
    
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

          <fieldset class='settings spectrum-settings peaks-settings'>
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

            <label title='Create distance labels by dragging from one peak to another, delete them by clicking on them, delete all with the button here.'>Distance labels <button id='distance-labels-clear'>Clear</button></label>
          </fieldset>
    
          <fieldset class='settings spectrum-settings'>
            <legend>Label settings</legend>
            <label class='has-slider label row align'>
              <span>Show labels for top:</span>
              <input id='spectrum-label' type='range' min='0' max='100' value='90'/>
              <input id='spectrum-label-value' type='number' min='0' max='100' value='90'/>
              %
            </label>
            
            <label class='has-slider mz row align'>
              <span>Show m/z for top:</span>
              <input id='spectrum-m-z' type='range' min='0' max='100' value='0'/>
              <input id='spectrum-m-z-value' type='number' min='0' max='100' value='0'/>
              %
            </label>

            <label class='has-slider mz row align'>
              <span>Show glycan for top:</span>
              <input id='spectrum-glycan' type='range' min='0' max='100' value='90'/>
              <input id='spectrum-glycan-value' type='number' min='0' max='100' value='90'/>
              %
            </label>
    
            <div class='row'>
              <span class='title'>Manually force show</span>
              <div class='select-box' id='force-show' style='min-width: fit-content;'>
                <label tabindex='0'><input type='radio' name='force-show' value='none' id='force-show-none' checked>None</label>
                <label tabindex='0'><input type='radio' name='force-show' value='label' id='force-show-label'>Label</label>
                <label tabindex='0'><input type='radio' name='force-show' value='m-z' id='force-show-m-z'>m/z</label>
                <label tabindex='0'><input type='radio' name='force-show' value='glycan' id='force-show-glycan'>Glycan</label>
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
        <h2>Search for modifications</h2>
        <p>Search by mass with a given tolerance, by formula (e.g. <code>Formula:O-1</code>), or by glycan composition  (e.g. <code>Glycan:Hex4HexNAc2</code>) to search for matching modification. Additionally searching for a name or index for an existing modification will show all details of that modification.</p>
        <div class="flex-input collapsible-content">
          <label for="search-modification">Search modification</label>
          <input id="search-modification" type="text" title="Search for a ProForma modification, it will display its details. If you give a mass, formula, or glycan composition modification it instead will return all modifications with that composition."></input>
          <label for="search-modification-tolerance" title="If searching for a mass modification, the tolerance on the mass">Tolerance (in Da)</label>
          <input id="search-modification-tolerance" type="number" value="0.1" min="0"></input>
          <button id="search-modification-button">Search</button>
        </div>
        <output class="error collapsible-content hidden" id="search-modification-error"></output>
        <output class="collapsible-content" id="search-modification-result"></output>

        <h2>Isotopic distribution</h2>
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
              <label for='custom-mod-type-linker' tabindex='0'><input type='radio' name='custom-mod-type' value='linker' id='custom-mod-type-linker'>Cross-linker</label>
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
              <ouput class="error" hidden></output>
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
                <label for="custom-mod-linker-stubs" title="Define the ways in which this crosslinker can break apart. These formulas are used as the full formula when breakage is allowed, so the original cross-linker definition is not used in this case. For example a crosslinker of formula 'Ca+bHa+bOa+b' which can break into 'CaHaOa' and 'CbHbOb' has to be defined as 'CaHaOa:CbHbOb'. Commonly the summation of the two halves after breaking will result in the same formual as the full cross-linker, but does not have to be the case.">Breakage</label>
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
              <ouput class="error" hidden></output>
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
      <input type="checkbox" id="collapsible-custom-models">
      <fieldset class="collapsible" data-linked-item="collapsible-custom-models" id="custom-models-collapsible">
        <legend>Custom models</legend>
        <p>Path to configuration file: <span style='-webkit-user-select:all;user-select:all;' id='custom-models-path'>Not loaded</span></p>
        <dialog class="custom-model" id="custom-model-dialog">
          <legend>Custom model</legend>
          <p>Ion</p>
          <p>Location</p>
          <p>Neutral Loss/Gain</p>
          <p>Charge range</p>
          <p>Variant</p>"#).unwrap();
    for ion in ["a", "b", "c", "d", "v", "w", "x", "y", "z"] {
        write!(writer, "<label>{ion}</label>").unwrap();
        if ["d", "v", "w"].contains(&ion) {
            let other_ion = match ion {
                "d" => "a",
                "v" => "y",
                "w" => "z",
                _ => "",
            };
            write!(writer, r#"<div class="satellite-ion-location" title="{ion} ions are formed by secondary fragmentation from {other_ion}, so if turned on these will be generated for any location that produces {other_ion}. The maximal distance is the maximal number of side chains between the fragmenting side chain and the parent ion breakage, when the side chain immediatly adjecent to the parent ion breakage breaks this distance is 0, so to only get standard satellite ions the distance needs to be set to 0. If no distance is specified the satellite ions are not generated.">{}<label>Base maximal distance<input type="number" min="0" max="255" value="" id="model-{ion}-base-distance"></input></label></div>"#,
              HtmlElement::separated_input(format!("model-{ion}-location"), "Amino acid codes followed by a colon and the maximal distance", "satellite_ion"),
            ).unwrap();
        } else {
            write!(
              writer,
             r#"<div id="model-{ion}-location" class="location select-input">
                <select onchange="this.className=this.options[Number(this.value)].dataset.cls;">
                  <option value="0" data-cls="arg-0" data-value="All" title="All backbone bonds produce fragments ions">All</option>
                  <option value="1" data-cls="arg-0" data-value="None" selected title="No fragments are generated">None</option>
                  <option value="2" data-cls="arg-1" data-value="SkipN" title="Select a number of amino acids from the N terminal that do not produce a fragment, the rest does produce fragments.">Disallow x from N terminal</option>
                  <option value="3" data-cls="arg-1" data-value="SkipC" title="Select a number of amino acids from the C terminal that do not produce a fragment, the rest does produce fragments.">Disallow x from C terminal</option>
                  <option value="4" data-cls="arg-1" data-value="TakeN1" title="Select a number of amino acids from the N terminal that do produce a fragment, the rest does not produce fragments.">Allow x from N terminal</option>
                  <option value="5" data-cls="arg-1" data-value="TakeC" title="Select a number of amino acids from the C terminal that do produce a fragment, the rest does not produce fragments.">Allow x from C terminal</option>
                  <option value="6" data-cls="arg-2" data-value="TakeN" title="Select an offset from the N terminal that do not produce fragments, then select a number of amino acids that do.">Disallow x from N and allow y</option>
                  <option value="7" data-cls="arg-2" data-value="SkipNC" title="Select an offset from the N terminal that do not produce fragments, and select an offset from the C terminal that does not produce fragments, the middle left over section does">Disallow x from N and disallow y from C</option>
                </select>
                <input type="number" value="1" min="1">
                <input type="number" value="1" min="1">
                </div>"#
          ).unwrap();
        }
        write!(writer, "{}", create_loss_modal(ion)).unwrap();
        write!(
            writer,
            "{}",
            create_charge_range_fields(ion, ChargeRange::OneToPrecursor)
        )
        .unwrap();
        write!(
                writer,
                r#"<div class="variant" id="variant-{ion}">
                  <label title="The normal {ion} ion with two less hydrogens"><input type="checkbox" id="variant-{ion}-2"></input>''</label>
                  <label title="The normal {ion} ion with one less hydrogen"><input type="checkbox" id="variant-{ion}-1"></input>'</label>
                  <label title="The normal definition of the {ion} ion as defined by Biemann"><input type="checkbox" id="variant-{ion}0" checked></input>{ion}</label>
                  <label title="The normal {ion} ion with one more hydrogen"><input type="checkbox" id="variant-{ion}+1"></input>·</label>
                  <label title="The normal {ion} ion with two more hydrogens"><input type="checkbox" id="variant-{ion}+2"></input>··</label>
                </div>"#,
              )
              .unwrap();
    }
    write!(writer, r#"<label>precursor</label>
          <div class="empty"></div>
          {}{}
          <div class="grid-row">
            <label>glycan</label>
            <label><input id='model-glycan-enabled' type='checkbox' switch/>Enable fragments from structure (GNO)</label>
            {}{}
            <span>(Y)</span>
          </div>
          <div class='grid-row'>
            <label>glycan</label>
            <span style="grid-column: span 2">Enable fragments from compositions between <input id='model-glycan-composition-min' type='number' min='0' value='0'/> — <input id='model-glycan-composition-max' type='number' min='0' value='0'/> monosaccharides</span>
            {}
            <span>(B)</span>
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
            <label title="Allow modification specific neutral losses, as defined by the database">modification neutral losses</label>
            <label><input id='model-modification-neutral-enabled' type='checkbox' switch/>Enable</label>
          </div>
          <div class="grid-row">
            <label title="Allow MS cleavable cross-links to be cleaved">MS cleavable cross-links</label>
            <label><input id='model-cleave-cross-links-enabled' type='checkbox' switch/>Enable</label>
          </div>
        </dialog>
        <button id="custom-model-create">Create new</button>
        <ul id="custom-models"></ul>
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
    create_charge_range_fields("precursor", ChargeRange::Precursor),
    create_loss_modal("glycan"),
    create_charge_range_fields("glycan-other", ChargeRange::OneToPrecursor),
    create_charge_range_fields("glycan-oxonium", ChargeRange::One),
    create_charge_range_fields("diagnostic", ChargeRange::One),
    create_charge_range_fields("immonium", ChargeRange::One),
    version
        )
        .unwrap();
    tauri_build::build()
}
