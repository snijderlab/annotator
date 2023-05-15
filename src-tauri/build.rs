use std::io::Write;
use std::time::Instant;
use std::{fs::File, io::BufWriter};

fn main() {
    let file = File::create("../src/index.html").unwrap();
    println!("{:?}", file.metadata());
    let mut writer = BufWriter::new(file);
    let now = Instant::now();
    write!(
        writer,
r#"<!DOCTYPE html>
<html lang="en">

<head>
  <meta charset="UTF-8" />
  <link rel="stylesheet" href="stitch-assets/styles.css" />
  <link rel="stylesheet" href="style.css" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Stitch+Ox</title>
  <script type="module" src="/main.js" defer></script>
  <script src="stitch-assets/script.js"></script>
  <script src="/basic.js" defer></script>
</head>

<body>
<p>Made at {now:?}</p>
  <details>
    <summary>Align</summary>
    <div class="input-settings">
      <h2>Load reads from structure file</h2>
      <input id="load-path" placeholder="Enter path" value="C:\Users\douwe\src\pdbtbx\example-pdbs\1yyf.pdb" />
      <input id="load-min-length" type="number" value="5" />
      <button id="load-button" type="button">Load</button>
      <pre id="error-log"></pre>
      <h2>Settings for alignment</h2>
      <input id="greet-input-a" placeholder="Enter template sequence"
        value="EVQLVESGGGLVQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARXXXXXXXXXXXXXXXXXXXX" />
      <textarea id="greet-input-b" placeholder="Enter reads (each on their own line)">
DLQLVESGGGLVGAKSPPGTLSAAASGFNL
DLQLVESGGGLVGAKSPPGTLSAAASGFNL
EVQLVESGGGLVQPGGSLSGAKYHSGFNL
EVVQLVESGGGLVQPGGSLGVLSCAASGF
DLQLVESGGGLVQPGGSLGVLSCAASGF
DLQLVESGGGLVQPGTPLYWNAASGFNL
DLQLVESGGGLVQPGGSLRLSCAASGF
QVQLVESGGGLVQPGGSLRLSCAASGF
EVQLVESGGGLPVQGGSLRLSCAADGF
EVQLVESGGGLVQPGGSLRLSCAASGF
EVQLVSGEGGLVQPGGSLRLSCAASGF
QVELVESGGGLVQPGGSLRLSCAASGF
TLSADTSKNTAYLQMNSLRAEDTAVY
RFTLSADTSKNTAYLQMNSLRAEDTA
QLVESGGGLVQPGGSLTHVAGAGHSGF
SADTSKNTAYLQMNSLRAEDTAVYY
LMLTDGYTRYADSVKGRFTLSADTS
QLVESGGGLVQPGGSLRLSCAASGF
QLVESGGGLVQPGGSLRLSCQTGF
LVESGGGLVQPNSLRLSCAASGF
      </textarea>
      <select id="input-alignment-type">
        <option value="1">Local</option>
        <option value="2">Global for B</option>
        <option value="3">Global</option>
      </select>
      <button id="greet-button" type="button">Align</button>
    </div>

    <div class="results">
      <h2>Resulting alignment</h2>
      <div class="alignment">
        <div class="alignment-wrapper" id="reads-alignment">
        </div>
      </div>
    </div>
  </details>
  <details open>
    <div class="input-settings">
      <summary>Spectra</summary>
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
      <label for="peptide">Peptide sequence </label>
      <input id="peptide" value="VAEINPSNGGTTFNEKFKGGKATJ" />
      <label for="mass-system">Mass type </label>
      <select id="mass-system">
        <option value="monoisotopic">MonoIsotopic</option>
        <option value="averageweight">AverageWeight</option>
        <option value="hecklib">Hecklib</option>
      </select>
      <label for="spectrum-model">Model </label>
      <select id="spectrum-model">
        <option value="all">All</option>
        <option value="ethcd">Ethcd</option>
        <option value="etcid">Etcid</option>
        <option value="cidhcd">CidHcd</option>
      </select>
      <label for="spectrum-charge">Max charge </label>
      <input type="number" id="spectrum-charge" value="" placeholder="Empty takes peptide charge from raw data" />
      <button id="annotate-button" type="button">Annotate</button>
    </div>
    <div id="spectrum-results-wrapper"></div>
    <details>
      <summary>Logs</summary>
      <pre id="spectrum-fragments"></pre>
      <pre id="spectrum-error-log"></pre>
    </details>
  </details>
</body>

</html>"#
    ).unwrap();
    tauri_build::build()
}
