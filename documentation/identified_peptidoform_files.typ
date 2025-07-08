#import "functions.typ": aside, button

= Identified peptidoforms files

There are quite some programs that export peptidoforms with metadata based on spectral data. The annotator supports a selection of database matching and _de novo_ software file formats. To load these files use the button #button[Load identified peptidoforms file] or drag in a file from the file explorer. Any number of files can be loaded at the same time. See the table below for supported formats. Opened identified peptidoforms files can be closed with #button[Close file] in the peptidoforms details pane.

#table(
  columns: (auto, auto),
  table.header([*Open format*], [*Comment*]),
  [Fasta],
  [No header requirements],
  [mzTab],
  [v1.0 & Casanovo],
  [SSL],
  [Spectrum sequence list]
)

#table(
  columns: (auto, auto),
  table.header([*Software*], [*Versions/Formats*]),
  [Basic CSV],
  [A CSV file with the following columns: 'raw_file', 'sequence' (in ProForma format), 'z', 'scan_index', and possibly 'mode' with the fragmentation mode, this ignores any other columns.],
  [DeepNovoFamily],
  [DeepNovo, PointNovo, BiatNovo, PGPointNovo],
  [InstaNovo],
  [1.0.0],
  [MaxQuant],
  [msms, msms scans, novo msms scans, & silac],
  [MSFragger],
  [4.2 Fragpipe: 20, 21, & 22, & Philosopher],
  [NovoB],
  [0.0.1],
  [Novor],
  [Denovo and PSM],
  [OPair],
  [common version],
  [Peaks],
  [X, X+, 11, 12, 13 Dia de novo, Ab, DB peptidoforms, DB PSM, & DB protein-peptidoform],
  [PepNet],
  [1.0],
  [PowerNovo],
  [1.0.1],
  [pLink],
  [2.3],
  [PLGS],
  [3.0],
  [Sage],
  [0.14],
)

#aside[If you have data from an unsupported version of a supported program please open an issue on GitHub and give an example file so that the support can be extended.]

#aside[If you have data from an unsupported program that you think should be supported open an issue on GitHub to discuss.]

== Peptidoform details

Once a file is opened the peptidoform details pane opens. Select the right identified peptidoforms file from the drop down menu and select the right peptidoform by peptidoform index. Peptidoforms can be browsed efficiently by selecting the peptidoform index box and using the arrow up/down keys.

Once a peptidoform is selected it will show common metadata in a structured format, followed by an overview of the peptidoform, with the local confidence (if present in the file) depicted by blue squares, and modifications depicted by blue dots. Hover over a modified amino acid to see the full definition. Terminal modifications will be depicted as blue dots on their respective terminal.

Below the structured metadata follows in table format all other metadata from the file format.

Use #button[Load] to select the right spectrum (if the corresponding raw file is open) and load the details for annotation. This will load the sequence, charge, and method if these are available.

== Peptidoform search 

Once at least one identified peptidoforms file is open all peptidoforms can be searched for sequence patterns. Type the search pattern as a ProForma (see @proforma) peptidoform in the search box and hit #button[Search] to search. The search is based on #link("https://doi.org/10.1021/acs.jproteome.4c00188")[mass based alignment] so any peptidoform matching the mass pattern of the search will come up. For example, searching for 'WNA' matches 'R#text(fill: blue)[WGGA]PG'. By default it will show the 25 best matching peptidoforms, but this number can be changed. The search can be restricted with a minimal score for the searched peptidoforms, which also makes the search faster, or with a minimal score for the alignment. Both these minimal score have to be in range 0.0 to 1.0.

Once the search is complete all matching peptidoforms (up to the maximum) will be shown below. The index indicates from which identified peptidoforms file the peptidoform originated as well as the index in that file. Clicking on the index selects this peptidoform in the peptidoform details pane. The sequence column shows the sequence of the peptidoform, with in blue the section that matched the search term. The match score (normalised mass based alignment score between 0 and 1) as well as the peptidoform score is shown in the last two columns.