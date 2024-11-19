#import "functions.typ": aside, button

= Identified peptides files

There are quite some programs that export peptides with metadata based on spectral data. The annotator supports a selection of database matching and _de novo_ software file formats. To load these files use the button #button[Load identified peptides file] or drag in a file from the file explorer. Any number of files can be loaded at the same time. See the table below for supported formats. Opened identified peptides files can be closed with #button[Close file] in the peptides details pane.

#table(
  columns: (auto, auto),
  table.header([*Open format*], [*Comment*]),
  [Fasta],
  [No header requirements],
  [mzTab],
  [v1.0 & Casanovo],
)

#table(
  columns: (auto, auto),
  table.header([*Software*], [*Versions/Formats*]),
  [DeepNovoFamily],
  [DeepNovo, PointNovo, BiatNovo, PGPointNovo],
  [InstaNovo],
  [1.0.0],
  [MaxQuant],
  [msms, msms scans, novo msms scans, & silac],
  [MSFragger],
  [21 and 22],
  [Novor],
  [Denovo and PSM],
  [OPair],
  [common version],
  [Peaks],
  [X, X+, 11, 12, Ab, DB peptides, DB PSM, & DB protein-peptide],
  [PowerNovo],
  [1.0.1],
  [pLink],
  [2.3],
  [Sage],
  [0.14],
)

#aside[If you have data from an unsupported version of a supported program please open an issue on GitHub and give an example file so that the support can be extended.]

#aside[If you have data from an unsupported program that you think should be supported open an issue on GitHub to discuss.]

== Peptide details

Once a file is opened the peptide details pane opens. Select the right identified peptides file from the drop down menu and select the right peptide by peptide index. Peptides can be browsed efficiently by selecting the peptide index box and using the arrow up/down keys.

Once a peptide is selected it will show common metadata in a structured format, followed by an overview of the peptide, with the local confidence (if present in the file) depicted by blue squares, and modifications depicted by blue dots. Hover over a modified amino acid to see the full definition. Terminal modifications will be depicted as blue dots on their respective terminal.

Below the structured metadata follows in table format all other metadata from the file format.

Use #button[Load] to select the right spectrum (if the corresponding raw file is open) and load the details for annotation. This will load the sequence, charge, and method if these are available.

== Peptide search 

Once at least one identified peptides file is open all peptides can be searched for sequence patterns. Type the search pattern as a ProForma (see @proforma) peptide in the search box and hit #button[Search] to search. The search is based on #link("https://doi.org/10.1021/acs.jproteome.4c00188")[mass based alignment] so any peptide matching the mass pattern of the search will come up. For example, searching for 'WNA' matches 'R#text(fill: blue)[WGGA]PG'. It will show the 25 best matching peptides. The search can be restricted with a minimal score for the searched peptides, which also makes the search faster, or with a minimal score for the alignment. Both these minimal score have to be in range 0.0 to 1.0.

Once the search is complete all matching peptides (up to 25) will be shown below. The index indicates from which identified peptides file the peptide originated as well as the index in that file. Clicking on the index selects this peptide in the peptide details pane. The sequence column shows the sequence of the peptide, with in blue the section that matched the search term. The match score (normalised mass based alignment score between 0 and 1) as well as the peptide score is shown in the last two columns.