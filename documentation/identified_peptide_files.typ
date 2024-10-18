#import "functions.typ": aside, button

= Identified peptides files

There are quite some programs that export peptides with metadata based on spectral data. Both database matching and de novo software. These files can be loaded using the button #button[Load identified peptides file] or by dragging in a file from the file explorer. Any number of files can be loaded at the same time. See the table below for supported formats. Opened identified peptides files can be closed with #button[Close file] in the peptides details pane.

#table(
  columns: (auto, auto),
  table.header([*Software*], [*Versions/Formats*]),
  [Fasta],
  [common version],
  [MaxQuant],
  [msms, msms scans, novo msms scans, & silac],
  [MSFragger],
  [21 and 22],
  [Novor],
  [Denovo and PSM],
  [OPair],
  [common version],
  [Peaks],
  [X, X+, 11, 12, & Ab],
  [Sage],
  [0.14],
)

#aside[If you have data from an unsupported version of a supported program please open an issue on GitHub and give an example file so that the support can be extended.]

#aside[If you have data from an unsupported program that you think should be supported open an issue on GitHub to dicuss.]

== Peptide details

Once a file is opened the peptide details pane opens. Select the right identified peptides file from the dropdown menu and select the right peptide by peptide index. Peptides can be browsed efficiently by selecting the peptide index box and using the arrow up/down keys.

Once a peptide is selected it will show common metadata in a structured format. Followed by an overview of the peptide, with the local confidence (if present in the file) depicted by blue squares, and modifications written on their respective amino acids. Terminal modifications will be depicted as modifications on the symbol '⚬' on the respective terminal.

Below the structured metadata follows in table format all other metadata from the file format.

Use #button[Load] to select the right spectrum (if the corresponding raw file is open) and load the details for annotationd. This will load the sequence, charge, and method if these are available.

== Peptide search 

Once at least on identified peptides file is open all peptides can be searched for sequence patterns. Type the search pattern as a ProForma (see @proforma) peptide in the search box and hit #button[Search] to search. The search is based on #link("https://doi.org/10.1021/acs.jproteome.4c00188")[mass based alignment] so any peptide matching the mass pattern of the search will come up. For example searching for 'WNA' matches 'R#text(fill: blue)[WGGA]PG'. It will show the 25 best matching peptides. The search can be restricted with a minimal score for the peptides that are searched in, which also makes the search faster. Additionally there can be a minimal score for the alignment specified that disregards peptides that score too low. Both these minimal score have to be in range 0.0 to 1.0.

Once the search is complete all matching peptides (up to 25) will be shown below. The index indicates from wich identified peptides file the peptide originated as well as the index in that file. Clicking on the index selects this peptide in the peptide details pane. The sequence column shows the sequence of the peptide, with in blue the section that matched the search term. The match score (normalised mass based alignment score between 0 and 1) as well as the peptide score is shown in the last two columns.