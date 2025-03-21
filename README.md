A simple tool to help you manually discover the depths of your spectra one spectrum at a time. It can load MGF/mzML and thermo raw (see note at the end) files (only in centroid mode, also deconvolution must be done separately if needed). Once loaded you can select a scan and add you annotation while tweaking the exact settings for generating the annotation. The annotation itself is interactive to help you discover what the spectrum means. Which you can then export as nice images for use in other environments. 

## Peptide sequence

It uses the [ProForma](https://github.com/HUPO-PSI/ProForma) specification to specify the sequence. Here are some examples of valid sequences:

* `VAEINPSNGGTTFNEKFKGGKATJ` Normal aminoacids
* `EM[L-methionine sulfoxide]EVEES[UNIMOD:21]PEK` Modifications using [unimod](http://www.unimod.org) and [PSI-MOD](https://www.ebi.ac.uk/ols/ontologies/mod)
* `TFNEKF[+15.9949]KGGKATJ` Modifications using raw masses
* `TFNEKF[Formula:O]KGGKATJ` Modifications using elemental formula
* `TFNEKF[Glycan:HexNAc1Hex2]KGGKATJ` Modifications glycan compositions
* `[+16]-TFNEKFKGGKATJ-[Methyl]` Terminal modifications
* `<15N>TFNEKFKGGKATJ` Global isotope modifications (all Nitrogen is 15N)
* `<[S-carboxamidomethyl-L-cysteine]@C>AVYYCSRWGGDGFYAMDYWGQG` Global modifications (all C are carboxamidomethylated)
* `[UNIMOD:374]?TFNEKFCKGGCKATJ` Modifications where the location is unknown
* `TFNEKFC[UNIMOD:374#g1]KGGC[#g1]KATJ` (identical to the one above)
* `TFNEKF(CKGGCK)[UNIMOD:374#g1]ATJ` (identical to the one above)
* `VAEINPSNGGTT+FNEKFKGGKATJ` Multimeric spectra, meaning two separate peptides are in your spectrum at the same time
* `VAEINPSNGGTT/2[1Na+,1H+]` Defined charge and adduct ions
* `VAEINK[X:DSSO#XL1]SNGGTT//WAK[#XL1]INK` A DSSO cross-link between two lysines on two peptides (note the use of `//` versus `+` to indicate cross-linked peptides)
* `VAEINK[X:DSSO#XL1]SNGGTT` A hydrolysed DSSO cross-linker

## Note on custom modifications

Custom modifications can be defined, these allow diagnostic ions and neutral losses to be defined. Additionally custom cross-linkers can be defined to have certain cleavage patterns that can then be searched for in the annotator. The custom modifications are stored in a separate json file on your computer. Updating the annotator will not remove any previously defined modifications. Additionally copying this file to another computer will copy the whole database to that computer. This can be used to move all your definitions to a new computer when upgrading or to aid a colleague with your definitions.

## Installing

### Using winget

On windows use `winget install --id Snijderlab.Annotator`.

### From binary 

See [releases](https://github.com/snijderlab/annotator/releases) for the latest release, here you will also find the prebuilt binaries for your architecture.

### From source

To build from source. Clone the repository. And build with cargo. Make sure you have installed [Rust](https://www.rust-lang.org/tools/install) and [Tauri](https://tauri.app/) beforehand.

## Thermo RAW files

The .NET 8.0 runtime is needed to open Thermo RAW files. [Which can be downloaded here.](https://dotnet.microsoft.com/en-us/download/dotnet/8.0) Additionally on windows you can use `winget install Microsoft.DotNet.Runtime.8` for a quick install.