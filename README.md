A simple tool to help you manually discover the depths of your spectra one spectrum at a time. It can load MGF files (only in centroid mode, also do not forget to deconvolute if you have TD data). Once loaded you can select a scan and add you annotation while tweaking the exact settings for generating the annotation. The annotation itself is interactive to help you discover what the spectrum means. Which you can then export as nice images for use in other environments. 

## Peptide sequence

It uses the [ProForma](https://github.com/HUPO-PSI/ProForma) specification to specify the sequence, it does not handle every last detail of this specification yet, for details see [rustyms](https://github.com/snijderlab/rustyms). Here are some examples of valid sequences:

* `VAEINPSNGGTTFNEKFKGGKATJ` Normal aminoacids
* `EM[L-methionine sulfoxide]EVEES[UNIMOD:21]PEK` Modifications using [unimod](http://www.unimod.org) and [PSI-MOD](https://www.ebi.ac.uk/ols/ontologies/mod)
* `TFNEKF[+15.9949]KGGKATJ` Modifications using raw masses
* `TFNEKF[Formula:O]KGGKATJ` Modifications using elemental formula
* `TFNEKF[Glycan:HexNAc1Hex2]KGGKATJ` Modifications glycan compositions
* `[+16]-TFNEKFKGGKATJ-[Methyl]` Terminal modifications
* `<15N>TFNEKFKGGKATJ` Global isotope modifications (all Nitrogen is 15N)
* `<[S-carboxamidomethyl-L-cysteine]@C>AVYYCSRWGGDGFYAMDYWGQG` Global modifications (all C are carboxamidomethylated)

## Installing

See [releases](https://github.com/snijderlab/annotator/releases) for the latest release, here you will also find the prebuilt binaries for your architecture.