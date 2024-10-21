#import "functions.typ": aside, button, peptide

= Annotate

The annotate sections allow control over the annotation of theoretical spectra on the selected spectrum. The following sections are present:
+ The tolerance for matching theoretical peaks to experimental peaks can controlled and set to ppm or Thompson (mz). 
+ The maximum charge for precursors in the theoretical spectrum can be set. If this is not set it takes the charge from the spectrum. If this is not set it takes +1.
+ The noise filter can be controlled, the noise floor is automatically determined and the noise filter disregards any peak below the factor times the noise floor. Setting this to 0.0 fully removes noise filtering.
+ The match mode indicates the method of determining the mz for theoretical peaks. Set to mono isotopic, average weight, or most abundant isotope.
+ A mz range for the theoretical peaks can be imposed. Setting only one side sets no bounds on the other side.
+ The model can be set to any predefined model. `All` allows all known fragmentation schemes. `None` only allows generation of the precursors. `Custom` allows creating a custom model, see @custom-model.
+ The peptide sequence contains the sequence for the peptide in ProForma notation, see @proforma.

Hitting #button[Annotate] generates the theoretical fragmentation for the given peptide with the given settings. The annotated spectrum will be shown below.

== Custom model <custom-model>

In the custom model section many aspects of the theoretical fragment generation can be controlled. 

For all main ion series (`abcdvwxyz`) control is given over the section of the peptide that produces these fragments. Neutral losses or gains can be #button[Select]ed from a set of common losses/gains and additionally custom losses/gains can be specified. The charge range for any series can be controlled. Both sides of the charge range can be set to an absolute value or a value relative to the precursor. A charge range of 1 to the precursor charge is the default.

For precursors the neutral losses/gain and charge range are available.

For glycans the neutral losses/gains can be set for all glycan fragment types. The charge range can be set specifically for Y and Oxonium ions. Additionally the fragmentation of structures from GNO can be turned on/off as well as the generation of fragmentation from compositional glycans (e.g. Hex1HexNac2) can be controlled.

Some modifications generate specific diagnostic ions, these can be turned on or off. Additionally some modifications generate specific neutral losses, which can also be turned on or off.

Immonium ions are internal fragments that break on both sides of a single amino acid. These can be turned on/off. Turning these on also allows common immonium related ions, which are common neutral losses/gains often seen for immonium ions.

In some fragmentation techniques sidechains of amino acids can be lost as neutral loss. This can be turned on for precursors.

If there are cross-links in the peptidoform it can be controlled if these are allowed to cleave in theoretical fragmentation. This only work for cross-link modifications that have defined cleavage rules, see @custom-modifications.

#aside[There are some reference sheets available to help keep on overview of all these fragmentation chemistry. #link("https://github.com/douweschulte/reference-sheets")[See douweschulte/reference-sheets on GitHub.]]

== ProForma <proforma>

It uses the #link("https://github.com/HUPO-PSI/ProForma")[ProForma] specification to specify the sequence. Here are some examples of valid sequences:

+ Normal amino acids \ #peptide("VAEINPSNGGTTFNEKFKGGKATJ")
+ Modifications using #link("http://www.unimod.org")[UNIMOD], #link("https://www.ebi.ac.uk/ols/ontologies/mod")[PSI-MOD], RESID, XL-MOD, and GNO \ #peptide("EM[L-methionine sulfoxide]EVEES[UNIMOD:21]PEK")
+ Modifications using raw masses \ #peptide("TFNEKF[+15.9949]KGGKATJ") 
+ Modifications using elemental formula \ #peptide("TFNEKF[Formula:O]KGGKATJ")
+ Modifications glycan compositions \ #peptide("TFNEKF[Glycan:HexNAc1Hex2]KGGKATJ") 
+ Terminal modifications \ #peptide("[+16]-TFNEKFKGGKATJ-[Methyl]")
+ Global isotope modifications (all Nitrogen is 15N) \ #peptide("<15N>TFNEKFKGGKATJ")
+ Global modifications (all C are carboxamidomethylated) \ #peptide("<[S-carboxamidomethyl-L-cysteine]@C>AVYYCSRWGGDGFYAMDYWGQG")
+ Modifications where the location is unknown \ #peptide("[UNIMOD:374]?TFNEKFCKGGCKATJ")
+ Modification of unknown position specified on two positions \ #peptide("TFNEKFC[UNIMOD:374#g1]KGGC[#g1]KATJ")
+ Modification of unknown position specified on a subset of the peptide \ #peptide("TFNEKF(CKGGCK)[UNIMOD:374#g1]ATJ")
+ Chimeric spectra, meaning two separate peptides are in your spectrum at the same time \ #peptide("VAEINPSNGGTT+FNEKFKGGKATJ")
+ Defined charge, especially good for chimeric cases that have different charges \ #peptide("VAEINPSNGGTT/2")
+ Defined charge and adduct ions \ #peptide("VAEINPSNGGTT/2[1Na+,1H+]")
+ A DSSO cross-link between two lysines on two peptides (note the use of `//` versus `+` to indicate cross-linked peptides) \ #peptide("VAEINK[X:DSSO#XL1]SNGGTT//WAK[#XL1]INK")
+ A hydrolysed DSSO cross-linker \ #peptide("VAEINK[X:DSSO#XL1]SNGGTT")

== Custom modifications <custom-modifications>

Custom modifications can be defined, open the custom modifications section and hit #button[Create new]. Add the chemical formula, if only the monoisotopic mass is known that can also be defined. The numerical Id is preset but the name must be set. For a custom modification metadata can be set, with description, synonyms, and identifiers for other identification systems (cross IDs).

Modifications can be defined as 'Modification' or 'Cross-linker'. For the first category placement rules can be defined that list the positions where this modification can be placed together with place specific neutral losses and diagnostic ions.

For cross-linkers the length of the cross-linker can be defined as additional metadata item. Cross-linker placement rules can be defined as a symmetric cross-linker, which binds two identical positions on both sides of the cross-linker, or assymetrical, where both sides bind different sides. For custom cross-linkers MS cleavage patterns can be defined. These are defined as two molecular formulas separated by `:`. If the theoretical fragmentation model allows MS cleavable cross-linkers these are allowed as breakage patterns. Additionally diagnostic ions can be added.

The custom modifications are stored in a separate JSON file on your computer, the path will be shown in the custom modification section. Updating the annotator will not remove any previously defined modifications. Additionally copying this file to another computer will copy the whole database to that computer. This can be used to move all your definitions to a new computer when upgrading or to aid a colleague with your definitions.