#import "functions.typ": aside, button, peptide

= Annotate

The annotate sections allow control over the annotation of theoretical spectra on the selected spectrum. The following sections are present:
+ The tolerance for matching theoretical peaks to experimental peaks can controlled and set to ppm or Thompson (mz). 
+ The maximum charge for precursors in the theoretical spectrum can be set. If this is not set it takes the charge from the spectrum. If this too is not set it takes +1.
+ The noise filter can be controlled, the noise floor is automatically determined and the noise filter disregards any peak below the factor times the noise floor. Setting this to 0.0 fully removes noise filtering.
+ The match mode indicates the method of determining the mz for theoretical peaks. Set to mono isotopic, average weight, or most abundant isotope.
+ An mz range for the theoretical peaks can be imposed. Setting only one side sets no bounds on the other side.
+ The model can be set to any predefined model. `All` allows all known fragmentation schemes. `None` only allows generation of the precursors. Custom models can also be created see @custom-model.
+ The peptide sequence contains the sequence for the peptide in ProForma notation, see @proforma.

Hitting #button[Annotate] generates the theoretical fragmentation for the given peptide with the given settings. The annotated spectrum will be shown below.

== ProForma <proforma>

The Annotator uses the #link("https://github.com/HUPO-PSI/ProForma")[ProForma 2.0] specification to specify the sequence. Here are some examples of valid sequences:

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
+ A DSSO cross-link between two lysines on two peptides (note the use of `//` versus `+` to indicate cross-linked peptides) \ #peptide("VAEINK[X:DSSO#XL1]SNGGTT//WAK[#XL1]INK")
+ A hydrolysed DSSO cross-linker \ #peptide("VAEINK[X:DSSO#XL1]SNGGTT")
+ An antibody Fab, encoding the disulfide bridge as `L-cystine (cross link)` from PSI-MOD \ #peptide("EVQLVESGGGLVQPGGSLRLSC[L-cystine (cross link)#XL1]AASGFNIKDTYIHWVRQAPGKGL\u{200B}EWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC[#XL1]SRWGGDGFYAMDYW​G\u{200B}QGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGC[L-cystine (cross link)#XL2]LVKDYFPEP\u{200B}VTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYIC[#XL2]NVNHKPSNTKVDKKVEPKS\u{200B}C[L-cystine (cross link)#XL3]DKT//DIQMTQSPSSLSASVGDRVTIT\u{200B}C[L-cystine (cross link)#XL4]\u{200B}RASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQ\u{200B}PEDFATYYC[#XL4]QQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVV\u{200B}C[L-cystine (cross link)#XL5]LLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHK\u{200B}VYAC[#XL5]EVTHQ\u{200B}GLSSPVTKSFNRGEC[#XL3]")

== Custom modifications <custom-modifications>

Custom modifications can be defined by opening the custom modifications section and hitting #button[Create new]. Add the chemical formula, if only the monoisotopic mass is known that can also be defined. The numerical ID is preset but the name must be set. For a custom modification metadata can be set, with description, synonyms, and identifiers for other identification systems (cross IDs).

Modifications can be defined as 'Modification' or 'Cross-linker'. For the first category placement rules can be defined that list the positions where this modification can be placed together with place-specific neutral losses and diagnostic ions.

For cross-linkers the length of the cross-linker can be defined as additional metadata item. Cross-linker placement rules can be defined as a symmetric cross-linker, which binds two identical positions on both sides of the cross-linker, or asymmetrical, where both sides bind different sides. For custom cross-linkers MS cleavage patterns can be defined. These are defined as two molecular formulas separated by `:`. If the theoretical fragmentation model allows MS cleavable cross-linkers these are allowed as breakage patterns. Additionally diagnostic ions can be added.

The custom modifications are stored in a separate JSON file on your computer, the path will be shown in the custom modification section. Updating the annotator will not remove any previously defined modifications. Additionally copying this file to another computer will copy the whole database to that computer. This can be used to move all your definitions to a new computer when upgrading or to aid a colleague with your definitions.

== Custom model <custom-model>

Custom models can be defined by opening the custom model section and hitting #button[Duplicate] on any existing model. For all main ion series (`abcxyz`) control is given over the section of the peptide that produces these fragments. 

Satellite ions series ('dvw') are formed at all location where the parent fragment ion are formed when turned on. The set of amino acid that gives rise to satellite ions can be controlled, if this field is left empty all amino acids will give rise to satellite ions. Additionally the maximal distance from the parent ion cleavage can be controlled. This is the number of side chains between the parent cleavage and the side chain that fragments. For example on a given peptide #peptide("HKSLG") the z3 fragment would contain #peptide("SLG"). The normal satellite ion w3 would be the side chain of Serine falling off (with a loss of 16 Da). The non-standard satellite ion 1w3 contains the same amino acids #peptide("SLG") but experiences a loss of the Leucine side chain in reference to the z3 fragment, with a mass of 43 Da. 1w3 means 1 side chain between the z3 cleavage and the w fragment site.

Neutral losses or gains can be #button[Select]ed from a set of common losses/gains and custom losses/gains can be specified. If there are losses that only occur when a certain amino acid is present in the fragment this can also be specified. For example COOH loss is seen from Asparagine in ETD. Additionally, side chain losses can sometimes be seen, primarily in ETD, where any side chain in the fragment can be lost. The maximal number of side chains lost can be controlled, commonly this should be kept at 1 unless there is good evidence for losing multiple side chains from the same fragment. As well as the amino acids that give rise to side chain losses can be controlled, if this is not specified all amino acids are assumed to be able to loose their side chains.

The charge range for any series can be controlled. Both sides of the charge range can be set to an absolute value or a value relative to the precursor. A charge range of 1 to the precursor charge is the default for most ions.

For all peptide fragmentation series ('abcdvwxyz') the existence of variant ions can be controlled. Variant ions are fragments that do not follow the common pathway for that ion but end up with some hydrogen difference. For all ions all options between -2 hydrogens and + 2 hydrogens can be turned on. Variant ions are depicted using a single quote #peptide(sym.quote.single) for each hydrogen lost or a middle dot #peptide(sym.dot.c) for each hydrogen gained. For example #peptide([z#sym.dot.c]) indicates a z ion with one hydrogen gained and #peptide([c#sym.quote.single#sym.quote.single]) a c ion with two hydrogens lost.

For precursors the neutral losses/gain and charge range are available.

For glycans the neutral losses/gains can be set for all glycan fragment types. The expected losses from diagnostic glycan fragments of one monosaccharide can be set. For glycans that reside on peptide fragments the glycan can be fragmented further depending on the fragmentation technique. This can be controlled with custom rules. Each rules applies to a set of amino acids, for example Asparagine for N-glycans, and potentially a set of fragment kinds ('abcdvwxyz'). Each rule defines if contained glycans do not undergo fragmentation ('full') or if they do ('core'), and in the latter case a range can be given. This range is the minimal and maximal depth in the glycan structure for structural glycans, so 0-1 would indicate an absence of the full glycan or the presence of the first monosaccharide. For compositional glycans this range indicates the number of monosaccharides expected, so 0-1 would indicate complete loss of the glycan or the inclusion of one the possible monosaccharides, all options will be generated. The charge range can be set separately for Y and B ions. Additionally, the fragmentation of structures from GNO can be turned on/off as well as the generation of fragmentation from compositional glycans (e.g. Hex1HexNac2) can be controlled.

Some modifications generate specific diagnostic ions, these can be turned on or off. Additionally, some modifications generate specific neutral losses, which can also be turned on or off.

Immonium ions are internal fragments that break on both sides of a single amino acid. These can be turned on/off. For immonium ions there are a lot of common neutral losses and gains, these can be controlled separately for each amino acid. 

If there are cross-links in the peptidoform it can be controlled if these are allowed to cleave in theoretical fragmentation. This only work for cross-link modifications that have defined cleavage rules, see @custom-modifications.

#aside[There are some reference sheets available and the end of this manual to help keep an overview of all fragmentation chemistry. #link("https://github.com/douweschulte/reference-sheets")[Or see douweschulte/reference-sheets on GitHub.]]