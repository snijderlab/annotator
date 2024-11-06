= Tools

== Search modification

This can be used to find modifications in any of the ontologies (modification databases) or find more information on a single modification. The search box can be used as follows:

+ When a mass is given any modification that fits within the tolerance is returned. For example searching for `16` gives the following results: #table(columns: (auto, auto, auto, auto), [Name],	[Id],	[monoisotopic mass],	[Formula],
[U:Oxidation],	[UNIMOD:35], [15.995 Da / 15.999 Da / 15.995 Da],	[O1],
[U:Methyl:2H(2)],	[UNIMOD:284],	[16.028 Da / 16.039 Da / 16.028 Da],	[C12H2],
[U:Carboxy->Thiocarboxy],	[UNIMOD:420],	[15.977 Da / 16.068 Da / 15.977 Da],	[O-1S1],
[U:Ala->Ser],	[UNIMOD:540],	[15.995 Da / 15.999 Da / 15.995 Da],	[O1])
+ When a formula is given any modification that has the same molecular formula is returned. For example searching for `Formula:O` gives the following results: #table(columns: (auto, auto), [Name], [Id], [U:Oxidation],	[UNIMOD:35],
[U:Ala->Ser],	[UNIMOD:540],
[U:Phe->Tyr],	[UNIMOD:569],
[M:(2S,3R)-3-hydroxyasparagine], [MOD:35],)
+ When a glycan composition is given any glycan with the same composition is returned (all isomeric information is ignored), which makes it possible to find topologies for given glycan compositions. Another way of finding glycan compositions based on composition is using the #link("https://glycosmos.org/glycans/gnome")[GlyCosmos GNOme structure browser].
+ When a modification is given the details for that modification are displayed.

== Isotopic distribution

Here the isotopic distribution for any molecular formula can be generated. All weights of the formula will be displayed as well as a graph with all isotopic peaks. Hovering over the peaks gives additional details on that specific peak. For the generation of the isotopic distribution an averagine model is used that slightly overestimates the prevalence of higher weight isotopes, especially for elements with multiple isotopes. Any element with a defined monoisotopic weight, so any element that is stable and naturally occurring, can be used in the formula.

== GNOme structure browser

The GNOme structure browser is a tool developed by GlyCosmos to find GNOme structures based on glycan composition. This tools is embeded into the Annotator for ease of use. Find more information on this tool here: #link("https://glycosmos.org/glycans/gnome")[https://glycosmos.org/glycans/gnome].