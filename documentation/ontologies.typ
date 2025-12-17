#import "functions.typ": aside, button, peptidoform

= Ontologies

In the last section the current state of the ontologies can be seen. For each ontology the version and last updated date are displayed (if set properly by the ontology). The SHA256 of the loaded file is shown alongside the number of contained modifications. Most ontologies can be updated by clicking on the #button("Internet"), this automatically downloads the ontologies from the canonical locations, parses the files and updates the data. Any error will be shown below the table. For more advanced use cases or when the ontologies have moved (or for RESID which cannot be automatically downloaded) a local file can be specified to update the ontologies with. Note that it depends on the ontology how many files and of which file type have to be supplied. Any of these files can also be supplied in gzipped form when supplied with the `.<ext>.gz`. If the update was successful the ontology is directly updated and can immediately be used and will be in the exact same state when the annotator is closed and opened again. Note that during the updating itself no other calculations can be done with the Annotator and that depending on the internet speed it might take a bit, especially for GNOme as that contains quite a big collection of glycans.

#table(
  columns: (auto, auto),
  [Ontology], [File(s)],
  [Unimod], [#link("https://unimod.org/xml/unimod.xml")[unimod.xml] (not unimod_tables.xml)],
  [PSI-MOD], [#link("https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/refs/heads/master/PSI-MOD.obo")[PSI-MOD.obo]],
  [GNOme], [#link("https://purl.obolibrary.org/obo/gno.obo")[gno.obo] and #link("https://glycosmos.org/download/glycosmos_glycans_list.csv")[glycosmos_glycans_list.csv]],
  [XLMOD], [#link("https://raw.githubusercontent.com/HUPO-PSI/xlmod-CV/refs/heads/main/XLMOD.obo")[XLMOD.obo]],
  [RESID], [#link("ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/RESIDUES.XML")[RESIDUES.XML]]
)