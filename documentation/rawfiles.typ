#import "functions.typ": aside, button, peptide

= Rawfiles

The annotator supports mgf, mzML, IndexedMzML, Bruker .tdf, and Thermo raw files (needs setup see @thermo-raw). Rawfiles can be loaded using the button #button[Load raw data file] or by dragging in a file from the file explorer. Any number of files can be loaded at the same time. Once opened the file will be added to the list of opened files as 'RX:filename' with two input fields to select spectra and a #button[Close file] button to close the file. Spectra can be selected on index, 0 based index of the spectrum in the file, or on native id. Native id is the vendor specific textual identifier for the spectrum, #link("https://docs.rs/mzdata/latest/mzdata/meta/enum.NativeSpectrumIdentifierFormatTerm.html#variants")[for more detailed info see the documentation in mzdata].

Files in profile mode can be loaded and automatically peak picked in the annotator. For top and middle down data it is recommended to deconvolute the data before loading it in the annotator.

== Other formats

Other raw file formats can be converted using available converters. For convenience some common converters are listed below.

#table(
  columns: (auto, auto),
  table.header([*Converter*], [*Supported formats*]),
  link("https://proteowizard.sourceforge.io/download.html")[Proteowizard MSConvert],
  link("https://proteowizard.sourceforge.io/doc_users.html#SupportedFormats")[AB/Sciex T2D, Agilent MassHunter, Bruker BAF/Data Exchange/FID/TDF/U2/YEP, Mascot Generic, MS1, MS2, MZ5, mzML, mzXML, Sciex WIFF/WIFF2, Shimadzu LCD, Thermo raw, UIMF, Waters raw],
  link("https://openms.readthedocs.io/en/latest/about/installation.html")[TOPP FileConverter],
  link("https://openms.readthedocs.io/en/latest/getting-started/types-of-topp-tools/file-handling.html#converting-your-files-to-mzml")[mzData, mzXML, ANDI/MS],
  link("https://openms.readthedocs.io/en/latest/about/installation.html")[TOPP DTAExtractor],
  link("https://openms.readthedocs.io/en/latest/getting-started/types-of-topp-tools/file-handling.html#converting-between-dta-and-mzml")[Sequest DTA],
)

== Selected spectra

Once a spectrum is selected it will be listed below the respective raw file. It will show its index and native id as well as a #button[Unselect] button to undo the selection.

=== Multiple spectra

If multiple spectra are selected at the same time these spectra will be merged before being annotated.

== Thermo raw files <thermo-raw>

The .NET 8.0 runtime is needed to open Thermo raw files, #link("https://dotnet.microsoft.com/en-us/download/dotnet/8.0")[which can be downloaded here.] Additionally on windows you can use `winget install Microsoft.DotNet.Runtime.8` for a quick install. Once this is installed Thermo raw files can be loaded as any other file.

== Clipboard

Some programs allow copying a spectrum into the clipboard, use the #button[Load Clipboard] button to load such a spectrum from the clipboard. Currently spectra from selected Bruker, Stitch, Sciex, and Thermo programs are supported. 

#aside[If you find another program that allows this behaviour please open an issue on GitHub for the Annotator, with the name and version of the program in question along with an example of the format.]

== Universal Spectrum Identifier (USI)

A #link("https://www.psidev.info/usi")[USI] is a text based format that identifies any spectrum in any publicly accessible dataset. It is build using the following parts.

#table(columns: (auto, auto, auto),
table.header([*Part*], [*Example*], [*Comment*]),
[Prefix], [mzspec], [required],
[Collection], [PXD043489], [required],
[Dataset], [20201103_F1_UM5_Peng0013\_-?SA_139H2_InS_Elastase.raw], [required],
[Index type], [scan], [required, can be scan/nativeId/index],
[Index], [11809], [required, content depens on Index type],
[Interpretation], [VSLFPPSSEQLTSNASVV], [optional, can be any ProForma sequence],
[Provenance], [PR-G47], [optional]
)

These parts have to be strung together using colons ':', resulting in the following string:

#h(12pt)#peptide[mzspec:-?PXD043489:-?20201103_F1_UM5_Peng0013_SA_139H2_InS_Elastase.raw:-?scan:-?11809:-?VSLFPPSSEQLTSNASVV]

Pasting a USI into the Annotator downloads that spectrum from the internet and loads the interpretation (if given) into the peptide field.