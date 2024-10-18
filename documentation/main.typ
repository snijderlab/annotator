#import "functions.typ": crate, version

#set document(
    title: [Annotator documentation - #version],
    author: crate.package.authors,
    date: auto,
)
#set par(
    justify: true,
)
#set page(
  header: align(
    right + horizon,
    [#box(image("../src-tauri/icons/annotator_icon_export_grey.svg", height: 0.7em))nnotator documentation - #version - #link(crate.package.repository)[Open source at GitHub]]
  ),
  numbering: "1"
)

#import "functions.typ": aside, button

#grid(
    columns: (auto, auto),
    image("../src-tauri/icons/annotator_icon_export.svg"),
    [= Annotator

    A simple tool to help you manually discover the depths of your spectra, one spectrum at a time. Once your rawfiles are loaded you can select a scan and add you annotation with full control over theoretical fragment genration. The interactive spectrum will help you discover what the spectrum means, with full control to export gorgeous images.

    GitHub: #link(crate.package.repository)[#crate.package.repository] \
    License: #crate.package.license \
    Version: #version 
    ],
)

#outline(indent: auto)

= Formatting guide

Buttons present in the annotator are displayed as #button[Button].

#aside[Any sidenote will be show like this.]

#set heading(numbering: "1.1")

#pagebreak()
#include "installing.typ"
#include "rawfiles.typ"
#include "identified_peptide_files.typ"
#include "annotate.typ"
#include "spectrum.typ"
#include "tools.typ"

