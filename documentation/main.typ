#import "functions.typ": crate, version, key, button, aside

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
    [#box(image("../backend/icons/annotator_icon_export_grey.svg", height: 0.7em))nnotator documentation - #version - #link(crate.package.repository)[Open source at GitHub]]
  ),
  numbering: "1"
)
#set text(1em, weight: "regular", font: "Roboto", hyphenate: true)
#show heading: it => [
  #set text(1em, weight: "regular", font: "Roboto Slab")
  #block(it.body)
]
#set heading(outlined: false)

#grid(
    columns: (auto, auto),
    image("../backend/icons/annotator_icon_export.svg"),
    [= Annotator

    #crate.package.description

    GitHub: #link(crate.package.repository)[#crate.package.repository] \
    License: #crate.package.license \
    Version: #version 
    ],
)

#outline(indent: auto, depth: 2)

#pagebreak()

= Formatting guide

Buttons are displayed as #button[Button]. Keys used in key combinations as #key[Key].

#aside[Any side note will be shown like this.]

#set heading(numbering: "1.1", outlined: true)


#include "installing.typ"
#include "rawfiles.typ"
#include "identified_peptide_files.typ"
#include "annotate.typ"
#include "spectrum.typ"
#include "tools.typ"
#include "exporting.typ"

#{
  show heading: none
  heading()[Amino acids cheat sheet]
  page(header: none, margin: 0pt, image(width: 100%, height: 100%, "amino_acids_cheat_sheet_1.2_print.svg"))
  heading()[Fragmentation cheat sheet]
  page(header: none, margin: 0pt, image(width: 100%, height: 100%, "fragmentation_cheat_sheet_1.1_print.svg"))
  heading()[ProForma cheat sheet]
  page(header: none, margin: 0pt, image(width: 100%, height: 100%, "proforma_cheat_sheet_1.0_print.svg"))
}