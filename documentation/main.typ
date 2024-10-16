#set document(
    title: [Annotator documentation],
    author: ("Douwe Schulte"),
    date: auto,
)
#set par(
    justify: true,
)

#import "functions.typ": aside, button

= Introduction

A simple tool to help you manually discover the depths of your spectra one spectrum at a time. It can load MGF/mzML and thermo raw files. Once loaded you can select a scan and add you annotation while tweaking the exact settings for generating the annotation. The annotation itself is interactive to help you discover what the spectrum means. Which you can then export as nice images for use in other environments. 


#outline()

= Formatting guide

Buttons present in the annotator are displayed as #button[Button].

#aside[Any sidenote will be show like this.]

#set heading(numbering: "1.1")

#pagebreak()
#include "installing.typ"
#pagebreak()
#include "rawfiles.typ"
#pagebreak()
#include "identified_peptide_files.typ"
#pagebreak()
#include "annotate.typ"

