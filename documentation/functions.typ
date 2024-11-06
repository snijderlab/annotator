#let lightgrey = 192;
#let grey = 128;
#let darkgrey = 32;

#let aside(text) = {
  block(
    fill: color.linear-rgb(grey, grey, grey, 64),
    stroke: 2pt + color.linear-rgb(grey, grey, grey, 255),
    inset: 8pt,
    width: 100%,
    text,
  )
}

#let button(text) = {
  h(4pt)
  box(
    fill: color.linear-rgb(71, 153, 217, 255),
    inset: 0pt,
    radius: 8pt,
    outset: 4pt,
    text
  )
  h(4pt)
}

#let key(it) = {
  h(3pt)
  box(
    fill: color.linear-rgb(lightgrey, lightgrey, lightgrey, 255),
    inset: 0pt,
    radius: 2pt,
    outset: 3pt,
    [#text(font: "Roboto Mono", size: 0.9em, fill: color.linear-rgb(darkgrey, darkgrey, darkgrey, 255), it)]
  )
  h(3pt)
}

#let peptide(it) = {
  h(3pt)
  box(
    fill: color.linear-rgb(lightgrey, lightgrey, lightgrey, 255),
    inset: 0pt,
    radius: 0pt,
    outset: 3pt,
    text(font: "Roboto Mono", size: 0.9em)[#it]
  )
  h(3pt)
}

#let crate = toml("../backend/Cargo.toml")
#let version = [v#crate.package.version]