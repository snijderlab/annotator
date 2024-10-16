#let aside(text) = {
  block(
    fill: color.linear-rgb(128, 128, 128, 64),
    stroke: 2pt + color.linear-rgb(128, 128, 128, 255),
    inset: 8pt,
    width: 100%,
    text,
  )
}

#let button(text) = {
  box(
    fill: color.linear-rgb(71, 153, 217, 255),
    inset: 0pt,
    radius: 8pt,
    outset: 4pt,
    text,
  )
}