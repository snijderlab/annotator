#import "functions.typ": crate

= Installing

== Using winget

On windows use `winget install --id Snijderlab.Annotator`.

== From binary 

See #link(crate.package.repository + "/releases")[GitHub releases] for prebuilt binaries for your architecture as well as information on the releases.

== From source

To build from source, clone the repository, and build with cargo. Make sure you have installed #link("https://www.rust-lang.org/tools/install")[Rust] and #link("https://tauri.app/")[Tauri] beforehand.

#raw("git clone " + crate.package.repository + " " + crate.package.name + "\ncd " + crate.package.name + "\ncargo tauri dev", lang: "sh")