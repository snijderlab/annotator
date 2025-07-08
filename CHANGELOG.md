# Version 1.2.0 - Neutral losses & bug hunting

- Fixed manual not opening & updated manual
- Fixed links to open in the default browser
- Fixed some identified peptidoform issues (amongst other handle MSFragger files with missing columns)
- Fixed custom modification neutral losses not being shown (they were always applied)
- Fixed EAD glycan fragmentation
- Fixed custom glycan peptide fragments not being saved
- Improved peptide search speed, first search still takes a while but every subsequent search will be a lot quicker
- Improved handling of error in custom modifications and models files (errors are shown at the top of the list)
- Added support for PEAKS v13 Dia _de novo_
- Added support for factors for neutral losses (-2H<sub>2</sub>O instead of -H<sub>4</sub>O<sub>2</sub>)
- Added support for neutral losses from cross-linkers
 
Because of these last two the custom modifications and models JSON files have changed. The old files will still be read fine, but once edited the files will change to the updated JSON format. This is the reason why this update is listed as a minor release.

# Version 1.1.3

 - Fixed errors in parsing FragPipe files + added support for MSFragger
 - Updated mzdata containing fixes for Bruker data
 - Fixed some more smaller bugs

# Version 1.1.2

 - Fixed glycan related bugs
   - Naming of Y ions was inverted, meaning Y6 instead of Y1
   - The composition listed for Y composition ions was inverted, listing `Hex1` when 1 hexose was lost instead of when 1 hexose was retained

# Version 1.1.1

 - Fixed USI related bugs
 - Fixed theme related bugs

# Version 1.1.0 - Glycan visualisation & more complex fragmentation

 - Added glycan visualisations both full structures and fragments, in the spectrum, on peptides, and when searching / displaying modifications
 - Custom models can be saved (and backed up and shared between machines)
 - Added support for a basic CSV file to open as identified peptide (needing 'scan_index', 'z', 'raw_file', 'sequence', and possibly 'mode' as columns)
 - Added support for more complex fragmentation model
    - Added amino acid specific neutral losses to peptide fragments
    - Added side chain losses from peptide fragments
    - Added variant ions, where ions lose or gain hydrogen, z· means z + H, z' means z - H
    - Added control over glycan fragments on peptide fragments
    - Added control to neutral losses of glycan diagnostic ions and immonium ions
    - Added satellite ions distant from the parent breakage, <sup>n</sup>w indicates a w ion where the nth side chain (starting at 0) from the parent z ion fragments off
 - MS2 annotations are ordered on crude likelihood
 - General bug fixes, amongst others fixed high memory usage on merging spectra
