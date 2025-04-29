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
