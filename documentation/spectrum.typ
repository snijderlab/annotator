#import "functions.typ": key, button

= Interactive spectrum

== Legend

Once an annotation is done the annotated spectrum will be shown below the annotation settings. It starts with a legend of the colours. For any of these legend elements if hovered over it highlights all peaks that match. If clicked this highlighting will stay until clicked again. Multiple highlights can be stacked which will show any peaks matching any of the criteria.

== Peptidoform legend

Below the legend is an overview of the peptidoform(s) that are annotated. All main ion series are displayed as flags in the corner between the amino acids. Hovering highlights any fragment from that position similar to highlighting in the legend. Clicking makes the highlighting permanent in the same way. If there are multiple peptidoforms they are prefixed by the peptidoform index (peptidoform index followed by a dot followed by the peptidoform index). Hovering/clicking highlights all peptidoform fragments.

In the settings section the peptidoform can be set to use a more compact representation. Additionally, parts of the peptidoform (and any matching fragment) can be highlighted in a different colour.

== Spectrum

In the spectrum hovering over a peak shows the mz value as well as all annotations if there are multiple annotations for that peak. Zooming in can be done by clicking at one point to select one mz bound and move the mouse and release it at another point to set the other mz bound as well as the intensity bound. A selection window will appear to indicate the section shown after zooming in. Additionally zooming can also be done with the scroll wheel. Scrolling zooms in or out on the spot where the mouse is located. Scrolling with #key[Shift] pans the spectrum to the left or right on the mz axis. Scrolling with #key[Control] zooms in on the intensity axis. Lastly, using #key[Control] with #key[+] and #key[-] zooms in or out on the middle of the spectrum, and #key[Control] with #key[0] zooms to the original zoom level of the spectrum. The distance between peaks can be annotated by clicking on one peak, dragging to the next and releasing on another peak. Such labels can be removed by clicking on them, or all labels can be removed at once in the settings.

=== Peaks settings 

The theoretical spectrum can be shown below the x axis. Displaying the unassigned peaks can be toggled. The colouration can be set to the ion type (default), to the peptidoform, the peptidoform, or grey (none). Also there is an option to remove all distance labels.

=== Label settings 

By default the labels for peaks are shown for any peak within the 90% of intensity, meaning any peak having an intensity of at least $(1.0 - 0.9) * "max_intensity"$. This threshold can be controlled in the label settings section. Showing the mz value can be controlled in the same way for peaks. The 'Manually force show' allows to show either the label or mz value for any peak. Select the mode of interest and any peak clicked will show that information. The 'Hide' mode allows to hide the information (label and mz) for any peak. The button #button[Clear] clears all manually forced labels. The checkboxes 'Show in label' allow the hiding/showing of certain parts of the label to slim down the information shown. Lastly, the labels can be set to be 90 degrees rotated to make overlap less likely. Any of these changes will also show up in the exported image.

=== Spectrum settings

This section allows zooming to precise numbers. Additionally, the number of tickmarks for the x and y axis can be set. The peaks can be set to use square root intensity instead of linear intensity to see a bigger dynamic range. The y axis can be changed to show percentages instead of raw intensities as well.

== Error graph

The graph below the spectrum by default shows the ppm error for all theoretical peaks to their experimental matches. The points are at the same location on the x axis and any zooming in the spectrum will zoom in the error graph as well. In the settings section the points can be set to reflect the peak intensity with the point diameter and the range of y axis can be controlled.

Using the settings the y axis can be changed from relative (ppm) to absolute (Thompson). Additionally, the x axis can be set to reflect errors of unassigned peaks instead of the default assigned peaks. The errors of the unassigned peaks are calculated as being the error to the closest theoretical fragment for all unassigned peaks. In this mode at least one ion series needs to be picked to work as reference series. Commonly the unassigned mode is used together with absolute error mode.

Using the unassigned mode is most helpful to detect a series of peaks that is a constant offset from the reference peaks. This indicates that the assignment is off just before the constant offset series. See @fig-error-graph for an example where an annotation from the N terminus matched the reference series but has a constant offset from a certain point. This indicates that at that point the actual peaks are 16 Thompson off from the reference series likely indicating an oxidation. This way of using the error graph is most helpful when annotating _de novo_ top down or middle down sequences to determine where the annotation deviates from the experimental evidence.

#figure(image("error_graph.svg"), caption: [Example error graph in absolute error mode with the errors for unassigned peaks relative to some N terminal ion series. Visible is that four peaks match perfectly but right after that four additional peaks have a constant offset of 16 Thompson to the reference series.] ) <fig-error-graph>

== Statistics overview

The statistic overview gives detailed statistics about the match. Listed are:
+ The precursor mass for each peptidoform. If a peptidoform is part of a peptidoform no mass is listed. 
+ The number of fragments found out of all theoretical fragments.
+ The number of peaks annotated out of all peaks.
+ The total intensity annotated out of the sum of the intensity of all peaks.
+ The sequence positions covered by at least one matched fragment. If some positions do not generate theoretical fragments, for example with cross-linking loop links, this statistic will be split into one regarding the full peptidoform and one only taking the possible locations into account.
+ Peaks false match chance a false discovery metric at the peaks level. Estimating the FDR by permutation. The spectrum is shifted by -25 to +25 Da plus π (to have non integer offsets) and the number of matches is counted. The percentage is the number of actual matches divided by the average found for the shifted spectra. The number between brackets denotes the number of standard deviations the actual matches are from the shifted matches. This number is the average chance for a single peak annotation to be based on chance. For some rough guidelines, 10% or less is commonly seen for bottom up data while up to 50% can be seen for top or middle down data.
+ Intensity false match chance. This is based on the same calculation as the peaks false match chance, However this does not count the number of peaks matched but the sum intensity. This number should be similar or lower compared to the peaks false match chance.

The toggle switch above the table allows to open ion series specific statistics. These are exactly the same as above but split per ion series.

== Fragment table

This table contains all details on the spectrum in table format. It can display unassigned peaks, annotated peaks, and missing fragments. For each peak/fragment it displays all data. This whole table can be copied to other software for other analysis. Additionally, above the table is a normalised output for the ProForma definition, which removes any implementation specific notation and returns a fully specification compliant ProForma sequence.

== Settings

This section below the general statistics allows fine tuning the annotated spectrum and related sections. The graphic section allows fine control over the graphics. These boxes allow sizes set in any CSS unit. The other settings sections are detailed in their related sections above.