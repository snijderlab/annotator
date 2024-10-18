#import "functions.typ": key, button

= Interactive spectrum

== Legend

Once an annotation is done the annotated spectrum will be shown below the annotation settings. It starts with a legend of the colours. For any of these legend elements if hovered over it highlights all peaks that match. If clicked this highlighting will stay until click again. Multiple highlights can be stacked which will show any peaks matching any of the criteria.

== Peptide legend

Below the legend is an overview of the peptide(s) that are annotated. All main ion series are displayed as flags in the corner between the amino acids. Hovering highlights any fragment from that position similar to highlighting in the legend. Clicking makes the highlighting permanent in the same way. If there is multiple peptides they are prefixed by the peptide index (peptidoform index followed by a dot followed by the peptide index). Hovering/clicking highlights all peptide fragments.

In the settings section the peptide can be set to use a more compact representation. Additionally parts of the peptide (and any matching fragment) can be highlighted in a different colour.

== Spectrum

In the spectrum hovering over a peak shows the mz value as well as all annotations if there are multiple annotations for that peak. Zooming in can be done by clicking at one point to select one mz bound and move the mouse and release it at another point to set the other mz bound as well as the intensity bound, a selection window will appear to indicate the section shown after zooming in. Additionally zooming can also be done with the scroll wheel. Scrolling zooms in or out on the spot where the mouse is located. Scrolling with #key[Shift] pans the spectrum to the left or right on the mz axis. Scrolling with #key[Control] zooms in on the intensity axis. Lastly using #key[Control] with #key[+] and #key[-], zooms in or out on the middle of the spectrum, and #key[Control] with #key[0] zooms to the original zoom level of the spectrum. The distance between peaks can be annotated by clicking on one peak dragging to the next and releasing on another peak. Such labels can be removed by clicking on them or all labels can be removed at once in the settings.

=== Peaks settings 

The theoretical spectrum can be shown below the x axis. Displaying the unassigned peaks can be toggled. The colouration can be set to the ion type (default), to the peptide id, the peptidoform id, or grey (none). Also there is an option to remove all distance labels.

=== Label settings 

By default the labels for peaks are shown for any peak within the 90% of intensity, that means any peak having at least $(1.0 - 0.9) * "max_intensity"$ intensity. This threshold can be controlled in the label settings section. Showing the mz value can be controlled in the same way for peaks. The 'Manually force show' allows to show either the label or mz value for any peak. Select the mode of interest and any peak clicked will show that information. The 'Hide' mode allow to hide the information (label and mz) for any peak. The button #button[Clear] clears all manually forced labels. The checkboxes 'Show in label' allows the hiding/showing of certain parts of the label to slim down the information show. Lastly the labels can be set to be 90degrees rotated to make overlap less likely. Any of these changes will also show up in the exported image.

=== Spectrum settings

This section allows zooming to precise numbers. Additionally the number of tickmarks for the x and y axis can be set. The peaks can be set to use square root intensity instead of linear intensity to see a bigger dynamic range. The y axis can be changed to show percentages instead of raw intensities as well.

== Error graph

The graph below the spectrum shows 

== Statistics overview

== Fragment table

== Settings

This section below the general statistics allows fine tuning the annotated spectrum and related sections. The graphics section allows fine control over the graphics, these boxes allow sizes set in any CSS unit. The other settings sections are detailed in their related sections above.