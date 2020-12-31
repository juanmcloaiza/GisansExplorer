.. include:: Figures.rst

How to...
=========

Open a first data file 
**********************

 - By navigating to it on the bottom left panel:
    |NavigateToFile|

 - By browsing to it using the `open...` button on the bottom rigth:
    |OpenButton|
    |FileBrowser|

Select the region of interest (ROI)
***********************************

- By using the ROI selector on the left plot:
    - Click a corner and drag **to resize**.  
    - Click the center and drag **to move**.

    |SquareROISelector|

- By setting the limits on the information table on the right:
    - Set the pixel limit values `x0`, `y0`, `xf`, `xf` (left, bottom, right, top).
    - After they turn green --wrong entries will turn red, click `Update`.

    |TableROISettings|

Take a closer look to the ROI plots
***********************************

- By double clicking on any of the ROI plots, |DoubleClickRegion| ,three pop-up windows will appear:

    - :math:`Q_z` vs integration along :math:`Q_y`

    |IntegrationAlongQy|

    - :math:`Q_y` vs integration along :math:`Q_z`

    |IntegrationAlongQz|

    - 2D Intensity map as function of :math:`Q_y, Q_z`

    |GisansMap|

    - 3D Intensity surface as function of :math:`Q_y, Q_z`

    |GisansSurface|

Save individual figures
***********************

Refer to the step above and click the floppy disk icon in the desired figure.

Add or subtract intensities from different Gisans maps
******************************************************

    - After loading one or more files, click the corresponding maps on the gisans map list. To select several files, hold the keys ``shift`` or ``ctrl`` while clicking on each entry. This will automatically show the **addition** of the intensities of the selected gisans maps. 

    |AddIntensities|

    - To show the **subtraction** of the intensities, select any two entries on the list and tick the checkbox *Subtract intensities*

    |SubtractIntensities|

    N.B. Addition is calculated in the obvious way, :math:`I = I_A + I_B`, while subtraction is calculated as :math:`I = |I_A - I_B|`


Save the ROI as ascii column files
**********************************

    