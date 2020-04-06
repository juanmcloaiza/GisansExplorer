## NICOS data reduction [ongoing]:

![image](.\ScreenshotNICOSViewer.png)

  - Done:
    - Reproduce functionality from Alexandros script.
      - Fix bug when converting from detector bins to $Q_y$ and $Q_z$.
      - Optimize the plotting.
      - Optimize the integration along $x$ and $y$ axes.
      - Optimize the conversion from detector bins to $Q_y$ and $Q_z$.
    - Explored several plotting strategies: 
      - matplotlib (Easy to use, not very fast):
        - scatter
        - imshow
        - contourf
      - pyqtgraph (Faster and fancier but more complicated):
        - http://www.pyqtgraph.org/
        - Consider to use it at a later stage.
    - Show plot coordinates when hovering with the mouse over the plots.
    - Allow for interactively selecting a region of interest (ROI).
      - This automatically computes and shows the integration along both axes.
    - Allow to change the colorbar limits:
      - via the mouse wheel
      - via spinbox
    - Allow to toggle between Log and Linear scale for the colorbar.
    - Allow to load different files in the same session.
    - Allow to save the figure shown as png or pdf.
    - Allow to save the $Q_y$, $Q_z$, $I$ data as ascii tables (3 files).
        
  - To do:
    - Fix issues relating to poor drawing/scaling of the GUI on different systems.
    - Allow to sum / substract the intensities from two different files.
    - Correct the ticks and labels for the axes.
    - Consolidate a first decent version and upload it to GitLab.
    - Fix issue: resizing the window is incredibly slow (ask someone, who?).
    - Go for the next steps with Alexandros.