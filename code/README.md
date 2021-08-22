# Code
This folder contains all code used in my bachelor project. The folder layout looks like:

#### /classes
Contains the lightcurve class used as a framework for working around with lightcurves from different telescopes and observing campains.

#### /data
Contains lightcurves from RXTE and Swift in the folders **/RXTE**, **/Swift** and **/SwiftGC** (Galactic Center). The **/outbursts** folder includes csv files containing outburst fit parameters obtained within this project.

#### /notebook
Contains notebooks and light curves data files for fitting the exponetial to linear decay model. 

#### /output
All plots that were made are saved in this folder. The **/report** folder contains all final plots used in my report.

#### /plots
This folder contains all plot functions.

## User instructions
Before running any program make sure you have astropy and matplotlib installed. Go to the main **/code** folder and open you terminal. Run the following commands.

`python fit.py`
Analysis on a single source will be excecuted: lightcurve plotting, calculating background, region selection, binning, outburst selection, gaussian fitting expontial fitting and linear fitting. Currently XTEJ1728-295 is set as default but of course this can be changed into all other available light curves. Do not forget to set the telescope of the light curve.

`python analysis.py`
This will produce distributions for the duration and decay time of each outbursts per type (Black hole, Neutron star or unkown) and per source (if more than 2 outbursts).

`python fit_all.py`
For each source a lightcurve is made and all outbursts are fitted automaticly to a gaussian and exponential. 

