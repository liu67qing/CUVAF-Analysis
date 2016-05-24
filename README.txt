CUVAF Analysis version 1.0
Matlab © (R2011b) 1994 - 2011 The MathWorks, Inc.

Emily Huynh
Lions Eye Institute
emilyhuynh@lei.org.au

May 2016


Description
------------

Semi-automatic method for analyzing the region of sun-related changes in the eye in Conjunctival Ultraviolet 
Autofluorescence (CUVAF). GUI enables user to browse for images, mark the region of fluorescence and refine 
the boundary. Area measurements and parameters to enable retrospective review and post-analysis modification
are saved in an excel spreadsheet. 

The calibration factor is currently set at 3008 pixels = 24 mm which was determined by photographing a ruler. This
should be adjusted according to settings used to aquire CUVAF images.


License
-------

See file LICENSE.txt.



Files
-----

Files Included in CUVAF Analysis source code:
	CUVAF_GUI.m			Refine_UVAF.m				poi_library
	CUVAF_GUI.fig			RemoveRegions.m				CUVAF_Analysis.exe
	Analyze_CUVAF.m			Reproduce_Results.m
	LocalThresh.m			xlwrite.m
	

Data organization for CUVAF analysis
--------------------------------------

CUVAF images are organized into folders labelled with the subject ID. To ensure data is saved correctly
in the corresponding excel spreadsheet, CUVAF images are to be labelled as follows:
	OD Nasal UVAF
	OD Temporal UVAF
	OS Nasal UVAF
	OS Temporal UVAF
If colour images are available for the subject, these are also saved in the same folder with the labels:
	OD Iris Colour
	OD Nasal Colour
	OD Temporal Colour
	OS Iris Colour
	OS Nasal Colour
	OS Temporal Colour


Instructions for source code - for users with a Matlab license
----------------------------
1. Open CUVAF_GUI.m in Matlab.
2. Click 'run' to launch the GUI


Instructions for executable - for users without a Matlab license
---------------------------
* Users of the Application must accept the terms of the MCR Library License prior to installation

1. Download Executable.zip from https://github.com/hewittlab/CUVAF-Analysis/releases
2. Unzip the folder and run CUVAF_Analysis_pkg.exe to extract Matlab Compiler Runtime from package
3. Install MRCInstaller.exe 
4. Run CUVAF_Analysis.exe

