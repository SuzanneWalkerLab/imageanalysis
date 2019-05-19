# SaureusImageAnalysis
MATLAB scripts to analyze wide-field fluorescence images of Staphylococcus aureus cells

Use version MATLAB _R2018b

StaphSizer
	Main script to automatically segment cells and define cell boundaries using brightfield images. The user is prompted to sort cells according to septal phenotypes using fluorescence images. Identified cells are then fitted with an ellipse to calculate their dimensions.
	Input: Brightfield and wide-field fluorescence images saved in separate image stacks in .nd2 format
	Output: Dimensions of cells in a .mat file
	Additional files: Download the package “bfmatlab.zip” from Bioformats in order to open .nd2 files, and fit_ellipse.m.

StaphPlot

Supporting script that collects and averages dimensions of cells identified in separate StaphSizer runs to calculate the average cell volume and aspect ratio.

Input: Dimensions of cells in a .mat output file from StaphSizer

Output: Average volume and aspect ratio of cells in a .mat file

FreehandPeriSep

Script to calculate the ratio of fluorescence intensity at the cell periphery versus at the septum. The user is prompted to select cells with a single complete septum and to trace the outer boundary of the cell and its septum on a phase-contrast image.

Input: Phase-contrast and wide-field fluorescence images saved as an image stack in .nd2 format

Output: Ratio of fluorescence intensity at the periphery versus at the septum

Additional files: Download the package “bfmatlab.zip” from Bioformats in order to open .nd2 files

FreehandPCC

Script to calculate the Pearson correlation coefficient (PCC) between two fluorescence channels to measure the colocalization of two fluorescent fusion proteins. The user is prompted to trace the outer boundary of cells on a phase-contrast image. These identified cells are used for PCC calculations.

Input: An image stack in .nd2 format containing two wide-field fluorescence images from different channels and one phase-contrast image

Output: PCC between two fluorescence channels with the option to display a boxplot summary

Additional files: Download the package “bfmatlab.zip” from Bioformats in order to open .nd2 files
