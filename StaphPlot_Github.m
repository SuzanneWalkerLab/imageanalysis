%% This script calculates S. aureus dimensions from brightfield and fluorescence images.

%% ---- Routine 1 -----
% After all images for one strain have been analyzed, combine all widths and lengths for each septal category

% Initialize Categories in empty vectors: these are the most important
% variables

Full_comb = [];
Part_comb = [];
None_comb = [];
Mis_comb = [];

% load(.....) #1
load('FileName.mat') % load workspace for each image
Full_comb = vertcat(Full_comb, Full); 
Part_comb = vertcat(Part_comb, Part);
None_comb = vertcat(None_comb, None); 
Mis_comb = vertcat(Mis_comb, Mis);

% load(.....) #2
load('FileName.mat') % load workspace for each image
Full_comb = vertcat(Full_comb, Full); 
Part_comb = vertcat(Part_comb, Part);
None_comb = vertcat(None_comb, None); 
Mis_comb = vertcat(Mis_comb, Mis);

% load(.....) #...

% Once all .mat have been loaded and concatenized/combined, clear everything
% except the following variables.
clearvars -except Full_comb Part_comb None_comb Mis_comb 

% Save current workspace containing Full_comb, Part_comb, None_comb, Mis_comb 
filename = 'FileName.mat'; 
save(filename) 

%% ----- Routine 2 ------- 
% Commands for plotting histograms with fit

% Percentage of each class of cells
TotalCells = length(Full_comb)+length(Mis_comb)+length(None_comb)+length(Part_comb);
Full = length(Full_comb)/TotalCells*100;
Mis = length(Mis_comb)/TotalCells*100;
None = length(None_comb)/TotalCells*100;
Part = length(Part_comb)/TotalCells*100;

% Aspect Ratio
AR_None = None_comb(:,1)./None_comb(:,2); 
AR_Part = Part_comb(:,1)./Part_comb(:,2); 
% histfit(AR_Part) %to obtain figure showing frequency of cells vs. cell size
% distributionFitter(AR_None) %to obtain the mean and variance of aspect ratio

% To plot histogram

% Volume
Vol_None = ((4/3)*pi).*(None_comb(:,1)/2).*(None_comb(:,2)/2).^2;
Vol_Part = ((4/3)*pi).*(Part_comb(:,1)/2).*(Part_comb(:,2)/2).^2;
% distributionFitter(Vol_Part)
% distributionFitter(Vol_None)
% histfit(Vol_None)

save(filename);


%%
%%%%%%%% Written by: Ace George G. Santiago, Ph.D.
%%%%%%%% Date: 10/23/2018 
%%%%%%%% Walker Lab, Department of Microbiology 
%%%%%%%% Harvard Medical School, Boston, MA 02115
% Copyright (c) 2018, Ace Santiago / Harvard Medical School
% All rights reserved.