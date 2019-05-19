%% StaphSizer Code 
% Calculates S. aureus cell size from brightfield and fluorescent images.

% Input: Brightfield image and widefield fluorescent image file from .nd2
% format

% Output: lengths and widths of S. aureus in 4 different categories: full
% septa, partial septa, no septa, and misplaced septa. Saves workspace as
% .mat in the end.

%% ---------------------- STEP 1 -----------------------------------%%
% Make sure to first set the path of MATLAB to your folder.

% Open data file containing brightfield and fluorescent images:
filename1 = 'FileName'; % brightfield image raw data file
filename2 = 'FileName'; % fluorescent image raw data file
directory1 = ['URL/' filename1 '.nd2']; 
directory2 = ['URL/' filename2 '.nd2']; 
filename3 = [filename1 '.mat']; % output file of analyzed image
data1 = bfopen(directory1); % file directory 
data2 = bfopen(directory2);

F = data2{1, 1}{1, 1}; % Fluorescence image
P = data1{1, 1}{2, 1}; % Brightfield image

% Background correction of brightfield image
P2 = imtophat(P,strel('disk',15)); % set threshold - normalizes illumination, 15 is empirically determined threshold that works well
P3= imadjust(P2); % adjust contrast
P4 = imbinarize(P3,'adaptive','Sensitivity',0.01); % binarize brightfield image - 0.01 is empirically determined sensitivity that works well

%% ----- Perform Cell Segmentation --------- %%%

% Clean up noise a bit

bw2 = bwareaopen(P4, 200); % small objects - small objects that have less than 200 white pixels nearby will be removed
D = -bwdist(~P4); % distance transform of binary image - detects the rough center of a collection of pixels
mask = imextendedmin(D,0.8); 

% Extended minima transform produces small spots that are roughly in the middle of the cells...
% to be segmented. Filter out tiny local minima using imextendedmin and...
% then modify the distance transform so that no minima occur at the filtered-out locations. 
% "minima imposition" is implemented via the function imimposemin.

D2 = imimposemin(D,mask); % modified distance transform D2
Ld = watershed(D2); % watershed transform of the modified distance transform - contains segmentation information
bw3 = bw2;
bw3(Ld == 0) = 0; % applies segmentation Ld onto the binary image bw2 to produce a segmented binary image

% Get boundary coordinates of each segmented object
B = bwboundaries(bw3,'noholes'); % B is the collection of all detected boundaries.

% Initialize cell array for sorted cells - intialize empty arrays to be
% filled as each cell is analyzed
full_sep = {};
some_sep = {};
no_sep = {};
mis_sep = {};
b = 1; % b, h, l, and a variables are initalized values for the for loops
h = 1;
l = 1;
a = 1;

% Visually inspect and sort cells according to presence or absence of septa.

F2 = F; % assign fluorescence image F to another variable F2. - need F2 variable so the original F file is not overwritten
Br = P; % assign brightfield image P to another variable Br. - need Br variable so the original P file is not overwritten
B_length = length(B); % number of B indices for the for-loop. - detects the number of cells so that we know how many times to go through the for loop

for n = 1:B_length
 
 % ---- Get coordinates of segmented image boundary for 1 cell ---- 
    
 % Boundary of cell n - determines coordinates of the boundary of each cell
 boundary = B{n}; % two-element vector (x,y)
 bound_length = length(boundary); % length of boundary
 
 % Initialize x-y coordinates of the boundary of cell n
 x_coor = zeros(); % makes a vector temporarily filled with zeros to be filled in later
 y_coor = zeros();
 
 % Get x-y coordinates of the boundary of cell n
  for i = 1:bound_length
      x_coor(i) = boundary(i,2); % vector with all the x-coordinates (2nd column) of the boundary of the cell
      y_coor(i) = boundary(i,1);
  end
  
  % Final x-y coordinates
  x_coor2 = x_coor'; % transpose - switches column and row because only column vectors are recognized
  y_coor2 = y_coor'; % transpose
  
  
  % draw boundary around brightfield and fluorescent images
  for m = 1:length(x_coor2)
    F2(y_coor2(m),x_coor2(m)) = 1; % changes pixel to black to display the trace of the boundary onto the fluorescent and brightfield images
    Br(y_coor2(m),x_coor2(m)) = 1;
  end

%%%%  Magnify cell n   %%%%%%%
  
wd = size(F2, 2); % Image width
ht = size(F2, 1); % Image height
zoom = 25; % zoom level - higher values zoom into images more

 % Get median pixel value, x-y coordinate of approximate middle of the cell.
 % Determines center of each cell to zoom into it for ease of visualization
 Y_med = median(B{n,1}(:,1));
 X_med = median(B{n,1}(:,2));

% Compute displacement:
x0 = wd/2 - zoom*X_med;
y0 = ht/2 - zoom*Y_med;

% Build transformation matrix T. - convert x-y coordinates of cell boundary
% from original image to zoomed in image
 T = [zoom   0      0; ...
      0      zoom   0; ...
      x0     y0     1];
  
tform = affine2d(T); %Needed by imwarp 
 
% J1 and J2 are the zoomed in images  
J1 = imwarp(F2, tform, 'OutputView', imref2d(size(F2)), 'Interp', 'cubic'); 
J2 = imwarp(Br, tform, 'OutputView', imref2d(size(Br)), 'Interp', 'cubic'); 
% Selected cubic interpolation. - interpolates to help you see a smooth as
% opposed to pixelated image when zoomed in

% Sort and categorize cells -----------------------------------------------

 imshowpair(J1, J2, 'montage'); % shows zoomed in images side-by-side
 title('Fluorescence                                              Brightfield');
 hold on % keeps the figure displayed during analysis
 
 i=input('1 = full, 2 = partial, 3 = none, 4 = misplaced, enter = discard: \n','s');
 
 if i == '1'
     full_sep{1,b} = B{n};
     full_sep{2,b} = num2str(n); % index of the cell in the original B variable.
     b = b + 1;
 
 elseif i == '2'
     some_sep{1,h} = B{n};
     some_sep{2,h} = num2str(n);
     h = h + 1;
 
 elseif i == '3'
     no_sep{1,l} = B{n};
     no_sep{2,l} = num2str(n);
     l = l + 1;
 
 elseif i == '4'
     mis_sep{1,a} = B{n};
     mis_sep{2,a} = num2str(n);
     a = a + 1;

 elseif i == '5'
     
     break
 
 else
     
 end
 
hold off
  
F2 = F; % re-initialize to the original fluorescent image.
Br = P; % re-initialize to the original brightfield image.

end

%% ------------------ STEP 2 ---------------------------------%%%%%%%%%%%%
% Obtain cell dimensions (width, length) using fit_ellipse.m, a Least
% Squares Fitting Algorithm.

%Full septa
if ~ isempty(full_sep) % if the full_sep category is not empty, do the following
    for f = 1:length(full_sep(1,:))
        
        cellfitXY = full_sep{1,f};
        
        ellipse_t = fit_ellipse(cellfitXY(:,1),cellfitXY(:,2)); % contains short axis and long axis measurements
        width_fulsep(f) = (ellipse_t.short_axis)*(1/15); % 1/15 is the micron:pixel ratio of the 100X oil objective used for imaging
        length_fulsep(f) = (ellipse_t.long_axis)*(1/15);
        
    end
    
else % otherwise, consider this category empty.
    width_fulsep = [];
    length_fulsep = [];
end

% No septa

if ~ isempty(no_sep)
    for f = 1:length(no_sep(1,:))
        
        cellfitXY = no_sep{1,f};
        
        ellipse_t = fit_ellipse(cellfitXY(:,1),cellfitXY(:,2));
        width_nosep(f) = (ellipse_t.short_axis)*(1/15);
        length_nosep(f) = (ellipse_t.long_axis)*(1/15);
     
    end
    
else
    width_nosep = [];
    length_nosep = [];
end

% Some Septa

if ~ isempty(some_sep)
    for f = 1:length(some_sep(1,:))
        
        cellfitXY = some_sep{1,f};
        
        ellipse_t = fit_ellipse(cellfitXY(:,1),cellfitXY(:,2));
        width_somsep(f) = (ellipse_t.short_axis)*(1/15);
        length_somsep(f) = (ellipse_t.long_axis)*(1/15);
     
    end
   
else
    width_somsep = [];
    length_somsep = [];
end

% Misplaced Septa

if ~ isempty(mis_sep)
    for f = 1:length(mis_sep(1,:))
        
        cellfitXY = mis_sep{1,f};
        
        ellipse_t = fit_ellipse(cellfitXY(:,1),cellfitXY(:,2));
        width_missep(f) = (ellipse_t.short_axis)*(1/15);
        length_missep(f) = (ellipse_t.long_axis)*(1/15);
             
    end
    
else
    width_missep = [];
    length_missep = [];
end

% Consolidate widths and lengths

if ~ isempty(length_fulsep)
    Full = [length_fulsep(1,:)' width_fulsep(1,:)']; % Full Septa
else
    Full =[];
end

if ~ isempty(length_somsep)
    Part = [length_somsep(1,:)' width_somsep(1,:)']; % Partial Septa
else
    Part = [];
end

if ~ isempty(length_nosep)
    None = [length_nosep(1,:)' width_nosep(1,:)']; % No Septa
else
    None = [];
end

if ~ isempty(length_missep)
    Mis = [length_missep(1,:)' width_missep(1,:)']; % Misplaced Septa
else
    Mis = [];
end

save(filename3); % Save workspace for this image.

clc

% Once you confirm that the .mat files are generated, then you can clear
% the workspace

% Then repeat this for each image of the same strain before combining the
% data to calculate statistics in StaphPlot.m

%% File Exchange Function Citation: 
% Ohad Gal (2003), fit_ellipse (https://www.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse),... 
% MATLAB Central File Exchange. Retrieved October 18, 2018.

%%
%%%%%%%% Written by: Ace George G. Santiago, Ph.D.
%%%%%%%% Date: 10/23/2018 
%%%%%%%% Walker Lab, Department of Microbiology 
%%%%%%%% Harvard Medical School, Boston, MA 02115
% Copyright (c) 2018, Ace Santiago / Harvard Medical School
% All rights reserved.