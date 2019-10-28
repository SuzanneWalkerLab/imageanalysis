% Compares the absolute fluorescence intensity at the septum and at the
% the periphery between two strains that were differentially labeled with
% DAPI nuclear dye, mixed, and imaged simultaneously. Capturing both
% strains in the same image allows for a fair comparison and avoids the
% need to calculate septal/peripheral fluorescence ratio.

% First prompts the user to adjust contrast of phase-contrast image to help
% with visualization of cell boundaries and septum.
% Then prompts the user to choose an ROI to zoom into that field of cells.
% Asks the user to indicate the strain background based on DAPI signal.
% Asks the user to indicate the type of cell.
% Prompts the user to free-hand outline the outer cell wall and septum.
% Erosion is used to estimate the inner cell wall.
% Allows the user to redo outline of outer cell wall boundary and septum
% until it looks good.
% By manually drawing the outer cell wall boundary, the user can select the
% best-looking cells with full septa and not have boundaries automatically
% cut off by imperfect segmenting. Drawing on the phase-contrast image
% allows the user to better define the cell boundary. But the fluorescence
% intensities are calculated using background-corrected fluorescence image.

% After the user finishes choosing the desired cells, a for loop at the end
% of the program retrieves the user-defined boundaries to extract all pixel
% intensity values and calculates the mean and median of these values. The
% retrieved coordinates of the outer cell wall boundary are also used to
% fit an ellipse onto each cell. The dimensions of the major and minor axes
% of the fitted cell are determined and used to calculate the cell volume
% and aspect ratio by approximating each cell as a prolate spheroid.
% Fluorescence intensity can be binned according to cell type as well.

% Input: .nd2 file with FDL-DAPI-phase image stack
% Output: Pixel intensities at the septum and at the periphery; cell volume
% and aspect ratio
% Additional files required to run this script: fit_ellipse.m; bfmatlab

%%
% Make FileName to save after analysis is finished.

filename1 = 'FileName'; % composite FDL-DAPI-phase image to be analyzed
directory1 = ['URL/' filename1 '.nd2']; % URL of folder containing image to be analyzed
filename2 = [filename1 '_TriPeriSep.mat'];

data1 = bfopen(directory1); % file directory 

F = data1{1, 1}{1, 1}; % Fluorescence image - FDL label
FD = data1{1, 1}{2, 1}; % Fluorescence image - DAPI label
P = data1{1, 1}{3, 1}; % Phase-contrast image

% Initialize 
Isep = []; % matrix containing septum intensity
Iperi = []; % matrix containing peripheral intensity
n = 1; % index for cell count

% For Iperi/Isep variables: 
% 1st column = cell count
% 2nd column = peripheral/septal intensity mean of all non-zero intensity pixels
% 3rd column = peripheral/septal intensity median of all non-zero intensity pixels 
% 4th column = peripheral/septal intensity mean of top 25% intensity pixels
% 5th column = peripheral/septal intensity median of top 25% intensity pixels
% 6th column = peripheral/septal intensity mean of all intensity pixels
% 7th column = peripheral/septal intensity median of all intensity pixels
% 8th column = peripheral/septal intensity mean of top 75% intensity pixels
% 9th column = peripheral/septal intensity median of top 75% intensity pixels
% 10th column = peripheral/septal intensity mean of top 50% intensity pixels
% 11th column = peripheral/septal intensity median of top 50% intensity pixels

% Initialize cell array to store drawfreehand outline
% 1st column = cell count
% 2nd column = type of cell (1:no septum, 2:partial, 3:complete, 4:misplaced)
% 3rd column = outer cell wall outline
% 4rd column = inner cell wall outline
% 5th column = septum outline
% 6th column = all pixel intensity values for cell periphery
% 7th column = all pixel intensity values for septum
% 8th column = cell size estimated as volume of prolate spheroid
% 9th column = aspect ratio of cell
% 10th column = strain background identified by DAPI stain (1:DAPI-stained, 2:no stain)
% each row = different cell
Ioutline = cell(0,10); 

% Subtract background from FDL image
se = strel('disk',15); % for background calculation
background = imopen(F,se);
F2 = F - background; 

% Subtact background from DAPI image
backgroundD = imopen(FD,se);
FD2 = FD - backgroundD;

% Merge the background-corrected FDL and DAPI images
C = imfuse(F2,FD2);

% For identification of inner cell wall boundary by erosion from
% user-defined outer cell wall boundary
se2 = strel('disk',3,6); 

set(0,'DefaultFigureWindowStyle','docked'); % set default to docked display

% First prompts user to adjust contrast of phase-contrast image to
% help with visualization of cell boundary and septum.
P2 = rescale(P); % rescale pixels to within 0-1 range for imadjust function
imshow(P2,[])
disp('Adjust contrast of image to view cell boundaries and septum clearly')
imcontrast

% Prompt user to enter in new maximum intensity value
new_max = input('Enter new maximum intensity for window: ');
% Prompt user to enter in new minimum intensity value
new_min = input('Enter new minimum intensity for window: ');

% Make contrast-adjusted phase image for display in subsequent analysis
P3 = imadjust(P2,[new_min new_max], [0 1]);

% Initialize P4 as P3 - each counted cell will be displayed via imfuse
P4 = P3; 

close(gcf); % close current figure

% Loop to interrogate each field of cell cropped by user-defined ROI
% until user decides to stop analysis

q = 'y'; % initialize q to start image analysis

while q == 'y'

    % Ask user to define ROI for cropping
    % Displays contrast-adjusted phase image with analyzed cells marked
    imshow(P4,[]) 
    Faxes = gca;
    Faxes.Visible = 'On';
    disp('Select ROI to crop and zoom into cells. Shift key + drag for square ROI.')
    croprect = getrect;
    topleftX = floor(croprect(1)); % round so that x-y translation is not shifted
    topleftY = floor(croprect(2)); % and also contains integer number of pixels
    Xlength = ceil(croprect(3));
    Ylength = ceil(croprect(4));
    croprect = [topleftX, topleftY, Xlength, Ylength]; % reset cropping rectangle
    
    % Crop the images
    croppedImageF = imcrop(F2, croprect);
    croppedImageP = imcrop(P3, croprect);
    croppedImageC = imcrop(C, croprect);
    
    % Close full phase-contrast image
    close(gcf);

    % Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Give a name to the title bar.
    set(gcf, 'Name', '', 'NumberTitle', 'Off')
    format long g;
    format compact;
    fontSize = 20;

    % Plot cropped fluorescence image on the left
    subplot(1,2,1)
    imshow(croppedImageF,[])
    axis on;
    title('Fluorescence', 'FontSize', fontSize, 'Interpreter', 'None');
    xlabel('X', 'FontSize', fontSize);
    ylabel('Y', 'FontSize', fontSize);

    % Plot cropped FDL-DAPI composite image on the right
    subplot(1,2,2)
    imshow(croppedImageC,[])
    axis on;
    title('FDL-DAPI composite', 'FontSize', fontSize, 'Interpreter', 'None');
    xlabel('X', 'FontSize', fontSize);
    ylabel('Y', 'FontSize', fontSize);

    j = input('Keep field of view? (y/n): ','s');
    
    while j == 'y'
        
        % Add a row of 10 cells in boundary array to store outlines
        Ioutline = vertcat(Ioutline,cell(1,10));
        
        % Ask user to indicate the type of cell to be analyzed
        disp('1:no septum, 2:partial septum, 3:complete septum, 4:misplaced septa');
        Ioutline{n,2} = input('Indicate the type of cell to be analyzed (1-4): ');
        
        % Ask user to indicate the strain background of the cell
        disp('1:DAPI-stained, 2:no stain');
        Ioutline{n,10} = input('Indicate if this cell is stained with DAPI (1-2): ');
        
        % Initiate conditional loop to repeat until user is satisfied with
        % outline of outer cell wall boundary
        
        ploop = 'y';
        
        while ploop == 'y'
        
            % Set up figure properties:
            % Enlarge figure to full screen.
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            % Give a name to the title bar.
            set(gcf, 'Name', '', 'NumberTitle', 'Off')
            format long g;
            format compact;
            fontSize = 20;

            % For selection of outer cell wall boundary
            % Plot cropped fluorescence image on the left
            subplot(1,2,1)
            imshow(croppedImageF,[])
            axis on;
            title('Fluorescence', 'FontSize', fontSize, 'Interpreter', 'None');
            xlabel('X', 'FontSize', fontSize);
            ylabel('Y', 'FontSize', fontSize);

            % Plot cropped phase image on the right
            subplot(1,2,2)
            imshow(croppedImageP,[])
            axis on;
            title('Phase-contrast', 'FontSize', fontSize, 'Interpreter', 'None');
            xlabel('X', 'FontSize', fontSize);
            ylabel('Y', 'FontSize', fontSize);

            % Get fluorescence at the cell periphery

            % Ask user to outline outer cell wall boundary
            disp('Outline outer cell wall boundary on phase-contrast image');

            % Get outer cell wall boundary
            outperi = drawfreehand('FaceAlpha', 0, 'LineWidth', 1, 'InteractionsAllowed','none');

            ploop = input('Redo outer cell wall boundary? (y/n): ','s');
            
        end
        
        outperi_x = outperi.Position(:,1) + topleftX - 1; % account for cropped window 
        outperi_y = outperi.Position(:,2) + topleftY - 1; % due to imcrop with rect function
        
        outperi.Position(:,1) = outperi_x; % actual outlined pixels in image accounting for sliding window
        outperi.Position(:,2) = outperi_y;
        
        % Store outline of outer cell wall boundary
        Ioutline{n,3} = [outperi_x,outperi_y];
      
        outperi_mask = createMask(outperi,F2);
        % outperi_image = bsxfun(@times, F2, cast(outperi_mask, class(F2)));
        
        % Automatic erosion to identify inner cell wall boundary
        inperi_mask = imerode(outperi_mask,se2);
        % inperi_image = bsxfun(@times, F2, cast(inperi_mask, class(F2)));
        
        % Identify x,y coordinates of inner cell wall boundary
        inperi_bound = bwboundaries(inperi_mask,'noholes'); 
        inperi_x = inperi_bound{1,1}(:,2); % second column is x coordinate
        inperi_y = inperi_bound{1,1}(:,1); % first column is y coordinate
        
        % Store outline of inner cell wall boundary
        Ioutline{n,4} = [inperi_x,inperi_y];
  
        Iperi(n,1) = n; % stores index of counted cell
        
        % Initiate conditional loop to repeat until user is satisfied with
        % outline of septum
        
        sloop = 'y';
        
        while sloop == 'y'
        
            % Set up figure properties:
            % Enlarge figure to full screen.
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            % Give a name to the title bar.
            set(gcf, 'Name', '', 'NumberTitle', 'Off')
            format long g;
            format compact;
            fontSize = 20;

            % For selection of septum
            % Plot cropped phase image on the left
            subplot(1,2,1)
            imshow(croppedImageP,[])
            axis on;
            title('Phase-contrast', 'FontSize', fontSize, 'Interpreter', 'None');
            xlabel('X', 'FontSize', fontSize);
            ylabel('Y', 'FontSize', fontSize);

            % Plot cropped fluorescence image on the right
            subplot(1,2,2)
            imshow(croppedImageF,[])
            axis on;
            title('Fluorescence', 'FontSize', fontSize, 'Interpreter', 'None');
            xlabel('X', 'FontSize', fontSize);
            ylabel('Y', 'FontSize', fontSize);

            % Get fluorescence at the septum

            % Ask user to outline septum
            disp('Outline septum on fluorescence image');

            % Get outline of septum
            sept = drawfreehand('FaceAlpha', 0, 'LineWidth', 1, 'InteractionsAllowed','none');

            sloop = input('Redo septum boundary? (y/n): ','s');
            
        end
     
        sept_x = sept.Position(:,1) + topleftX - 1; % account for cropped window
        sept_y = sept.Position(:,2) + topleftY - 1; % due to imcrop with rect function
        
        sept.Position(:,1) = sept_x; % actual outlined pixels in image accounting for sliding window
        sept.Position(:,2) = sept_y;
        
        % Store outline of septum
        Ioutline{n,5} = [sept_x,sept_y];
  
        Isep(n,1) = n; % stores index of counted cell
        
        Ioutline{n,1} = n; % stores index of cell in cell boundary array
        
        % Highlight cells that have now been counted on image P4
        P4 = bsxfun(@plus, P4, cast(outperi_mask, class(P4)));
        
        n = n + 1; % increase cell counter

        % Plot cropped FDL-DAPI composite image again on the left
        % So that the user can decide if they want to analyze more cells in
        % this field of view
        subplot(1,2,1)
        imshow(croppedImageC,[])
        axis on;
        title('FDL-DAPI composite', 'FontSize', fontSize, 'Interpreter', 'None');
        xlabel('X', 'FontSize', fontSize);
        ylabel('Y', 'FontSize', fontSize);      
             
        % Ask whether to keep field of view to analyze more cells
        j = input('Keep field of view? (y/n): ','s');
        
    end
    
    q = input('Continue analysis and select new ROI? (y/n): ','s');
    close(gcf); % close subplots
    
end

% If at least one cell was analyzed (i.e. size(Isep,1) > 0)
% Loop to retrieve all pixel intensity values of cell periphery and septum
% and to calculate cell size (volume) and aspect ratio (major axis/minor
% axis) by approximating cells as prolate spheroids

if size(Isep,1) > 0
   
    for a = 1:size(Isep,1)

        % Retrieve coordinates of outline of cell periphery and septum
        % Store as logical arrays called peri_mask and septum_mask
        outer_outline = roipoly(F2,Ioutline{a,3}(:,1),Ioutline{a,3}(:,2));
        inner_outline = roipoly(F2,Ioutline{a,4}(:,1),Ioutline{a,4}(:,2));
        peri_mask = logical(outer_outline-inner_outline);
        septum_mask = roipoly(F2,Ioutline{a,5}(:,1),Ioutline{a,5}(:,2));

        % Retrieve all pixel intensity values for cell periphery and septum
        Periphery = F2(peri_mask);
        Septum = F2(septum_mask);

        % Vectorize pixel intensity values
        Peri_vector = double(Periphery(:));
        Sept_vector = double(Septum(:));

        % Sort from high to low pixel intensity values
        Peri_vector_sorted = sort(Peri_vector,'descend');
        Sept_vector_sorted = sort(Sept_vector,'descend');

        % Calculate mean and median of all nonzero intensity pixels
        Peri_vector_nz = nonzeros(Peri_vector);
        Iperi(a,2) = mean(Peri_vector_nz);
        Iperi(a,3) = median(Peri_vector_nz);

        Sept_vector_nz = nonzeros(Sept_vector);
        Isep(a,2) = mean(Sept_vector_nz);
        Isep(a,3) = median(Sept_vector_nz);

        % Calculate mean and median of top 25% intensity pixels
        qt_peri = ceil(0.25*(length(Peri_vector_sorted)));
        Peri_vector_qt = Peri_vector_sorted(1:qt_peri,:);
        Iperi(a,4) = mean(Peri_vector_qt);
        Iperi(a,5) = median(Peri_vector_qt);

        qt_sep = ceil(0.25*(length(Sept_vector_sorted)));
        Sept_vector_qt = Sept_vector_sorted(1:qt_sep,:);
        Isep(a,4) = mean(Sept_vector_qt);
        Isep(a,5) = median(Sept_vector_qt); 

        % Calculate the mean and median of all pixel intensity values
        Iperi(a,6) = mean(Peri_vector);
        Iperi(a,7) = median(Peri_vector);

        Isep(a,6) = mean(Sept_vector);
        Isep(a,7) = median(Sept_vector);

        % Calulate the mean and median of top 75% intensity pixels
        qt_peri = ceil(0.75*(length(Peri_vector_sorted)));
        Peri_vector_qt = Peri_vector_sorted(1:qt_peri,:);
        Iperi(a,8) = mean(Peri_vector_qt);
        Iperi(a,9) = median(Peri_vector_qt);

        qt_sep = ceil(0.75*(length(Sept_vector_sorted)));
        Sept_vector_qt = Sept_vector_sorted(1:qt_sep,:);
        Isep(a,8) = mean(Sept_vector_qt);
        Isep(a,9) = median(Sept_vector_qt); 

        % Calculate the mean and median of top 50% intensity pixels 
        qt_peri = ceil(0.50*(length(Peri_vector_sorted)));
        Peri_vector_qt = Peri_vector_sorted(1:qt_peri,:);
        Iperi(a,10) = mean(Peri_vector_qt);
        Iperi(a,11) = median(Peri_vector_qt);

        qt_sep = ceil(0.50*(length(Sept_vector_sorted)));
        Sept_vector_qt = Sept_vector_sorted(1:qt_sep,:);
        Isep(a,10) = mean(Sept_vector_qt);
        Isep(a,11) = median(Sept_vector_qt);

        % Store all pixel intensity values as a matrix within Ioutline array
        % for cell periphery and septum
        Ioutline{a,6} = Peri_vector;
        Ioutline{a,7} = Sept_vector; 

        % Fit cells with an ellipse to calculate cell size
        % First extract outline of each cell
        outer_outlineX = Ioutline{a,3}(:,1);
        outer_outlineY = Ioutline{a,3}(:,2);

        % The micron:pixel ratio of the 100X oil obj on Screeny is 1:15
        ellipse_cell = fit_ellipse(outer_outlineX,outer_outlineY);
        major_axis = (ellipse_cell.long_axis)*(1/15);
        minor_axis = (ellipse_cell.short_axis)*(1/15);

        % Calculate and store the cell volume within Ioutline array
        Ioutline{a,8} = (4/3)*pi*(0.5*major_axis)*(0.5*minor_axis)^2;

        % Calculate and store the aspect ratio within Ioutline array
        Ioutline{a,9} = major_axis/minor_axis;  

    end
    
end

save(filename2); 

clear all
clc

% To view individual cells again
% n = cell number
% y = 3 for outer wall, 4 for inner wall, 5 for septum wall
% outline = roipoly(F2,Ioutline{n,y}(:,1),Ioutline{n,y}(:,2));
% outline_cell = bsxfun(@times, F2, cast(outline, class(F)));
% imshow(outline_cell,[])
% To outline the cells on the entire image with a dashed line
% boundary = [Ioutline{n,y}(:,1),Ioutline{n,y}(:,2)];
% imshow(F,[])
% hold on
% plot(boundary(:,1),boundary(:,2),'--r','LineWidth',3)

%% File Exchange Function Citation: 
% Ohad Gal (2003), fit_ellipse (https://www.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse),... 
% MATLAB Central File Exchange. Retrieved October 18, 2018.

%% 
%%%%%%%% Written by: Truc Do and Ace George G. Santiago, Ph.D.
%%%%%%%% Date: 09/01/2019
%%%%%%%% Walker Lab, Department of Microbiology and Immunobiology
%%%%%%%% Harvard Medical School, Boston, MA 02115
% Copyright (c) 2019, Truc Do and Ace Santiago / Harvard Medical School
% All rights reserved.