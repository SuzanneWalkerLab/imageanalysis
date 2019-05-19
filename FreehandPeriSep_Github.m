% Calculates the ratio of fluorescence intensity at the septum vs. at 
% the periphery for Staphylococcus aureus cells using free-hand input
% from user to identify boundaries. 

% First prompts the user to adjust contrast of phase-contrast image to help
% with visualization of cell boundaries and septum.
% Then prompts the user to choose an ROI to zoom into that field of cells.
% Prompts the user to free-hand outline the septum and outer cell wall.
% Dilation is used to estimate the inner cell wall.
% Allows the user to redo outline of outer cell wall boundary and septum
% until it looks good.
% By manually drawing the outer cell wall boundary, the user can select the
% best-looking cells with full septa and not have boundaries automatically
% cut off by imperfect segmenting. Drawing on the phase-contrast image
% allows the user to better define the cell boundary. But the fluorescence
% intensities are calculated using background-corrected fluorescence image.

% Because I did not originally calculate the mean and median of all
% intensity pixels without discarding zero pixels, at the end of the code
% after all cells are analyzed, a for loop is run to calculate these
% values. This for loop will also calculate the mean and median for the top
% 75% and 50% pixels.This for loop will also retrieve and store all
% intensity pixels for the cell periphery and septum in an array. 

% Input: .nd2 file with fluorescence-phase image stack
% Output: Ratio of Mean Septal Fluorescence over Peripheral Fluorescence

%%
% Make FileName to save after analysis is finished.

filename1 = 'FileName'; % composite fluorescence-phase image to be analyzed
directory1 = ['URL/' filename1 '.nd2']; %URL of folder containing image to be analyzed
filename2 = [filename1 '_FreehandPeriSep.mat'];

data1 = bfopen(directory1); % file directory 

F = data1{1, 1}{1, 1}; % Fluorescence image
P = data1{1, 1}{2, 1}; % Phase-contrast image

% Initialize 
Isep = []; % matrix containing septum intensity
Iperi = []; % matrix containing peripheral intensity
n = 1; % index for cell count

% For Iperi/Isep variables: 
% 1st column=cell count
% 2nd column= peripheral/septal intensity mean of all non-zero intensity pixels
% 3rd column=peripheral/septal intensity median of all non-zero intensity pixels 
% 4th column=peripheral/septal intensity mean of top 25% intensity pixels
% 5th column=peripheral/septal intensity median of top 25% intensity pixels
% 6th column=peripheral/septal intensity mean of all intensity pixels
% 7th column=peripheral/septal intensity median of all intensity pixels
% 8th column=peripheral/septal intensity mean of top 75% intensity pixels
% 9th column=peripheral/septal intensity median of top 75% intensity pixels
% 10th column=peripheral/septal intensity mean of top 50% intensity pixels
% 11th column=peripheral/septal intensity median of top 50% intensity pixels

% Initialize cell array in which each cell contains drawfreehand outline
% 1st column = cell count
% 2nd column = outer cell wall outline
% 3rd column = inner cell wall outline
% 4th column = septum outline
% 5th column = all pixel intensity values for cell periphery
% 6th column = all pixel intensity values for septum
% each row = different cell
Ioutline = cell(0,4); 

% Subtract background
se = strel('disk',15); % for background calculation
background = imopen(F,se);
F2 = F - background; 

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
    
    % Crop the image
    croppedImageF = imcrop(F2, croprect);
    croppedImageP = imcrop(P3, croprect);
    
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

    %Plot cropped fluorescence image on the left
    subplot(1,2,1)
    imshow(croppedImageF,[])
    axis on;
    title('Fluorescence', 'FontSize', fontSize, 'Interpreter', 'None');
    xlabel('X', 'FontSize', fontSize);
    ylabel('Y', 'FontSize', fontSize);

    %Plot cropped phase image on the right
    subplot(1,2,2)
    imshow(croppedImageP,[])
    axis on;
    title('Phase-contrast', 'FontSize', fontSize, 'Interpreter', 'None');
    xlabel('X', 'FontSize', fontSize);
    ylabel('Y', 'FontSize', fontSize);

    j = input('Keep field of view? y/n: ','s');
    
    while j == 'y'
    
        % close(gcf); % close subplots
        
        % Add a row of 4 cells in boundary array to store outlines
        Ioutline = vertcat(Ioutline,cell(1,4));
        
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
            disp('Outline outer cell wall boundary on phase-contrast image')

            % Get outer cell wall boundary
            outperi = drawfreehand('FaceAlpha', 0, 'LineWidth', 1, 'InteractionsAllowed','none');

            ploop = input('Redo outer cell wall boundary? y/n: ','s');
            
        end
        
        outperi_x = outperi.Position(:,1) + topleftX - 1; % account for cropped window 
        outperi_y = outperi.Position(:,2) + topleftY - 1; % due to imcrop with rect function
        
        outperi.Position(:,1) = outperi_x; % actual outlined pixels in image accounting for sliding window
        outperi.Position(:,2) = outperi_y;
        
        % Store outline of outer cell wall boundary
        Ioutline{n,2} = [outperi_x,outperi_y];
      
        outperi_mask = createMask(outperi,F2);
        outperi_image = bsxfun(@times, F2, cast(outperi_mask, class(F2)));
        
        % Automatic erosion to identify inner cell wall boundary
        
        % Erode to obtain inner cell wall boundary
        inperi_mask = imerode(outperi_mask,se2);
        inperi_image = bsxfun(@times, F2, cast(inperi_mask, class(F2)));
        
        % Identify x,y coordinates of inner cell wall boundary
        inperi_bound = bwboundaries(inperi_mask,'noholes'); 
        inperi_x = inperi_bound{1,1}(:,2); % second column is x coordinate
        inperi_y = inperi_bound{1,1}(:,1); % first column is y coordinate
        
        % Store outline of inner cell wall boundary
        Ioutline{n,3} = [inperi_x,inperi_y];
        
        % Obtain intensity of cell periphery by subtracting boundaries
        peri_image = imsubtract(outperi_image,inperi_image);
        peri_mask = logical(outperi_mask-inperi_mask);
        Periphery = F2(peri_mask); % returns pixel intensities where logical = 1
  
        Iperi(n,1) = n; % stores index of counted cell
       
        % Matrix containing intensity pixels for the cell periphery
        Peri_vector = double(Periphery(:)); %vectorize
        Peri_vector_sorted = sort(Peri_vector,'descend'); %sort from high to low intensity
        Peri_vector_nz = nonzeros(Peri_vector_sorted);
        % To calculate mean and median of all nonzero intensity pixels
        Iperi(n,2) = mean(Peri_vector_nz);
        Iperi(n,3) = median(Peri_vector_nz);
        % To calculate mean and median of top 25% intensity pixels
        qt_peri = ceil(0.25*(length(Peri_vector_sorted)));
        Peri_vector_qt = Peri_vector_sorted(1:qt_peri,:);
        Iperi(n,4) = mean(Peri_vector_qt);
        Iperi(n,5) = median(Peri_vector_qt);
        
        % close(gcf); % close subplots before replotting for septum selection
        
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
            disp('Outline septum on fluorescence image')

            % Get outline of septum
            sept = drawfreehand('FaceAlpha', 0, 'LineWidth', 1, 'InteractionsAllowed','none');

            sloop = input('Redo septum boundary? y/n: ','s');
            
        end
     
        sept_x = sept.Position(:,1) + topleftX - 1; % account for cropped window
        sept_y = sept.Position(:,2) + topleftY - 1; % due to imcrop with rect function
        
        sept.Position(:,1) = sept_x; % actual outlined pixels in image accounting for sliding window
        sept.Position(:,2) = sept_y;
        
        % Store outline of septum
        Ioutline{n,4} = [sept_x,sept_y];
  
        septum_mask = createMask(sept,F2);
        septum_image = bsxfun(@times, F2, cast(septum_mask, class(F2)));
        Septum = F2(septum_mask); % returns pixel intensities where logical = 1
        
        Isep(n,1) = n; % stores index of counted cell
        
        %Get mean septum fluorescence
        % Matrix containing intensity pixels for the septum
        Septum_vect = double(Septum(:)); % vectorize
        Septum_vect_sorted = sort(Septum_vect,'descend');
        Septum_vect_nz = nonzeros(Septum_vect_sorted);
        % To calculate mean and median of all nonzero intensity pixels
        Isep(n,2) = mean(Septum_vect_nz);
        Isep(n,3) = median(Septum_vect_nz);
        % To calculate mean and median of top 25% intensity pixels
        qt_sep = ceil(0.25*(length(Septum_vect_sorted)));
        Septum_vect_qt = Septum_vect_sorted(1:qt_sep,:);
        Isep(n,4) = mean(Septum_vect_qt);
        Isep(n,5) = median(Septum_vect_qt);
        
        Ioutline{n,1} = n; % stores index of cell in cell boundary array
        
        % Highlight cells that have now been counted on image P4
        P4 = bsxfun(@plus, P4, cast(outperi_mask, class(P4)));
        
        n = n + 1; % increase cell counter

        % Ask whether to keep field of view if there are more cells in the
        % same field to count
        j = input('keep field of view? y/n: ','s');
        
    end
    
    q = input('Continue analysis and select new ROI? y/n: ','s');
    close(gcf); % close subplots
    
end

% Loop to retrieve all pixel intensity values of cell periphery and septum

% Add 2 columns to cell array to store all pixel intensity values
Ioutline = horzcat(Ioutline,cell(size(Isep,1),2));

for a = 1:size(Isep,1)
    
    % Retrieve coordinates of outline of cell periphery and septum
    % Store as logical arrays called mask2
    outer_outline = roipoly(F2,Ioutline{a,2}(:,1),Ioutline{a,2}(:,2));
    inner_outline = roipoly(F2,Ioutline{a,3}(:,1),Ioutline{a,3}(:,2));
    peri_mask2 = logical(outer_outline-inner_outline);
    septum_mask2 = roipoly(F2,Ioutline{a,4}(:,1),Ioutline{a,4}(:,2));
    
    % Retrieve all pixel intensity values for cell periphery and septum
    Periphery2 = F2(peri_mask2);
    Septum2 = F2(septum_mask2);
    
    % Calculate the mean and median of all pixel intensity values
    Peri_vector2 = double(Periphery2(:)); %vectorize
    Septum_vect2 = double(Septum2(:)); % vectorize
    
    Iperi(a,6) = mean(Peri_vector2);
    Iperi(a,7) = median(Peri_vector2);
    
    Isep(a,6) = mean(Septum_vect2);
    Isep(a,7) = median(Septum_vect2);
    
    % Calulate the mean and median of top 75% intensity pixels
    
    Peri_vector2_sorted = sort(Peri_vector2,'descend'); %sort from high to low intensity
    qt_peri2 = ceil(0.75*(length(Peri_vector2_sorted)));
    Peri_vector2_qt = Peri_vector2_sorted(1:qt_peri2,:);
    Iperi(a,8) = mean(Peri_vector2_qt);
    Iperi(a,9) = median(Peri_vector2_qt);
    
    Septum_vect2_sorted = sort(Septum_vect2,'descend'); %sort from high to low intensity
    qt_sep2 = ceil(0.75*(length(Septum_vect2_sorted)));
    Septum_vect2_qt = Septum_vect2_sorted(1:qt_sep2,:);
    Isep(a,8) = mean(Septum_vect2_qt);
    Isep(a,9) = median(Septum_vect2_qt); 
    
    % Calculate the mean and median of top 50% intensity pixels
    
    qt_peri2 = ceil(0.50*(length(Peri_vector2_sorted)));
    Peri_vector2_qt = Peri_vector2_sorted(1:qt_peri2,:);
    Iperi(a,10) = mean(Peri_vector2_qt);
    Iperi(a,11) = median(Peri_vector2_qt);
 
    qt_sep2 = ceil(0.50*(length(Septum_vect2_sorted)));
    Septum_vect2_qt = Septum_vect2_sorted(1:qt_sep2,:);
    Isep(a,10) = mean(Septum_vect2_qt);
    Isep(a,11) = median(Septum_vect2_qt);
    
    % Store all pixel intensity values as a matrix within Ioutline array
    % for cell periphery and septum
    Ioutline{a,5} = Peri_vector2;
    Ioutline{a,6} = Septum_vect2; 
    
end
    
F_ratio = mean(Isep(:,7)./Iperi(:,7));

save(filename2); 

% To view individual cells again
% n = cell number
% y = 2 for outer wall, 3 for inner wall, 4 for septum wall
% outline = roipoly(F2,Ioutline{n,y}(:,1),Ioutline{n,y}(:,2));
% outline_cell = bsxfun(@times, F2, cast(outline, class(F)));
% imshow(outline_cell,[])
% To outline the cells on the entire image with a dashed line
% boundary = [Ioutline{n,y}(:,1),Ioutline{n,y}(:,2)];
% imshow(F,[])
% hold on
% plot(boundary(:,1),boundary(:,2),'--r')

%% 
%%%%%%%% Written by: Truc Do and Ace George G. Santiago, Ph.D.
%%%%%%%% Date: 01/17/2019 
%%%%%%%% Walker Lab, Department of Microbiology
%%%%%%%% Harvard Medical School, Boston, MA 02115
% Copyright (c) 2019, Truc Do and Ace Santiago / Harvard Medical School
% All rights reserved.