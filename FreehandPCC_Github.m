% Calculates the Pearson Correlation Coefficient (R) between two
% fluorescence channels. 

% First prompts the user to indicate the type and name of each image in the
% image stack. Then prompts the user to choose which image to use for
% sorting cells in the analysis.
% Next prompts the user to adjust contrast of phase-contrast image to help
% with visualization of cell boundaries.
% Then prompts the user to choose an ROI to zoom into that field of cells.
% Prompts the user to free-hand outline the outer boundary of the cell
% using a phase-contrast image.
% Allows the user to redo outline of cell until it looks good.
% By manually drawing the cell boundary, the user can select the best
% cells and not have boundaries be cut off by imperfect segmentation.
% Then a loop calculates the correlation coefficient between two
% background-corrected fluoresence images for each cell after its boundary
% is defined by the user.

% Input: .nd2 file with image stack containing two fluorescence images and
% one phase-contrast image
% Output: R statistics and box-plot showing Q1 to Q3 and median

%%
% Make FileName to save after analysis is finished.

% Indicate composite fluorescence-phase stack to be analyzed
filename1 = 'FileName';
% Indicate URL of folder containing image stack to be analyzed
directory1 = ['URL/' filename1 '.nd2'];
% Data will be saved in filename2
filename2 = [filename1 '_FreehandPCC.mat'];

% Open file directory
data1 = bfopen(directory1);

% Assign ImageInfo array to store information of image stack
% Columns 1-3 = images 1-3, respectively, in the image stack
% Row 1 = type of image (1: fluorescence, 2: phase)
% Row 2 = name of image
% Row 3 = image extracted from data1
ImageInfo = cell(3,3);

% Store images in third row of ImageInfo array
ImageInfo{3,1} = data1{1,1}{1,1};
ImageInfo{3,2} = data1{1,1}{2,1};
ImageInfo{3,3} = data1{1,1}{3,1};

% Loop to interrogate user for type and name of each image in the stack
for image_counter = 1:3
    display(horzcat('Image ', num2str(image_counter)));
    ImageInfo{1,image_counter} = input('Indicate type (1:fluorescence, 2:phase): ');
    ImageInfo{2,image_counter} = input('Indicate name: ', 's');  
end

% Assign fluorescence images to temporary variables F1temp and F2temp
% Assign phase-contrast image to permanent variable P
% Also determine prompt for display to ask user how to sort cells
if ImageInfo{1,1} == 1
    F1temp = ImageInfo{3,1};
    F1nametemp = ImageInfo{2,1};
    if ImageInfo{1,2} == 1
        F2temp = ImageInfo{3,2};
        F2nametemp = ImageInfo{2,2};
        P = ImageInfo{3,3};
        prompt_display = horzcat('1:', F1nametemp, ' 2:', F2nametemp);
    else
        F2temp = ImageInfo{3,3};
        F2nametemp = ImageInfo{2,3};
        P = ImageInfo{3,2};       
        prompt_display = horzcat('1:', F1nametemp, ' 2:', F2nametemp);
    end
else
    F1temp = ImageInfo{3,2};
    F1nametemp = ImageInfo{2,2};
    F2temp = ImageInfo{3,3};
    F2nametemp = ImageInfo{2,3};
    P = ImageInfo{3,1};
    prompt_display = horzcat('1:', F1nametemp, ' 2:', F2nametemp);
end

% Ask user to indicate the image to be used for sorting cells
disp(prompt_display);
sortID = input('Indicate image to use for sorting cells: ');

% F1 is the image to be displayed for use in sorting cells
if sortID == 1
    F1 = F1temp;
    F2 = F2temp;
    F1name = F1nametemp;
    F2name = F2nametemp;
else
    F1 = F2temp;
    F2 = F1temp;
    F1name = F2nametemp;
    F2name = F1nametemp;
end

% Subtract background from both fluorescence images
se = strel('disk',15); % for background calculation
background1 = imopen(F1,se);
F1cor = F1 - background1; % background-corrected fluorescence channel 1
background2 = imopen(F2,se);
F2cor = F2 - background2; % background-corrected fluorescence channel 2

% Initialize 
R = []; % master matrix to store data
n = 1; % index for cell count

% For R matrix: 
% 1st column = cell count
% 2nd column = cell type (1: dividing or 2: non-dividing)
% 3rd column = correlation coefficient R 
% each row = different cell

% Initialize cell array containing coordinates of cell boundary
% 1st column = cell count
% 2nd column = coordinates of cell boundary
% each row = different cell
Ioutline = cell(0,2); 
  
set(0,'DefaultFigureWindowStyle','docked'); % set default to docked display

% First prompts user to adjust contrast of phase-contrast image to
% help with visualization of cell boundary.
P2 = rescale(P); % rescale pixels to within 0-1 range for imadjust function
imshow(P2,[])
disp('Adjust contrast of image to view cell boundaries clearly')
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
    % Use the user-defined fluoresence channel F1 to sort cells
    croppedImageF = imcrop(F1cor, croprect);
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

    % Plot cropped F1 image on the left
    subplot(1,2,1)
    imshow(croppedImageF,[])
    axis on;
    title(F1name, 'FontSize', fontSize, 'Interpreter', 'None');
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
        
        % Add a row of 2 cells in boundary array to store coordinates of
        % cell boundary
        Ioutline = vertcat(Ioutline,cell(1,2));
        
        % Ask user to indicate type of cell based on sorting phenotype
        % Store in second column of matrix R
        R(n,2) = input('1:dividing cell, 2:non-dividing cell ');
        
        % Initiate conditional loop to repeat until user is satisfied with
        % outline of cell boundary
        
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

            % For selection of boundary of cell
            % Plot cropped F1 image on the left
            subplot(1,2,1)
            imshow(croppedImageF,[])
            axis on;
            title(F1name, 'FontSize', fontSize, 'Interpreter', 'None');
            xlabel('X', 'FontSize', fontSize);
            ylabel('Y', 'FontSize', fontSize);

            % Plot cropped phase-contrast image on the right
            subplot(1,2,2)
            imshow(croppedImageP,[])
            axis on;
            title('Phase-contrast', 'FontSize', fontSize, 'Interpreter', 'None');
            xlabel('X', 'FontSize', fontSize);
            ylabel('Y', 'FontSize', fontSize);

            % Ask user to outline boundary of cell
            disp('Outline cell boundary on phase-contrast image')

            % Get boundary of cell
            outperi = drawfreehand('FaceAlpha', 0, 'LineWidth', 1, 'InteractionsAllowed','none');

            ploop = input('Redo boundary of cell? y/n: ','s');
            
        end
        
        outperi_x = outperi.Position(:,1) + topleftX - 1; % account for cropped window 
        outperi_y = outperi.Position(:,2) + topleftY - 1; % due to imcrop with rect function
        
        outperi.Position(:,1) = outperi_x; % actual outlined pixels in image accounting for sliding window
        outperi.Position(:,2) = outperi_y;
        
        % Store coordinates of cell boundary in Ioutline array
        Ioutline{n,2} = [outperi_x,outperi_y];
  
        % Store index of counted cell in matrix R
        R(n,1) = n; 
        
        % Store index of counted cell in Ioutline array
        Ioutline{n,1} = n; 
        
        % Highlight cells that have now been counted on image P4
        outperi_mask = createMask(outperi,P3);
        P4 = bsxfun(@plus, P4, cast(outperi_mask, class(P4)));
        
        n = n + 1; % increase cell counter

        % Ask whether to keep field of view if there are more cells in the
        % same field to count
        j = input('keep field of view? y/n: ','s');
        
    end
    
    q = input('Continue analysis and select new ROI? y/n: ','s');
    close(gcf); % close subplots
    
end

% Loop to calculate the correlation between the two fluorescence channels
% for each cell whose boundary was defined

for a = 1:size(R,1)
          
    % Retrieve coordinates of cell boundary to crop images
    x_coor = Ioutline{a,2}(:,1);
    y_coor = Ioutline{a,2}(:,2);

    topleftX = floor(min(x_coor));
    topleftY = floor(min(y_coor));
    Xlength = ceil(max(x_coor)) - topleftX;
    Ylength = ceil(max(y_coor)) - topleftY;  

    % First fluorescence image
    % Use coordinates of cell boundary to retrieve ROI for the cell
    outer_outline1 = roipoly(F1cor,x_coor,y_coor);
    outer_mask1 = bsxfun(@times, F1cor, cast(outer_outline1, class(F1cor)));

    % Second fluorescence image
    % Use coordinates of cell boundary to retrieve ROI for the cell
    outer_outline2 = roipoly(F2cor,x_coor,y_coor);
    outer_mask2 = bsxfun(@times, F2cor, cast(outer_outline2, class(F2cor)));

    % Crop the ROI for a cell in both fluorescence channels
    croppedF1 = imcrop(outer_mask1, [topleftX, topleftY, Xlength, Ylength]);
    croppedF2 = imcrop(outer_mask2, [topleftX, topleftY, Xlength, Ylength]);
    
    R(a,3) = corr2(croppedF1,croppedF2);

end

% Compile PCC statistics for each class of cell

R_dividing = [];
R_nondividing = [];

for a = 1:size(R,1)
          
    if R(a,2) == 1 % dividing class of cells
        R_dividing = vertcat(R_dividing,R(a,3));
    else % non-dividing class of cells
        R_nondividing = vertcat(R_nondividing,R(a,3));
    end

end

% Calculate PCC statistics for each cell class if it is not empty

if size(R_dividing ~= 0)
    Median_R_dividing = median(R_dividing);
    Mean_R_dividing = mean(R_dividing);
    stdev_R_dividing = std(R_dividing);
else
end

if size(R_nondividing ~= 0)
    Median_R_nondividing = median(R_nondividing);
    Mean_R_nondividing = mean(R_nondividing);
    stdev_R_nondividing = std(R_nondividing);
else
end

% Display box-plot for dividing class of cells
boxplot(R_dividing)
xlabel('dividing cells')
ylabel('correlation coefficient')

save(filename2); 

% To view individual cells again
% n = cell number
% I = fluorescence image (F1 or F2)
% outline = roipoly(I,Ioutline{n,2}(:,1),Ioutline{n,2}(:,2));
% outline_cell = bsxfun(@times, I, cast(outline, class(I)));
% imshow(outline_cell,[])
% To outline the cells on the entire image with a dashed line
% boundary = [Ioutline{n,2}(:,1),Ioutline{n,2}(:,2)];
% imshow(I,[])
% hold on
% plot(boundary(:,1),boundary(:,2),'--r','LineWidth',3)

%% 
%%%%%%%% Written by: Truc Do and Ace George G. Santiago, Ph.D.
%%%%%%%% Date: 03/17/2019 
%%%%%%%% Walker Lab, Department of Microbiology
%%%%%%%% Harvard Medical School, Boston, MA 02115
% Copyright (c) 2019, Truc Do and Ace Santiago / Harvard Medical School
% All rights reserved.