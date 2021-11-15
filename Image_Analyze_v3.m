%Leaf Image Segmentation Analysis Tool (LISAT) 
%Program Created by: Jose M. Lazaro
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 20;
captionFontSize = 14;
% Open an image.
% Browse for the image file. 
[baseFileName, folder] = uigetfile('*.*', 'Specify an image file'); 
fullImageFileName = fullfile(folder, baseFileName); 
if folder == 0
  return;
end
%create new folder
fprintf('Leaf Analysis: initiated...\n')

image_name= strrep(baseFileName,'.jpg','');
mkdir([image_name]);

cd([image_name]);

RGB = imread(fullImageFileName);
%%
% ADD pixels border to the image to correct border detection issue
fprintf('Adding pixels bordering to orignal image\n')
border_color = [55 122 200];
border_pixels_width = 5;
%%   
RGB = addborder(RGB, border_pixels_width, border_color, 'outer');
%%
RGB = imgaussfilt(RGB);
%figure()
BW1 = imbinarize(RGB);
BW1 = bwareaopen(BW1, 1500);
BW2 = bwperim(BW1);

imwrite(BW2,'1_Binarized.jpg','jpg');
fprintf('STEP1 finalized...\n')
%% 
folder = pwd; % Determine where working folder is (works with all versions).
baseFileName = '1_Binarized.jpg';
fullFileName = fullfile(folder, baseFileName); 
Binarized_image = imread(fullFileName);

% Read in original RGB image.
rgbImage = Binarized_image;
% Extract color channels.
blueChannel = rgbImage(:,:,3); % Blue channel
% Create an all black channel.
allBlack = zeros(size(rgbImage, 1), size(rgbImage, 2), 'uint8');
% Create color versions of the individual color channels.
just_blue = cat(3, allBlack, allBlack, blueChannel);

% To black and white
redChannel = just_blue(:,:,1);
greenChannel = just_blue(:,:,2);
blueChannel = just_blue(:,:,3);
blue_thresholdValue = 200; % Whatever you define white as.
BluePixels = redChannel < 1 & greenChannel < 1 & blueChannel > blue_thresholdValue;
redChannel(BluePixels) = 255;
greenChannel(BluePixels) = 255;
blueChannel(BluePixels) = 255;
newBWImage = cat(3, redChannel, greenChannel, blueChannel);

imwrite(newBWImage,'2_bw.jpg','jpg');
fprintf('STEP2 finalized...\n')
%%
originalImage = RGB;
grayImage = rgb2gray(newBWImage);
binaryImage = im2bw(grayImage);

% Display the grayscale image.
subplot(2, 3, 1);
imshow(originalImage);
% Maximize the figure window.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
drawnow;
caption = sprintf('Original "leafs" image showing');
title(caption, 'FontSize', captionFontSize);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.

%%
%Create a mask
BW_mask = binaryImage;
BW_mask(:,:,2) = BW_mask;
BW_mask(:,:,3) = BW_mask(:,:,1);

ROI = originalImage;
ROI(BW_mask == 1) = 255;
imwrite(ROI,'3_masked.jpg','jpg');
fprintf('STEP3 finalized...\n')

% Display the binary image.
subplot(2, 3, 2);
imshow(ROI); 
title('Image, obtained by masking', 'FontSize', captionFontSize); 
%%
%remove traces of blue from overfitting mask
% To black and white
redChannel = ROI(:,:,1);
greenChannel = ROI(:,:,2);
blueChannel = ROI(:,:,3);
    red_thresholdValue = 65; %average values for thresholding
    green_thresholdValue = 105;
    blue_thresholdValue = 80;

BluePixels = redChannel < red_thresholdValue & greenChannel < green_thresholdValue & blueChannel > blue_thresholdValue;
redChannel(BluePixels) = 255;
greenChannel(BluePixels) = 255;
blueChannel(BluePixels) = 255;
masked_Image_blue_removed = cat(3, redChannel, greenChannel, blueChannel);
imwrite(masked_Image_blue_removed,'4_masked_overffitted_rm.jpg','jpg');
fprintf('STEP4 finalized...\n')
%%
%Turn into trinary image:
% Step1: starting by white
redChannel = masked_Image_blue_removed(:,:,1);
greenChannel = masked_Image_blue_removed(:,:,2);
blueChannel = masked_Image_blue_removed(:,:,3);
    red_thresholdValue = 254; %average values for thresholding
    green_thresholdValue = 254;
    blue_thresholdValue = 254;

WhitePixels = redChannel > red_thresholdValue & greenChannel > green_thresholdValue & blueChannel > blue_thresholdValue;
redChannel(WhitePixels) = 100;%Gray is 128 in RGB
greenChannel(WhitePixels) = 100;
blueChannel(WhitePixels) = 100;
masked_Image_Trinary_Step1 = cat(3, redChannel, greenChannel, blueChannel);

% Step2: Green to real green
redChannel = masked_Image_Trinary_Step1(:,:,1);
greenChannel = masked_Image_Trinary_Step1(:,:,2);
blueChannel = masked_Image_Trinary_Step1(:,:,3);
    red_thresholdValue = 65; %average values for thresholding
    green_thresholdValue = 100;
    blue_thresholdValue = 60;

WhitePixels = redChannel < red_thresholdValue & greenChannel < green_thresholdValue & blueChannel < blue_thresholdValue;
redChannel(WhitePixels) = 0;%Green is 0,128,0 in RGB
greenChannel(WhitePixels) = 128;
blueChannel(WhitePixels) = 0;
masked_Image_Trinary_Step2 = cat(3, redChannel, greenChannel, blueChannel);

% Step3: Green to real green
redChannel = masked_Image_Trinary_Step2(:,:,1);
greenChannel = masked_Image_Trinary_Step2(:,:,2);
blueChannel = masked_Image_Trinary_Step2(:,:,3);
    red_thresholdValue = 90; %average values for thresholding
    green_thresholdValue = 100;
    blue_thresholdValue = 90;

WhitePixels = redChannel > red_thresholdValue & greenChannel > green_thresholdValue & blueChannel > blue_thresholdValue;
redChannel(WhitePixels) = 255;%White is 255,255,255 in RGB
greenChannel(WhitePixels) = 255;
blueChannel(WhitePixels) = 255;
masked_Image_Trinary_Step3 = cat(3, redChannel, greenChannel, blueChannel);
imwrite(masked_Image_Trinary_Step3,'5_Trinary.jpg','jpg');
fprintf('STEP5 finalized...\n')

% Show trinary conversion.
subplot(2, 3, 3);
imshow(masked_Image_Trinary_Step3);  % Show the gray scale image.
title('Trinary Image', 'FontSize', captionFontSize);
%%
% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.

Boundaries_Img=subplot(2, 3, 4);
imshow(masked_Image_Trinary_Step3);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'red', 'LineWidth', 0.5);
end
hold off;
set(Boundaries_Img, 'Units','Normalized','OuterPosition',[0 0 1 1]);
imwrite(getframe(Boundaries_Img).cdata, '6_Boundaries_included.tif','tif')

subplot(2, 3, 4)
I = imread('6_Boundaries_included.tif');
title('Outlines, leaf boundaries', 'FontSize', captionFontSize); 
imshow(I);

close;

fprintf('STEP6 finalized...\n')

%saveas(Boundaries_Img,'1_Boundaries_included.tif','tif');
%%
b=bwboundaries(binaryImage);   % BW being your binary image
mask=false(size(binaryImage));
for i = 1:length(b)
    for j = 1:length(b{i})
        ind = b{i}(j,:);
        mask(ind(1),ind(2))=1;
    end
end
stats = regionprops(mask,'BoundingBox','Area','Centroid','Perimeter');

T = struct2table(stats); % convert the struct array to a table
sortedT = sortrows(T, 'Area','descend'); % sort the table by 'Area Size'
filterT = sortedT.Area > 8000 & sortedT.Area < 48000;
sortedT= sortedT(filterT,:);

%Sorting Leaves by y-centroid parameter
allCentroids = sortedT.Centroid;
xCentroids = allCentroids(:, 1);
yCentroids = allCentroids(:, 2);
sortedT.Centroid = yCentroids;
sortedT.Perimeter = xCentroids;
sortedT = sortrows(sortedT,'Centroid','ascend') % sort the table by 'Area Size'

sorted_stats = table2struct(sortedT); % change it back to struct array if necessary
numberOfLeafs = size(sorted_stats, 1);
%%
%signaling the leafs
textFontSize = 60;	% Used to control size of "number" labels put atop the image.
labelShiftX = -3;	% Used to align the labels in the centers of the leaves.
Boundaries_Img = figure('Menu','none','ToolBar','none');
    imshow(originalImage);

for k = 1 : numberOfLeafs           % Loop through all blobs.
	% Find the mean of each blob.  (R2008a has a better way where you can pass the original image
	% directly into regionprops.  The way below works for all versions including earlier versions.)
	leafCentroidy = sortedT.Centroid(k);		% Get centroid one at a time
    leafCentroidx = sortedT.Perimeter(k);
	% Put the "leaf number" labels on the "boundaries" of the original image.
	text(leafCentroidx + labelShiftX, leafCentroidy, num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','red');
end
fprintf('STEP6.1: Leaf identification finalized...\n')

saveas(Boundaries_Img,'6.1_Leaves_order.jpg','jpg');

%%
fprintf('STEP7: leaves segmentation initiated...\n')
mkdir leaves
cd leaves
%Writing images in leaves folder
if ~isempty(sorted_stats)
    for i=1:numberOfLeafs
            bb = sorted_stats(i).BoundingBox;
            obj = imcrop(masked_Image_Trinary_Step3,bb);
            %figure, imshow(obj);
            imwrite(obj,strcat('leaf_',string(i),'.tif'),'tif');
            
            obj = imcrop(originalImage,bb);
            %figure, imshow(obj);
            imwrite(obj,strcat('leaf_',string(i),'.jpg'),'jpg');

    end
end
cd ..
fprintf('STEP7: leaves segmentation done...\n')
%%
fprintf('STEP8: Results analysis initiated...\n')
mkdir results
% Read the directory to get a list of images.
imagesFolder = 'leaves';
filePattern = [imagesFolder, '\*.jpg'];
jpegFiles = dir(filePattern);

imageFiles = jpegFiles;

% Bail out if there aren't any images in that folder.
numberOfImagesProcessed = 0;
numberOfImagesToProcess = length(imageFiles);
if numberOfImagesToProcess <= 0
 	message = sprintf('No leaf images in the folder\n%s\nClick OK to Exit.', imagesFolder);
	uiwait(msgbox(message));
	return;
end

% Create a figure for our images.
figure;
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gcf,'name','Image Analysis Demo','numbertitle','off') 

% Preallocate arrays to hold the mean intensity values of all the images.
redChannel_Mean = zeros(numberOfImagesToProcess, 1);
greenChannel_Mean = zeros(numberOfImagesToProcess, 1);
blueChannel_Mean = zeros(numberOfImagesToProcess, 1);

% We'll be skipping monochrome and indexed images
% and just looking at true color images.
% Keep track of how many we actually look at.
numberOfImagesToProcess2 = numberOfImagesToProcess;

% Loop though all images, calculating and displaying the histograms.
% and then getting the means of the Red, green, and blue channels.
for k = 1 : numberOfImagesToProcess
	% Read in this one file.
	baseFileName = imageFiles(k).name;
	fullFileName = fullfile(imagesFolder, baseFileName);
	rgbImage = imread(fullFileName);
	[rows, columns, numberOfColorBands] = size(rgbImage);
    
    %read tif image
    tif_image_name= strrep(baseFileName,'.jpg','');
    tif_image_name = strcat(tif_image_name,'.tif');
    tif_fullFileName = fullfile(imagesFolder, tif_image_name);
    tifImage = imread(tif_fullFileName);
    tifImage = rgb2gray(tifImage);
    
    %Check image content of pixels

    numberOfPixels = numel(tifImage);
    whitePixels = tifImage > 150; % Produces a binary (logical) image.
        count_white = sum(whitePixels(:));
    backgroundPixels = tifImage == 100; % Produces a binary (logical) image.
         count_background = sum(backgroundPixels(:));
    greenPixels = tifImage < 100; % Produces a binary (logical) image.
        count_green = sum(greenPixels(:));
    %leaf area calculation
    leaf_area_pixels = numberOfPixels -  count_background;
    undetermined_color_area = leaf_area_pixels - count_white - count_green;
    %error rate determination
    error_rate = (undetermined_color_area/numberOfPixels) * 100;

    %Percentage calculation
    Green_percentage = (count_green/(count_white + count_green))*100;
    White_percentage = (count_white/(count_white + count_green))*100;

	% Check to see that it is a color image (3 dimensions).
	% Skip it if it is not true RGB color.
	if numberOfColorBands < 3
		% Skip monochrome or indexed images.
		fprintf('Skipped %s.  It is a grayscale or indexed image.\n', baseFileName);
		% Decrement the number of images that we'll report that we need to look at.
		numberOfImagesToProcess2 = numberOfImagesToProcess2 - 1;  
		continue;
	end
	
	% If we get to here, it's a true color image.
	subplot(3, 3, 1);
	imshow(rgbImage, []);
	
	% Create a title for the image.
	caption = sprintf('Original Color Image\n%s\n%d rows by %d columns by %d color channels', ...
		baseFileName, rows, columns, numberOfColorBands);
	% If there are underlines in the name, title() converts the next character to a subscript.
	% To avoid this, replace underlines by spaces.
	caption = strrep(caption, '_', ' ');
	title(caption, 'FontSize', fontSize);
	drawnow;  % Force it to update, otherwise it waits until after the conversion to double.
	
	% Extract the individual red, green, and blue color channels.
	redChannel = rgbImage(:, :, 1);
	greenChannel = rgbImage(:, :, 2);
	blueChannel = rgbImage(:, :, 3);
	
	% Red image:
	subplot(3, 3, 4);
	imshow(redChannel, []); % Display the image.
	% Compute mean
	redChannel_Mean(k) = mean(redChannel(:));
	caption = sprintf('Green Percentage = %6.2f', Green_percentage);
	title(caption, 'FontSize', fontSize);
	% Compute and display the histogram for the Red image.
	pixelCountRed = PlotHistogramOfOneColorChannel(redChannel, 7, 'Histogram of Red Image', 'r');
	drawnow;  % Force it to update, otherwise it waits until after the conversion to double.

	% Green image:
	subplot(3, 3, 5);
	imshow(greenChannel, []); % Display the image.
	% Compute mean
	greenChannel_Mean(k) = mean(greenChannel(:));
	caption = sprintf('White Percentage = %6.2f', White_percentage);
	title(caption, 'FontSize', fontSize);
	% Compute and display the histogram for the Green image.
	pixelCountGreen = PlotHistogramOfOneColorChannel(greenChannel, 8, 'Histogram of Green Image', 'g');
	drawnow;  % Force it to update, otherwise it waits until after the conversion to double.

	% Blue image:
	subplot(3, 3, 6);
	imshow(blueChannel, []); % Display the image.
	numberOfImagesProcessed = numberOfImagesProcessed + 1;
	% Compute mean
	blueChannel_Mean(k) = mean(blueChannel(:));
	caption = sprintf('Eror rate = %6.4f', error_rate);
	title(caption, 'FontSize', fontSize);
	% Compute and display the histogram for the Blue image.
	pixelCountBlue = PlotHistogramOfOneColorChannel(blueChannel, 9, 'Histogram of Blue Image', 'b');
	drawnow;  % Force it to update, otherwise it waits until after the conversion to double.
	
	% Plot all three histograms on the same plot.
	subplot(3, 3, 2:3);
	lineWidth = 2;
	hold off;
	plot(pixelCountRed, 'r', 'LineWidth', lineWidth);
	hold on;
	grid on;
	plot(pixelCountGreen, 'g', 'LineWidth', lineWidth);
	plot(pixelCountBlue, 'b', 'LineWidth', lineWidth);
	title('All the Color Histograms (Superimposed)', 'FontSize', fontSize);
	% Set the x axis range manually to be 0-255.
	xlim([0 255]); 
	
    %Save results file
    imwrite(getframe(gcf).cdata, strcat('results/','results_',tif_image_name),'tif')

	% Prompt user to continue, unless they're at the last image.
	if k < numberOfImagesToProcess
		promptMessage = sprintf('Currently displaying image #%d of a possible %d:\n%s\n\nDo you want to\nContinue processing, or\nCancel processing?',...
			numberOfImagesProcessed, numberOfImagesToProcess2, baseFileName);
		button = questdlg(promptMessage, 'Continue?', 'Continue', 'Cancel', 'Continue');
		if strcmp(button, 'Cancel')
			break;
		end
	end
end

% Crop off any unassigned values:
redChannel_Mean = redChannel_Mean(1:numberOfImagesProcessed);
greenChannel_Mean = greenChannel_Mean(1:numberOfImagesProcessed);
blueChannel_Mean = blueChannel_Mean(1:numberOfImagesProcessed);

% Print to command window
fprintf(1, '                Filename,   Red Mean, Green Mean, Blue Mean\n');
for k = 1 : length(redChannel_Mean)
	baseFileName = imageFiles(k).name;
	fprintf(1, '%24s    %6.2f,   %6.2f,     %6.2f\n', ...
		baseFileName, redChannel_Mean(k), greenChannel_Mean(k), blueChannel_Mean(k));
end

if numberOfImagesProcessed == 1
	caption = sprintf('Done with the analysis!\n\nProcessed 1 image.\nCheck out the command window for the results');
    else
	caption = sprintf('Done with analysis!\n\nProcessed %d images.\nCheck out the command window for the results', numberOfImagesProcessed);
end
msgbox(caption);
cd ..
fprintf('STEP8: Results analysis done...\n')
pause(10);
close;

%%
%----------------------------------------------------------
% Plots a bar chart of the histogram of the color channel.
function pixelCount = PlotHistogramOfOneColorChannel(oneColorChannel, subplotNumber, caption, color)
	try
		% Let's get its histogram into 256 bins.
		[pixelCount, grayLevels] = imhist(oneColorChannel, 256);

		subplot(3, 3, subplotNumber);
		bar(grayLevels, pixelCount, 'FaceColor', color, 'BarWidth', 1); 
		title(caption, 'FontSize', 16);
		grid on;
		% Set the x axis range manually to be 0-255.
		xlim([0 255]); 
	catch ME
		errorMessage = sprintf('Error in function PlotHistogramOfOneColorChannel.\n.\n\nError Message:\n%s', ME.message);
		uiwait(warndlg(errorMessage));
	end
	return;
end

%%
function img2 = addborder(img1, t, c, stroke)
% ADDBORDER draws a border around an image
%
%    NEWIMG = ADDBORDER(IMG, T, C, S) adds a border to the image IMG with
%    thickness T, in pixels. C specifies the color, and should be in the
%    same format as the image itself. STROKE is a string indicating the
%    position of the border:
%       'inner'  - border is added to the inside of the image. The dimensions
%                  of OUT will be the same as IMG.
%       'outer'  - the border sits completely outside of the image, and does
%                  not obscure any portion of it.
%       'center' - the border straddles the edges of the image.
%

    % Input data validation
    if nargin < 4
        error('MATLAB:addborder','ADDBORDER requires four inputs.');
    end
    
    if numel(c) ~= size(img1,3)
        error('MATLAB:addborder','C does not match the color format of IMG.');
    end

    % Ensure that C has the same datatype as the image itself.
    % Also, ensure that C is a "depth" vector, rather than a row or column
    % vector, in the event that C is a 3 element RGB vector.
    c = cast(c, class(img1));
    c = reshape(c,1,1,numel(c));

    % Compute the pixel size of the new image, and allocate a matrix for it.
    switch lower(stroke(1))
        case 'i'
            img2 = addInnerStroke(img1, t, c);
        case 'o'
            img2 = addOuterStroke(img1, t, c);
        case 'c'
            img2 = addCenterStroke(img1, t, c);
        otherwise
            error('MATLAB:addborder','Invalid value for ''stroke''.');
    end
    
    
end


% Helper functions for each stroke type
function img2 = addInnerStroke(img1, t, c)

    [nr1 nc1 d] = size(img1);

    % Initially create a copy of IMG1
    img2 = img1;
    
    % Now fill the border elements of IMG2 with color C
    img2(:,1:t,:)           = repmat(c,[nr1 t 1]);
    img2(:,(nc1-t+1):nc1,:) = repmat(c,[nr1 t 1]);
    img2(1:t,:,:)           = repmat(c,[t nc1 1]);
    img2((nr1-t+1):nr1,:,:) = repmat(c,[t nc1 1]);

end

function img2 = addOuterStroke(img1, t, c)

    [nr1 nc1 d] = size(img1);

    % Add the border thicknesses to the total image size
    nr2 = nr1 + 2*t;
    nc2 = nc1 + 2*t;
    
    % Create an empty matrix, filled with the border color.
    img2 = repmat(c, [nr2 nc2 1]);
    
    % Copy IMG1 to the inner portion of the image.
    img2( (t+1):(nr2-t), (t+1):(nc2-t), : ) = img1;

end

function img2 = addCenterStroke(img1, t, c)

    % Add an inner and outer stroke of width T/2
    img2 = addInnerStroke(img1, floor(t/2), c);
    img2 = addOuterStroke(img2, ceil(t/2), c);    

end