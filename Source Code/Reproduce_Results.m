function [UV_area_mm, edge_UV, im, overlap_thresh, subjective, user_initials, comments] = Reproduce_Results(handles)

% use the browse button to specify the pathname and then click reproduce results
% pathname = C:\Users\emily huynh\Desktop\160814\10331\OS temporal UV_AF.JPG

pathname = handles.pathname;

foldername = [];                        % ie. C:\Users\emily huynh\Desktop\160814\10331\
path_parts = regexp(pathname, handles.delimiter, 'split');
num_parts = length(path_parts);
for i = 1: num_parts-1
    foldername = [foldername path_parts{i} handles.delimiter];
end
imagename = path_parts{num_parts};

path_parts = regexp(imagename, '\.', 'split');
imagename = path_parts{1};              % ie OS temporal UV_AF

filename = handles.excel_filename;
path = [foldername filename];

%% Find image
if exist(path, 'file') ~= 0 
    win_ratio = xlsread(path, 1, 'B13');
    area_filt = xlsread(path, 1, 'B14');
    
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        ROI_coords = xlsread(path, 'ROI Coordinates', 'A2:B60');
        resize_factor = xlsread(path, 'Results', 'B2');
        overlap_thresh = xlsread(path, 'Results', 'C2');
        area = xlsread(path, 'Results', 'D2');
        [num notes data] = xlsread(path, 'Results', 'E2');
        remove_points = xlsread(path, 'Removed Regions', 'A2:B20');
        [num subjective data] = xlsread(path, 'Results', 'H2');
        [num user_initials data] = xlsread(path, 'Results', 'I2'); 
        [num comments data] = xlsread(path, 'Results', 'J2'); 

    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        ROI_coords = xlsread(path, 'ROI Coordinates', 'D2:E60');
        resize_factor = xlsread(path, 'Results', 'B3');
        overlap_thresh = xlsread(path, 'Results', 'C3');
        area = xlsread(path, 'Results', 'D3');
        [num notes data] = xlsread(path, 'Results', 'E3');
        remove_points = xlsread(path, 'Removed Regions', 'C2:D20');
        [num subjective data] = xlsread(path, 'Results', 'H3');
        [num user_initials data] = xlsread(path, 'Results', 'I3'); 
        [num comments data] = xlsread(path, 'Results', 'J3'); 

    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        ROI_coords = xlsread(path, 'ROI Coordinates', 'G2:H60');
        resize_factor = xlsread(path, 'Results', 'B4');
        overlap_thresh = xlsread(path, 'Results', 'C4');
        area = xlsread(path, 'Results', 'D4');
        [num notes data] = xlsread(path, 'Results', 'E4');
        remove_points = xlsread(path, 'Removed Regions', 'E2:F20');
        [num subjective data] = xlsread(path, 'Results', 'H4');
        [num user_initials data] = xlsread(path, 'Results', 'I4'); 
        [num comments data] = xlsread(path, 'Results', 'J4'); 

    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
        ROI_coords = xlsread(path, 'ROI Coordinates', 'J2:K60');
        resize_factor = xlsread(path, 'Results', 'B5');
        overlap_thresh = xlsread(path, 'Results', 'C5');
        area = xlsread(path, 'Results', 'D5');
        [num notes data] = xlsread(path, 'Results', 'E5');
        remove_points = xlsread(path, 'Removed Regions', 'G2:H20');
        [num subjective data] = xlsread(path, 'Results', 'H5');
        [num user_initials data] = xlsread(path, 'Results', 'I5'); 
        [num comments data] = xlsread(path, 'Results', 'J5'); 
    end

    if strcmp(notes,'No UVAF') == 1
        msgbox('No UVAF');
        UV_area_mm = area;
        edge_UV = [];
        im=imread(pathname);
        return;
    end
    if isempty(ROI_coords)
        msgbox('ROI has not been defined'); 
        return;
    end
else
    msgbox('CUVAF has not been measured!!!');
    return
end
    
calibfactor = handles.calibfactor;

%% Read Image
im=imread(pathname);

im = imresize(im, resize_factor);
[m,n,o] = size(im);
imG = im(:,:,2);
imB = im(:,:,3);
meanG = mean2(double(imG));

gamma = (1/50)*meanG + 0.6;              % min gamma = y-intercept
if gamma > 1.5
    gamma = 1.5;
end

tol = stretchlim(imB, [0.02 0.99]);
im2 = imadjust(im, tol, [0;1], gamma);
imG2 = im2(:,:,2);
imB2 = im2(:,:,3);

disp_im = cat(3, zeros(m,n), imG2, imB2);
image = rgb2gray(disp_im);

%% image stats 
% max_int_tophat = max(image(:));
% min_int_tophat = min(image(:));
% 
% im_mean_tophat = mean(double(image(:)));
im_sd_tophat = std(double(image(:)));

%% mask
UV_roi_mask = roipoly(image, ROI_coords(:,1), ROI_coords(:,2));
ROI_info = regionprops(UV_roi_mask, 'PixelList');
pix_list = [ROI_info.PixelList];
minR_UV = min(pix_list(:,2));
maxR_UV = max(pix_list(:,2));
minC_UV = min(pix_list(:,1));
maxC_UV = max(pix_list(:,1));

roi_area = sum(UV_roi_mask(:));    

inv_mask = imcomplement(UV_roi_mask);
inv_mask = 0.9*inv_mask;
UV_roi_mask2 = inv_mask+UV_roi_mask;

%% Apply mask to grayscale image
uv_image_roi = double(UV_roi_mask2).*double(image);
uv_image_roi = uint8(uv_image_roi);

%% Local thresholding to find UV region
win_width = round((sqrt(win_ratio*roi_area))); %200;
step = round(win_width/3);
uv_flag = zeros(m,n);
dilate_mask = bwmorph(UV_roi_mask, 'dilate', step);
dilate_mask = dilate_mask+UV_roi_mask;

for i = minR_UV-step:step:maxR_UV+step
    for j = minC_UV-step:step:maxC_UV+step
        if i-win_width > 0 && i+win_width <= m && j-win_width > 0 && j+win_width <= n && dilate_mask(i,j) > 0
            if dilate_mask(i,j) == 2 
                temp_window = uv_image_roi(i-win_width:i+win_width, j-win_width:j+win_width);
%                 [thresh, BCV, meanB, meanF] = Otsu(temp_window, 1);
                thresh = LocalThresh(temp_window, im_sd_tophat);
                thresh_im = temp_window>thresh;
                full_thresh_im = zeros(m,n);
                full_thresh_im(i-win_width:i+win_width, j-win_width:j+win_width) = thresh_im;      
                uv_flag = uv_flag+full_thresh_im;
            elseif dilate_mask(i,j) == 1
                win_width = round(win_width/2);
                temp_window = uv_image_roi(i-win_width:i+win_width, j-win_width:j+win_width);
%                 [thresh, BCV, meanB, meanF] = Otsu(temp_window, 1);
                thresh = LocalThresh(temp_window, im_sd_tophat);
                thresh_im = temp_window>thresh;
                full_thresh_im = zeros(m,n);
                full_thresh_im(i-win_width:i+win_width, j-win_width:j+win_width) = thresh_im;    
                uv_flag = uv_flag+full_thresh_im;
                win_width = round((sqrt(win_ratio*roi_area)));
            end
        end
    end
end

% figure; imshow(double(uv_flag), []); colormap jet; colorbar
% title('UV flag overlap values')

gaussfilt_size = step;
gaussian = fspecial('average', [gaussfilt_size gaussfilt_size]);
smoothed_flag = imfilter(double(uv_flag), gaussian); 
smoothed_flag = smoothed_flag.*double(dilate_mask);
% figure; imshow(smoothed_flag, []); colormap jet; colorbar
% title('smoothed overlap')

uv_flag2 = smoothed_flag >= overlap_thresh;
uv_flag2 = double(uv_flag2).*smoothed_flag;

UV_region = UV_roi_mask.*logical(uv_flag2);       
UV_region = imfill(UV_region, 'holes');

area_thresh = round(resize_factor*(area_filt/(calibfactor^2)));
UV_region = bwareaopen(UV_region, area_thresh);
% figure; imshow(UV_region);
% title('UV_region > 100');


%% remove regions
[imL,num_regs] = bwlabeln(UV_region);
if ~isempty(remove_points)
    for i = 1:size(remove_points,1)
        row = round(remove_points(i,1));
        col = round(remove_points(i,2));
        label = imL(row,col);
        temp = imL == label;            % binary image of region to remove
        imL(temp) = 0;
    end
end

filt_UV_region = logical(imL);

%% Delineate UV region
edge_UV = bwperim(filt_UV_region);

%% Calculate UVAF area
UV_pix_area = sum(filt_UV_region(:));
UV_area_mm = UV_pix_area*(calibfactor^2);


