function [UV_area_mm, im, disp_im, image, UV_roi_mask, smoothed_flag, edge_UV, overlap_thresh, UV_region] = Analyze_CUVAF(handles)

pathname = handles.pathname;
checkbox = get(handles.ROI_checkbox, 'Value');

foldername = [];                    %ie. C:\Users\emily huynh\Desktop\160814\10331\
path_parts = regexp(pathname, handles.delimiter, 'split');
num_parts = length(path_parts);
for i = 1: num_parts-1
    foldername = [foldername path_parts{i} handles.delimiter];
end
imagename = path_parts{num_parts};

path_parts = regexp(imagename, '\.', 'split');
imagename = path_parts{1};          %ie. OS temporal UV_AF

resize_factor = handles.resize_factor;
calibfactor = handles.calibfactor;

%% Read Image
im=imread(pathname);
% figure; imshow(im);

im = imresize(im, resize_factor);
[m,n,o] = size(im);

imG = im(:,:,2);
imB = im(:,:,3);

meanG = mean(double(imG(:)));
sdG = std(double(imG(:)));
meanB = mean(double(imB(:)));
sdB = std(double(imB(:)));

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
% figure; imshow(image);

%% image stats 
max_int_pp = max(image(:));
min_int_pp = min(image(:));

im_mean_pp = mean(double(image(:)));
im_sd_pp = std(double(image(:)));

image_stats = [im_mean_pp im_sd_pp min_int_pp max_int_pp meanG sdG meanB sdB];

%% UV region mask
if checkbox == 0 
    temp = figure;
    [UV_roi_mask, xi, yi] = roipoly(disp_im);    
    title('SELECT UV ROI');
    try
        close(temp);
    catch
        set(handles.measure_pushbutton,'Enable','on');
        UV_area_mm = 'NA';
        im = 'NA';
        disp_im = 'NA';
        image = 'NA';
        UV_roi_mask = 'NA';
        smoothed_flag = 'NA';
        edge_UV = 'NA'; 
        overlap_thresh = 'NA';
        UV_region = 'NA';
        return;
    end
    % make sure user selects ROI
    while isempty(xi)
        temp = figure;
        [UV_roi_mask, xi, yi] = roipoly(disp_im);     
        title('SELECT UV ROI');
        try 
            close(temp)
        catch
            set(handles.measure_pushbutton,'Enable','on');
            UV_area_mm = 'NA';
            im = 'NA';
            disp_im = 'NA';
            image = 'NA';
            UV_roi_mask = 'NA';
            smoothed_flag = 'NA';
            edge_UV = 'NA'; 
            overlap_thresh = 'NA';
            UV_region = 'NA';
            return;
        end
    end    
    ROI_coords = [xi yi];       % xi = R, yi = C 
   
elseif checkbox == 1
    % Use saved ROI coordinates
    filename = handles.excel_filename;
    path = [foldername filename];
    if exist(path, 'file')
        if isempty(regexpi(imagename, 'OD NASAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'A2:B60');                           
        elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'D2:E60');
        elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'G2:H60');
        elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'J2:K60');
        else
            msgbox('No ROI coordinates stored for this image')
            return
        end
    else
        msgbox('No data for this image')
        return
    end  
    if isempty(ROI_coords) == 1 
        msgbox('No ROI coordinates stored for this image')
        return
    end
end

UV_roi_mask = roipoly(disp_im, ROI_coords(:,1), ROI_coords(:,2)); 

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
% figure; imshow(UV_roi_mask);

%% Apply mask to grayscale image
uv_image_roi = double(UV_roi_mask2).*double(image);
uv_image_roi = uint8(uv_image_roi);
% figure; imshow(uv_image_roi);
% title('User-selected ROI');

%% Local thresholding to find UV region
win_ratio = handles.win_ratio;
win_width = round((sqrt(win_ratio*roi_area))); %200;
step = round(win_width/3);

UV_near_edge = 0;
if minR_UV-step < 0 || maxR_UV+step > m || minC_UV-step < 0 || maxC_UV+step > n
    UV_near_edge = 1;
    offset = step*2;  
    pad_im = padarray(uv_image_roi, [step*2, step*2], 0);
else
    offset = 0;
    pad_im = uv_image_roi;
end

[pm, pn] = size(pad_im);

uv_flag = zeros(pm,pn);
dilate_mask = bwmorph(UV_roi_mask, 'dilate', step);
dilate_mask = dilate_mask+UV_roi_mask;

for i = minR_UV-step:step:maxR_UV+step
    for j = minC_UV-step:step:maxC_UV+step
        if i-win_width+offset > 0 && i+win_width+offset <= pm && j-win_width+offset > 0 && j+win_width+offset <= pn && dilate_mask(i,j) > 0
            if dilate_mask(i,j) == 2 
                temp_window = pad_im(i-win_width+offset:i+win_width+offset, j-win_width+offset:j+win_width+offset);
%                 [thresh, BCV, meanB, meanF] = Otsu(temp_window, 1);
                thresh = LocalThresh(temp_window, im_sd_pp);
                thresh_im = temp_window>thresh;
                full_thresh_im = zeros(pm,pn);
                full_thresh_im(i-win_width+offset:i+win_width+offset, j-win_width+offset:j+win_width+offset) = thresh_im;      
                uv_flag = uv_flag+full_thresh_im;
            elseif dilate_mask(i,j) == 1
                win_width = round(win_width/2);
                temp_window = pad_im(i-win_width+offset:i+win_width+offset, j-win_width+offset:j+win_width+offset);
%                 [thresh, BCV, meanB, meanF] = Otsu(temp_window, 1);
                thresh = LocalThresh(temp_window, im_sd_pp);
                thresh_im = temp_window>thresh;
                full_thresh_im = zeros(pm,pn);
                full_thresh_im(i-win_width+offset:i+win_width+offset, j-win_width+offset:j+win_width+offset) = thresh_im;    
                uv_flag = uv_flag+full_thresh_im;
                win_width = round((sqrt(win_ratio*roi_area)));
            end
        end
    end
end

%remove padding
if UV_near_edge == 1
    uv_flag = uv_flag(offset+1:pm-offset, offset+1:pn-offset);
end

gaussfilt_size = round(step);
gaussian = fspecial('average', [gaussfilt_size gaussfilt_size]);
smoothed_flag = imfilter(double(uv_flag), gaussian); 

smoothed_flag = smoothed_flag.*double(dilate_mask);

overlap_thresh = 1; 
uv_flag2 = smoothed_flag >= overlap_thresh;
uv_flag2 = double(uv_flag2).*smoothed_flag;

UV_region = UV_roi_mask.*logical(uv_flag2);       
UV_region = imfill(UV_region, 'holes');
area_thresh = round(resize_factor*(handles.area_filt/(calibfactor^2)));
UV_region = bwareaopen(UV_region, area_thresh);

%% Delineate UV region
edge_UV = bwperim(UV_region);

%% Calculate UVAF area
UV_pix_area = sum(UV_region(:));
UV_area_mm = UV_pix_area*(calibfactor^2);

%% create excel spreadsheet - save ROI coordinates
if checkbox == 0 
    filename = handles.excel_filename;
    path = [foldername filename];
    clear_cells = nan(60,2);

    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, clear_cells, 'ROI Coordinates', 'A2');
        xlwrite(path, ROI_coords, 'ROI Coordinates', 'A2');
        xlwrite(path, roi_area, 'ROI Coordinates', 'B1');
        xlwrite(path, {datestr(now)}, 'ROI Coordinates', 'C1');
        xlwrite(path, image_stats, 'Results', 'N2');

    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0     
        xlwrite(path, clear_cells, 'ROI Coordinates', 'D2');
        xlwrite(path, ROI_coords, 'ROI Coordinates', 'D2');
        xlwrite(path, roi_area, 'ROI Coordinates', 'E1');
        xlwrite(path, {datestr(now)}, 'ROI Coordinates', 'F1');
        xlwrite(path, image_stats, 'Results', 'N3');

    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, clear_cells, 'ROI Coordinates', 'G2');
        xlwrite(path, ROI_coords, 'ROI Coordinates', 'G2');
        xlwrite(path, roi_area, 'ROI Coordinates', 'H1');
        xlwrite(path, {datestr(now)}, 'ROI Coordinates', 'I1');
        xlwrite(path, image_stats, 'Results', 'N4');

    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
        xlwrite(path, clear_cells, 'ROI Coordinates', 'J2');
        xlwrite(path, ROI_coords, 'ROI Coordinates', 'J2');
        xlwrite(path, roi_area, 'ROI Coordinates', 'K1');
        xlwrite(path, {datestr(now)}, 'ROI Coordinates', 'L1');
        xlwrite(path, image_stats, 'Results', 'N5');
    else
        msgbox('Error writing ROI coordinates to excel. Check spelling of image name.', 'warn');
        return;  
    end
end


