function [UV_area_mm, edge_UV, UV_region] = Refine_UVAF(handles)

overlap_thresh = handles.overlap_thresh;
smoothed_flag = handles.smoothed_flag;
UV_roi_mask = handles.UV_roi_mask;
resize_factor = handles.resize_factor;
calibfactor = handles.calibfactor;

uv_flag2 = smoothed_flag >= overlap_thresh;

uv_flag2 = double(uv_flag2).*smoothed_flag;

UV_region = logical(UV_roi_mask).*logical(uv_flag2);       
UV_region = imfill(UV_region, 'holes');
area_thresh = round(resize_factor*(handles.area_filt/(calibfactor^2)));
UV_region = bwareaopen(UV_region, area_thresh);

%% Delineate UV region
edge_UV = bwperim(UV_region);

%% Calculate UVAF area
UV_pix_area = sum(UV_region(:));
UV_area_mm = UV_pix_area*(calibfactor^2);


% UV_pix_area