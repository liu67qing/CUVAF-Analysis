function [UV_area_mm, edge_UV, new_UV_region] = RemoveRegions(handles)

UV_region = handles.UV_region;
selected_points = handles.selected_points;
calibfactor = handles.calibfactor;

%%
[imL,num_regs] = bwlabeln(UV_region);
regionAvailable = 1:num_regs;

% Get label of regions to remove
removeRegionLabel = [];
for i = 1:size(selected_points,1)
    row = round(selected_points(i,1));
    col = round(selected_points(i,2));
    label = imL(row,col);
    if label ~= 0 
        removeRegionLabel = [removeRegionLabel label];
    end
end

% Remove Regions
newAvailableRegion=[];
for i=1:length(regionAvailable)
    if sum(regionAvailable(i)==removeRegionLabel)==0         
        newAvailableRegion=[newAvailableRegion;regionAvailable(i)];
    end
end
        
new_UV_region = ismember(imL, newAvailableRegion);

%% Delineate UV region
edge_UV = bwperim(new_UV_region);

%% Calculate UVAF area
UV_pix_area = sum(new_UV_region(:));
UV_area_mm = UV_pix_area*(calibfactor^2);

