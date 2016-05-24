%% Local thresholding method

function thresh = LocalThresh(im, ROI_SD)

im = double(im);

mean_im = mean(im(:));
SD_im = std(im(:));

%% Local thresholding - Niblack
k = 0.5;
thresh = mean_im+k*SD_im;


% %% Sauvola
% R = 28;
% k = 0.5;
% thresh = mean_im * (1+k*((SD_im/R)-1));