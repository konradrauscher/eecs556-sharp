%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MAIN: REDUCED RESOLUTION VALIDATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%% Analyzed image choice

% sensor = 'WV3';
% im_tag = 'NY1';

sensor = 'GeoEye1';
im_tag = 'GeoEye1_July';

% sensor = 'IKONOS';
% im_tag = 'Toulouse';

%% Tools
addpath([pwd,'/Tools']);
addpath([pwd,'/Datasets']);

%% Quality Index Blocks
Qblocks_size = 32;

%% Interpolator
bicubic = 0;

%% Cut Final Image
flag_cut_bounds = 1;
dim_cut = 21;

%% Threshold values out of dynamic range
thvalues = 0;

%% Print Eps
printEPS = 0;

%% Resize Factor
ratio = 4;

%% Radiometric Resolution
L = 11;

%% %%%%%%%%%%%%%%%%%%%%%%%% Dataset load %%%%%%%%%%%%%%%%%%%%%%%%%%
switch im_tag
    case 'NY1'
        load('NY1_WV3_FR.mat','I_MS_LR','I_PAN');
        I_MS_LR = double(I_MS_LR);
        I_PAN = double(I_PAN);
    case 'GeoEye1_July'
        load('Collazzone_GeoEye_July_FR.mat','I_MS_LR','I_PAN');
        I_MS_LR = double(I_MS_LR);
        I_PAN = double(I_PAN);
    case 'Toulouse'
        load('Toulouse_IKONOS_FR.mat','I_MS_LR','I_PAN');
        I_MS_LR = double(I_MS_LR);
        I_PAN = double(I_PAN);
end

%% GT
I_GT = I_MS_LR;

%% %%%%%%%%%%%%%    Preparation of image to fuse            %%%%%%%%%%%%%%
[I_MS_LR, I_PAN] = resize_images(I_MS_LR,I_PAN,ratio,sensor);

%% Upsampling
if bicubic == 1
    H = zeros(size(I_PAN,1),size(I_PAN,2),size(I_MS_LR,3));    
    for idim = 1 : size(I_MS_LR,3)
        H(:,:,idim) = imresize(I_MS_LR(:,:,idim),ratio);
    end
    I_MS = H;
else
    I_MS = interp23tap(I_MS_LR,ratio);
end