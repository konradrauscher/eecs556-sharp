
%% Config
algorithms = {'RR'};
printEPS=false;%?
flag_cut_bounds = false;
dim_cut = 0;
thvalues = false;
L = 2^11;
ratio = 1;

NumIndexes = 3;
MatrixResults = zeros(numel(algorithms),NumIndexes);
alg = 0;

flagQNR = 0; % 1: QNR otherwise HQNR 

addpath([pwd,'/Tools']);

sensor = 'WV3';

%must be divisible by 4
AOI_SZ_PAN = 512;

%% Read specific data

I_MS = imread('C:\Users\Konrad\Documents\EECS556\proj\Rio_System-Ready_Stereo_8_Band_Bundle_50cm\056078906040\056078906040_01_P001_MUL\16JAN23133046-M1BS-056078906040_01_P001.TIF');
I_PAN = imread('C:\Users\Konrad\Documents\EECS556\proj\Rio_System-Ready_Stereo_8_Band_Bundle_50cm\056078906040\056078906040_01_P001_PAN\16JAN23133046-P1BS-056078906040_01_P001.TIF');
I_MS_LR = I_MS;

%Only process corner of image
AOI_SZ_MS = AOI_SZ_PAN/4;
I_MS_LR = I_MS_LR(1:AOI_SZ_MS,1:AOI_SZ_MS,:);
I_PAN = I_PAN(1:AOI_SZ_PAN,1:AOI_SZ_PAN);

%change dynamic range from 2^11 to 2^8
I_MS_LR = double(I_MS_LR);% / 2^3;
I_PAN = double(I_PAN);% / 2^3;
%% Show MS

if size(I_MS,3) == 4   
    showImage4LR(I_MS_LR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);    
else
    showImage8LR(I_MS_LR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);
end

%% Show PAN

%I_PAN_DS = I_PAN(1:2:end,1:2:end);

showPan(I_PAN,printEPS,2,flag_cut_bounds,dim_cut);

   

%% Actually Run RR

DELTA = 5;
addpath('../S2Sharp/irt/utilities');
addpath('../S2Sharp/irt/penalty');
potobj = potential_fun('cauchy', DELTA);
psidiff = @(x)(potobj.dpot(x));

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yim = cell(size(I_MS_LR,3),1);
Yim{1} = I_PAN; 
for ii = 2 : size(I_MS_LR,3) + 1
        Yim{ii} = I_MS_LR(:,:,ii-1);
end

mtf = 0.3 .* ones(1,size(I_MS_LR,3)+1);
d = [1, 4.*ones(1,size(I_MS_LR,3))]';

lambda = 1;
CDiter = 3;
tolgradnorm = 1e-3;

switch sensor
    case 'IKONOS'
        q = [0.000024842011904,0.000000031167283,0.000001472018618,2.912263429010321, 3.121296536058754];
    case {'GeoEye1','WV4'}
        q = [0.001308909327602,0.000009717342059,0.006023213549356, 0.047936473115674];
    case 'WV3'
        q = [0.000000603425034,0.002245365227623,0.000123865525036, 4.160835851678272,1.105372103739908];
end

r = numel(q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd RR
t2 = tic;
I_RR = RRpansharp(Yim,'r',r,'lambda',lambda,'q',q','CDiter',CDiter,'tolgradnorm',tolgradnorm,'d',d,'mtf',mtf,'W_method','sobel');%,'dpot',psidiff);
time_RR = toc(t2);
fprintf('Elaboration time RR: %.2f [sec]\n',time_RR);
cd ..

%% evaluation
[D_lambda_RR,D_S_RR,QNRI_RR] = indexes_evaluation_FS(I_RR,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

%MatrixResults(alg,:) = [D_lambda_RR,D_S_RR,QNRI_RR];
%MatrixImage(:,:,:,alg) = I_RR;


%% Show Sharpened


if size(I_MS,3) == 4   
    showImage4LR(I_RR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);    
else
    showImage8LR(I_RR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);
end