% An example illustrating how to use S2sharp using
% Sentinel-2 dataset simulated from Aviris.  For detail see 
% Sentinel-2 Sharpening Using a Reduced-Rank Method, M. Ulfarsson et al., 
% IEEE Transactions on Geoscience and Remote Sensing, 2019, 
% 
% https://www.researchgate.net/publication/332656193_Sentinel-2_Sharpening_Using_a_Reduced-Rank_Method

load Data/Aviris_cell_3.mat; 
r = 8; % subspace dimension / the rank
q  = [1, 0.3851, 6.9039, 19.9581, 47.8967, 27.5518, 2.7100, 34.8689]; % lam_i = lam*qi
ni = 5;  % number of iterations / one can select fewer iterations (e.g. ni=2) which yields similar result
          % but is much quicker
lam = 1.8998e-04;


%wmethods = {'original','original','prewitt','sobel'};
%dmethods = {'original','original_shift','prewitt','sobel'};
%wmethods = {'original','prewitt','sobel'};
%wmethods = {'sobel'};
%potnames = {'quad'};
%potnames = {'hyper3','huber','broken','cauchy'};
%potdeltas = {    0.5,   0.75,     0.5,     0.5};
%wmethods = {'nothing', 'original', 'sobel'};
%potnames = {'quad','cauchy'};
%potdeltas = {  nan,   0.5};
wmethods = {'sobel'};
potnames = {'quad'};
potdeltas = {nan};

%hopefully setup MIRT
addpath('irt/utilities');
addpath('irt/penalty');

%psidiff = @(x) m_huber_dpot(x,DELTA);

%dmethods = {'original', 'original', 'original'};
for ii = 1:length(wmethods)
    for jj = 1:length(potnames)
        w_method = wmethods{ii};
        potname = potnames{jj};
        delta = potdeltas{jj};

        potobj = potential_fun(potname, delta);
        psidiff = @(x)(potobj.dpot(x));
        [ Xhat_im, output_S2 ]=S2sharp(Yim,'Xm_im',Xm_im,'r',r,'lambda',lam,'q',q, ...
            'CDiter',ni,'W_method',w_method,'dpot',psidiff);

        % Output
        S2sharp_SRE = output_S2.SRE{end}([1,5,6,7,9:12]);
        S2sharp_SAM = output_S2.SAMm(end);
        S2sharp_aSRE= mean(S2sharp_SRE );    
        S2sharp_RMSE = output_S2.RMSE(end);
        S2sharp_NRMSE = output_S2.NRMSE(end);
        S2sharp_aSSIM = output_S2.aSSIM(end);
        S2sharp_ERGAS_60m=output_S2.ERGAS_60m(end);
        S2sharp_ERGAS_20m=output_S2.ERGAS_20m(end);
        S2sharp_time = output_S2.Time;

        disp(['WEIGHT METHOD = ' w_method]);
        disp(['POTENTIAL FUNCTION = ' potname]);
        disp(['POTENTIAL FUNCTION DELTA = ' num2str(delta)]);
        disp(['S2sharp: Best lambda = ' num2str(lam(end))]);
        disp(['S2sharp: SAM = ' num2str(S2sharp_SAM)])
        disp(['Average SRE = ' num2str(S2sharp_aSRE)]);
        disp(['S2sharp aSSIM = ' num2str(S2sharp_aSSIM)]);
        disp(['S2sharp RMSE = ' num2str(S2sharp_RMSE)]);
        disp(['S2sharp NRMSE = ' num2str(S2sharp_NRMSE)]);
        disp(['S2sharp time = ' num2str(S2sharp_time)]);
        disp(['S2sharp: SRE:'])
        disp(['B1    B5    B6    B7    B8a   B9    B11   B12'])
        fprintf('%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n\n', S2sharp_SRE)
    end
end

%% Display low-res/sharpened version of band 
BANDNUM = 1;%1 and 10 seem to give best results
limsub = 2;
Xm_disp = Xm_im(limsub+1:end-limsub,limsub+1:end-limsub,:);

B1_lowres = Yim{BANDNUM};
B1_truth = Xm_disp(:,:,BANDNUM);
B1_sharpened = Xhat_im(:,:,BANDNUM);

clim = [min(B1_truth(:)),max(B1_truth(:))];

B1_lowres_upsampled = imresize(B1_lowres,size(Xm_im(:,:,1)));
B1_lowres_upsampled = B1_lowres_upsampled(limsub+1:end-limsub,limsub+1:end-limsub,:);

axes = [];
figure; imagesc(B1_lowres_upsampled); colormap gray; caxis(clim); colorbar; title('Low resolution band'); axes(end+1) = gca;
figure; imagesc(B1_sharpened); colormap gray; caxis(clim);  colorbar; title('Sharpened band'); axes(end+1) = gca;
figure; imagesc(B1_truth); colormap gray; caxis(clim);  colorbar; title('Ground truth'); axes(end+1) = gca;
linkaxes(axes);
