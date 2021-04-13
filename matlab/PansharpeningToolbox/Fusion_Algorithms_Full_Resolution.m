%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  RUN AND FULL RESOLUTION VALIDATION  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select algorithms to run

algorithms = {'EXP','BT-H','BDSD','C-BDSD','BDSD-PC','GS','GSA','C-GSA','PRACS','AWLP',...
         'MTF-GLP','MTF-GLP-FS','MTF-GLP-HPM','MTF-GLP-HPM-H','MTF-GLP-HPM-R','MTF-GLP-CBD','C-MTF-GLP-CBD','MF','FE-HPM','SR-D','PWMBF','TV','RR','PNN','PNN-IDX','A-PNN','A-PNN-FT'};

%% Initialization of the Matrix of Results
NumIndexes = 3;
MatrixResults = zeros(numel(algorithms),NumIndexes);
alg = 0;

%% Flag QNR/HQNR
flagQNR = 0; % 1: QNR otherwise HQNR 

%% Tools
addpath([pwd,'/Tools']);

%% MS

if size(I_MS,3) == 4   
    showImage4LR(I_MS_LR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);    
else
    showImage8LR(I_MS_LR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);
end

%% PAN

showPan(I_PAN,printEPS,2,flag_cut_bounds,dim_cut);

%% EXP

if ismember('EXP',algorithms)
    alg = alg + 1;
    
    [D_lambda_EXP,D_S_EXP,QNRI_EXP] = indexes_evaluation_FS(I_MS,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);
    
    MatrixResults(alg,:) = [D_lambda_EXP,D_S_EXP,QNRI_EXP];
    MatrixImage(:,:,:,alg) = I_MS;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Component Substitution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BT-H

if ismember('BT-H',algorithms)
    alg = alg + 1;

    cd BT-H
    t2=tic;
    I_BT_H = BroveyRegHazeMin(I_MS,I_PAN,ratio);
    time_BT_H = toc(t2);
    fprintf('Elaboration time BT-H: %.2f [sec]\n',time_BT_H);
    cd ..

    %%% Quality indexes computation
    [D_lambda_BT_H,D_S_BT_H,QNRI_BT_H] = indexes_evaluation_FS(I_BT_H,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_BT_H,D_S_BT_H,QNRI_BT_H];
    MatrixImage(:,:,:,alg) = I_BT_H;
end

%% BDSD

if ismember('BDSD',algorithms)
    alg = alg + 1;
    
    cd BDSD
    t2=tic;
    I_BDSD = BDSD(I_MS,I_PAN,ratio,size(I_MS,1),sensor);
    time_BDSD = toc(t2);
    fprintf('Elaboration time BDSD: %.2f [sec]\n',time_BDSD);
    cd ..

    [D_lambda_BDSD,D_S_BDSD,QNRI_BDSD] = indexes_evaluation_FS(I_BDSD,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_BDSD,D_S_BDSD,QNRI_BDSD];
    MatrixImage(:,:,:,alg) = I_BDSD;
end

%% C-BDSD

if ismember('C-BDSD',algorithms)
    alg = alg + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter setting
    numclusters = 30; % number of clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd BDSD
    t2=tic;
    I_C_BDSD = C_BDSD(I_MS,I_PAN,ratio,sensor,numclusters);
    time_C_BDSD = toc(t2);
    fprintf('Elaboration time C-BDSD: %.2f [sec]\n',time_C_BDSD);
    cd ..

    [D_lambda_C_BDSD,D_S_C_BDSD,QNRI_C_BDSD] = indexes_evaluation_FS(I_C_BDSD,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_C_BDSD,D_S_C_BDSD,QNRI_C_BDSD];
    MatrixImage(:,:,:,alg) = I_C_BDSD;
end

%% BDSD-PC

if ismember('BDSD-PC',algorithms)
    alg = alg + 1;
    
    cd BDSD
    t2=tic;
    I_BDSD_PC = BDSD_PC(I_MS,I_PAN,ratio,sensor);
    time_BDSD_PC = toc(t2);
    fprintf('Elaboration time BDSD-PC: %.2f [sec]\n',time_BDSD_PC);
    cd ..

    [D_lambda_BDSD_PC,D_S_BDSD_PC,QNRI_BDSD_PC] = indexes_evaluation_FS(I_BDSD_PC,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_BDSD_PC,D_S_BDSD_PC,QNRI_BDSD_PC];
    MatrixImage(:,:,:,alg) = I_BDSD_PC;
end

%% GS

if ismember('GS',algorithms)
    alg = alg + 1;
    
    cd GS
    t2=tic;
    I_GS = GS(I_MS,I_PAN);
    time_GS = toc(t2);
    fprintf('Elaboration time GS: %.2f [sec]\n',time_GS);
    cd ..

    [D_lambda_GS,D_S_GS,QNRI_GS] = indexes_evaluation_FS(I_GS,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_GS,D_S_GS,QNRI_GS];
    MatrixImage(:,:,:,alg) = I_GS;
end

%% GSA

if ismember('GSA',algorithms)
    alg = alg + 1;

    cd GS
    t2=tic;
    I_GSA = GSA(I_MS,I_PAN,I_MS_LR,ratio);
    time_GSA = toc(t2);
    fprintf('Elaboration time GSA: %.2f [sec]\n',time_GSA);
    cd ..

    [D_lambda_GSA,D_S_GSA,QNRI_GSA] = indexes_evaluation_FS(I_GSA,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_GSA,D_S_GSA,QNRI_GSA];
    MatrixImage(:,:,:,alg) = I_GSA;
end

%% C-GSA

if ismember('C-GSA',algorithms)
    alg = alg + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PS_algorithm = 'GSA'; % Pansharpening algorithm 
    n_segm = 5; % Number of segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd GS

    %%% Pansharpening Algorithm
    t2=tic;
    I_C_GSA = GS_Segm(I_MS,I_PAN,gen_LP_image(PS_algorithm,I_MS,I_PAN,I_MS_LR,ratio,sensor), k_means_clustering(I_MS,n_segm));
    time_C_GSA = toc(t2);
    fprintf('Elaboration time C-GSA: %.2f [sec]\n',time_C_GSA);
    cd ..

    %%% Quality indexes computation
    [D_lambda_C_GSA,D_S_C_GSA,QNRI_C_GSA] = indexes_evaluation_FS(I_C_GSA,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_C_GSA,D_S_C_GSA,QNRI_C_GSA];
    MatrixImage(:,:,:,alg) = I_C_GSA;
end

%% PRACS

if ismember('PRACS',algorithms)
    alg = alg + 1;
    
    cd PRACS
    t2=tic;
    I_PRACS = PRACS(I_MS,I_PAN,ratio);
    time_PRACS = toc(t2);
    fprintf('Elaboration time PRACS: %.2f [sec]\n',time_PRACS);
    cd ..

    [D_lambda_PRACS,D_S_PRACS,QNRI_PRACS] = indexes_evaluation_FS(I_PRACS,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_PRACS,D_S_PRACS,QNRI_PRACS];
    MatrixImage(:,:,:,alg) = I_PRACS;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MultiResolution Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AWLP

if ismember('AWLP',algorithms)
    alg = alg + 1;
    
    cd AWLP
    t2=tic;
    I_AWLP = AWLP(I_MS,I_PAN,ratio);
    time_AWLP = toc(t2);
    fprintf('Elaboration time AWLP: %.2f [sec]\n',time_AWLP);
    cd ..

    [D_lambda_AWLP,D_S_AWLP,QNRI_AWLP] = indexes_evaluation_FS(I_AWLP,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_AWLP,D_S_AWLP,QNRI_AWLP];
    MatrixImage(:,:,:,alg) = I_AWLP;
end

%% MTF-GLP

if ismember('MTF-GLP',algorithms)
    alg = alg + 1;
    
    cd GLP
    t2=tic;
    I_MTF_GLP = MTF_GLP(I_MS,I_PAN,sensor,ratio);
    time_MTF_GLP = toc(t2);
    fprintf('Elaboration time MTF-GLP: %.2f [sec]\n',time_MTF_GLP);
    cd ..

    [D_lambda_MTF_GLP,D_S_MTF_GLP,QNRI_MTF_GLP] = indexes_evaluation_FS(I_MTF_GLP,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MTF_GLP,D_S_MTF_GLP,QNRI_MTF_GLP];
    MatrixImage(:,:,:,alg) = I_MTF_GLP;
end

%% MTF-GLP-FS

if ismember('MTF-GLP-FS',algorithms)
    alg = alg + 1;

    cd GLP
    t2=tic;
    I_MTF_GLP_FS = MTF_GLP_FS(I_MS,I_PAN,sensor,ratio);
    time_MTF_GLP_FS = toc(t2);
    fprintf('Elaboration time MTF-GLP-FS: %.2f [sec]\n',time_MTF_GLP_FS);
    cd ..

    [D_lambda_MTF_GLP_FS,D_S_MTF_GLP_FS,QNRI_MTF_GLP_FS] = indexes_evaluation_FS(I_MTF_GLP_FS,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MTF_GLP_FS,D_S_MTF_GLP_FS,QNRI_MTF_GLP_FS];
    MatrixImage(:,:,:,alg) = I_MTF_GLP_FS;
end

%% MTF-GLP-HPM

if ismember('MTF-GLP-HPM',algorithms)
    alg = alg + 1;
    
    cd GLP
    t2=tic;
    I_MTF_GLP_HPM = MTF_GLP_HPM(I_MS,I_PAN,sensor,ratio);
    time_MTF_GLP_HPM = toc(t2);
    fprintf('Elaboration time MTF-GLP-HPM: %.2f [sec]\n',time_MTF_GLP_HPM);
    cd ..

    [D_lambda_MTF_GLP_HPM,D_S_MTF_GLP_HPM,QNRI_MTF_GLP_HPM] = indexes_evaluation_FS(I_MTF_GLP_HPM,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MTF_GLP_HPM,D_S_MTF_GLP_HPM,QNRI_MTF_GLP_HPM];
    MatrixImage(:,:,:,alg) = I_MTF_GLP_HPM;
end

%% MTF-GLP-HPM-H

if ismember('MTF-GLP-HPM-H',algorithms)
    alg = alg + 1;
    
    cd GLP
    t2=tic;
    I_MTF_GLP_HPM_H = MTF_GLP_HPM_Haze_min(I_MS,I_PAN,sensor,ratio,1);
    time_MTF_GLP_HPM_H = toc(t2);
    fprintf('Elaboration time MTF-GLP-HPM-H: %.2f [sec]\n',time_MTF_GLP_HPM_H);
    cd ..

    %%% Quality indexes computation
    [D_lambda_MTF_GLP_HPM_H,D_S_MTF_GLP_HPM_H,QNRI_MTF_GLP_HPM_H] = indexes_evaluation_FS(I_MTF_GLP_HPM_H,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MTF_GLP_HPM_H,D_S_MTF_GLP_HPM_H,QNRI_MTF_GLP_HPM_H];
    MatrixImage(:,:,:,alg) = I_MTF_GLP_HPM_H;
end

%% MTF-GLP-HPM-R

if ismember('MTF-GLP-HPM-R',algorithms)
    alg = alg + 1;

    cd GLP
    t2=tic;
    I_MTF_GLP_HPM_R = MTF_GLP_HPM_R(I_MS,I_PAN,sensor,ratio);
    time_MTF_GLP_HPM_R = toc(t2);
    fprintf('Elaboration time MTF-GLP-HPM-R: %.2f [sec]\n',time_MTF_GLP_HPM_R);
    cd ..

    %%% Quality indexes computation
    [D_lambda_MTF_GLP_HPM_R,D_S_MTF_GLP_HPM_R,QNRI_MTF_GLP_HPM_R] = indexes_evaluation_FS(I_MTF_GLP_HPM_R,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MTF_GLP_HPM_R,D_S_MTF_GLP_HPM_R,QNRI_MTF_GLP_HPM_R];
    MatrixImage(:,:,:,alg) = I_MTF_GLP_HPM_R;
end

%% MTF-GLP-CBD

if ismember('MTF-GLP-CBD',algorithms)
    alg = alg + 1;
    
    cd GLP
    t2=tic;

    I_MTF_GLP_CBD = GS2_GLP(I_MS,I_PAN,ratio,sensor);

    time_MTF_GLP_CBD = toc(t2);
    fprintf('Elaboration time MTF-GLP-CBD: %.2f [sec]\n',time_MTF_GLP_CBD);
    cd ..

    [D_lambda_MTF_GLP_CBD,D_S_MTF_GLP_CBD,QNRI_MTF_GLP_CBD] = indexes_evaluation_FS(I_MTF_GLP_CBD,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MTF_GLP_CBD,D_S_MTF_GLP_CBD,QNRI_MTF_GLP_CBD];
    MatrixImage(:,:,:,alg) = I_MTF_GLP_CBD;
end

%% C-MTF-GLP-CBD

if ismember('C-MTF-GLP-CBD',algorithms)
    alg = alg + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PS_algorithm = 'GS2GLP'; % Pansharpening algorithm 
    n_segm = 5; % Number of segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd GS

    %%% Pansharpening Algorithm
    t2=tic;
    I_C_MTF_GLP_CBD = GS_Segm(I_MS,I_PAN,gen_LP_image(PS_algorithm,I_MS,I_PAN,I_MS_LR,ratio,sensor), k_means_clustering(I_MS,n_segm));
    time_C_MTF_GLP_CBD = toc(t2);
    fprintf('Elaboration time C-MTF-GLP-CBD: %.2f [sec]\n',time_C_MTF_GLP_CBD);
    cd ..

    %%% Quality indexes computation
    [D_lambda_C_MTF_GLP_CBD,D_S_C_MTF_GLP_CBD,QNRI_C_MTF_GLP_CBD] = indexes_evaluation_FS(I_C_MTF_GLP_CBD,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_C_MTF_GLP_CBD,D_S_C_MTF_GLP_CBD,QNRI_C_MTF_GLP_CBD];
    MatrixImage(:,:,:,alg) = I_C_MTF_GLP_CBD;
end

%% MF

if ismember('MF',algorithms)
    alg = alg + 1;
    
    cd MF
    t2=tic;
    I_MF= MF_HG_Pansharpen(I_MS,I_PAN,ratio);
    time_MF = toc(t2);
    fprintf('Elaboration time MF: %.2f [sec]\n',time_MF);
    cd ..

    [D_lambda_MF,D_S_MF,QNRI_MF] = indexes_evaluation_FS(I_MF,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_MF,D_S_MF,QNRI_MF];
    MatrixImage(:,:,:,alg) = I_MF;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variational Optimization Mathods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FE-HPM

if ismember('FE-HPM',algorithms)
    alg = alg + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tap = 25; % dimension of the support (at least 3*ratio) 
    num_iter_max = 5; % max number of iteration (at least 3; not sensitive)
    lambda = 10^5; % coefficient to weight the regularization term 
    mu = 10^5; % coefficient to weight the regularization term
    threshold = 10^(-3); % threshold on the kernel (it cuts to 0 values below threshold)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd FE-HPM
    t2 = tic;
    I_FE_HPM = FE_HPM(I_MS,I_PAN,ratio,max(tap,3*ratio),lambda,mu,threshold,max(num_iter_max,3),'Basic');
    time_FE_HPM = toc(t2);
    fprintf('Elaboration time FE-HPM: %.2f [sec]\n',time_FE_HPM);
    cd ..

    [D_lambda_FE_HPM,D_S_FE_HPM,QNRI_FE_HPM] = indexes_evaluation_FS(I_FE_HPM,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_FE_HPM,D_S_FE_HPM,QNRI_FE_HPM];
    MatrixImage(:,:,:,alg) = I_FE_HPM;
end

%% SR-D

if ismember('SR-D',algorithms)
    alg = alg + 1;
    
    %%%%%%%%%% Parameters setting
    TS = 7; % Tiling (dimensions of the patches are TS x TS)
    ol = 4; % Overlap (in pixels) between contiguous tile
    n_atoms = 10; % Max number of representation atoms (default value = 10)

    cd SR-D
    t2 = tic;
    I_SR_D = CS(I_MS,I_PAN,I_MS_LR,ratio,sensor,TS,ol,n_atoms);
    time_SR_D = toc(t2);
    fprintf('Elaboration time SR-D: %.2f [sec]\n',time_SR_D);
    cd ..

    [D_lambda_SR_D,D_S_SR_D,QNRI_SR_D] = indexes_evaluation_FS(I_SR_D,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_SR_D,D_S_SR_D,QNRI_SR_D];
    MatrixImage(:,:,:,alg) = I_SR_D;
end

%% PWMBF

if ismember('PWMBF',algorithms)
    alg = alg + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = 3;
    degrade=0;
    wavelet=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd PWMBF
    t2 = tic;
    I_PWMBF = PWMBF(I_PAN,I_MS_LR,ratio,r,wavelet,degrade);
    time_PWMBF = toc(t2);
    fprintf('Elaboration time PWMBF: %.2f [sec]\n',time_PWMBF);
    cd ..

    [D_lambda_PWMBF,D_S_PWMBF,QNRI_PWMBF] = indexes_evaluation_FS(I_PWMBF,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_PWMBF,D_S_PWMBF,QNRI_PWMBF];
    MatrixImage(:,:,:,alg) = I_PWMBF;
end

%% TV

if ismember('TV',algorithms)
    alg = alg + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch sensor
        case 'IKONOS'
            w=[0.1091    0.2127    0.2928    0.3854];
            c = 8;
            alpha=1.064;
            maxiter=10;
            lambda = 0.47106;
        case {'GeoEye1','WV4'}
            w=[0.1552, 0.3959, 0.2902, 0.1587];
            c = 8;
            alpha=0.75;
            maxiter=50;
            lambda = 157.8954;
        case 'WV3'
            w=[0.0657    0.1012    0.1537    0.1473    0.1245    0.1545    0.1338    0.1192];
            c = 8;
            alpha=0.75;
            maxiter=50;
            lambda = 1.0000e-03;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd TV
    t2 = tic;
    I_TV = TV_pansharpen(I_MS_LR,I_PAN,alpha,lambda,c,maxiter,w);
    time_TV = toc(t2);
    fprintf('Elaboration time TV: %.2f [sec]\n',time_TV);
    cd ..

    [D_lambda_TV,D_S_TV,QNRI_TV] = indexes_evaluation_FS(I_TV,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_TV,D_S_TV,QNRI_TV];
    MatrixImage(:,:,:,alg) = I_TV;
end

%% RR

if ismember('RR',algorithms)
    alg = alg + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Yim = cell(size(I_MS_LR,3),1);
    Yim{1} = I_PAN; 
    for ii = 2 : size(I_MS_LR,3) + 1
            Yim{ii} = I_MS_LR(:,:,ii-1);
    end
    
    mtf = 0.3 .* ones(1,size(I_MS_LR,3)+1);
    d = [1, 4.*ones(1,size(I_MS_LR,3))]';

    lambda = 1;
    CDiter = 1;
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
    I_RR = RRpansharp(Yim,'r',r,'lambda',lambda,'q',q','CDiter',CDiter,'tolgradnorm',tolgradnorm,'d',d,'mtf',mtf);
    time_RR = toc(t2);
    fprintf('Elaboration time RR: %.2f [sec]\n',time_RR);
    cd ..

    [D_lambda_RR,D_S_RR,QNRI_RR] = indexes_evaluation_FS(I_RR,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_RR,D_S_RR,QNRI_RR];
    MatrixImage(:,:,:,alg) = I_RR;
end

%% PNN

if ismember('PNN',algorithms)
    alg = alg + 1;
    
    cd PNN
    t2 = tic;
    I_PNN = PNN(I_MS_LR,I_PAN,sensor,L,false,I_MS);
    time_PNN = toc(t2);
    fprintf('Elaboration time PNN: %.2f [sec]\n',time_PNN);
    cd ..
    
    [D_lambda_PNN,D_S_PNN,QNRI_PNN] = indexes_evaluation_FS(I_PNN,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_PNN,D_S_PNN,QNRI_PNN];
    MatrixImage(:,:,:,alg) = I_PNN;
end

%% PNN-IDX

if ismember('PNN-IDX',algorithms)
    alg = alg + 1;
    
    cd PNN
    t2 = tic;
    I_PNN_IDX = PNN(I_MS_LR,I_PAN,sensor,L,true,I_MS);
    time_PNN_IDX = toc(t2);
    fprintf('Elaboration time PNN-IDX: %.2f [sec]\n',time_PNN_IDX);
    cd ..
    
    [D_lambda_PNN_IDX,D_S_PNN_IDX,QNRI_PNN_IDX] = indexes_evaluation_FS(I_PNN_IDX,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_PNN_IDX,D_S_PNN_IDX,QNRI_PNN_IDX];
    MatrixImage(:,:,:,alg) = I_PNN_IDX;
end

%% A-PNN

if ismember('A-PNN',algorithms)
    alg = alg + 1;
    
    cd PNN
    t2 = tic;
    I_A_PNN = PNNplus(I_MS_LR, I_PAN, sensor, 0, L, [], I_MS);
    time_A_PNN = toc(t2);
    fprintf('Elaboration time A-PNN: %.2f [sec]\n',time_A_PNN);
    cd ..
    
    [D_lambda_A_PNN,D_S_A_PNN,QNRI_A_PNN] = indexes_evaluation_FS(I_A_PNN,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_A_PNN,D_S_A_PNN,QNRI_A_PNN];
    MatrixImage(:,:,:,alg) = I_A_PNN;
end

%% A-PNN-FT

if ismember('A-PNN-FT',algorithms)
    alg = alg + 1;
    
    cd PNN
    t2 = tic;
    I_A_PNN_FT = PNNplus(I_MS_LR, I_PAN, sensor, 50, L, [], I_MS);
    time_A_PNN_FT = toc(t2);
    fprintf('Elaboration time A-PNN-FT: %.2f [sec]\n',time_A_PNN_FT);
    cd ..
    
    [D_lambda_A_PNN_FT,D_S_A_PNN_FT,QNRI_A_PNN_FT] = indexes_evaluation_FS(I_A_PNN_FT,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,ratio,flagQNR);

    MatrixResults(alg,:) = [D_lambda_A_PNN_FT,D_S_A_PNN_FT,QNRI_A_PNN_FT];
    MatrixImage(:,:,:,alg) = I_A_PNN_FT;
end

%% Print in LATEX
if flagQNR == 1
    matrix2latex(MatrixResults,'FR_Assessment.tex', 'rowLabels',algorithms,'columnLabels',[{'DL'},{'DS'},{'QNR'}],'alignment','c','format', '%.4f');
else
    matrix2latex(MatrixResults,'FR_Assessment.tex', 'rowLabels',algorithms,'columnLabels',[{'DL'},{'DS'},{'HQNR'}],'alignment','c','format', '%.4f');
end

%% View All

if size(I_MS,3) == 4
    vect_index_RGB = [3,2,1];
else
    vect_index_RGB = [5,3,2];
end

titleImages = algorithms;

figure, showImagesAll(MatrixImage,titleImages,vect_index_RGB,flag_cut_bounds,dim_cut,0);