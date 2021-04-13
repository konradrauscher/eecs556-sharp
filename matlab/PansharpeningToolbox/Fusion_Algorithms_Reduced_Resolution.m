%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  RUN AND REDUCED RESOLUTION VALIDATION  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select algorithms to run

algorithms = {'GT','EXP','BT-H','BDSD','C-BDSD','BDSD-PC','GS','GSA','C-GSA','PRACS','AWLP',...
        'MTF-GLP','MTF-GLP-FS','MTF-GLP-HPM','MTF-GLP-HPM-H','MTF-GLP-HPM-R','MTF-GLP-CBD','C-MTF-GLP-CBD','MF','FE-HPM','SR-D','PWMBF','TV','RR','PNN','PNN-IDX','A-PNN','A-PNN-FT'};

%% Initialization of the Matrix of Results
NumIndexes = 5;
MatrixResults = zeros(numel(algorithms),NumIndexes);
alg = 0;

%% Tools
addpath([pwd,'/Tools']);

%% MS

if size(I_GT,3) == 4   
    showImage4LR(I_MS_LR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);    
else
    showImage8LR(I_MS_LR,printEPS,1,flag_cut_bounds,dim_cut,thvalues,L,ratio);
end

%% PAN

showPan(I_PAN,printEPS,2,flag_cut_bounds,dim_cut);

%% GT

if ismember('GT',algorithms)
    alg = alg + 1;
    [Q_avg_GT, SAM_GT, ERGAS_GT, SCC_GT_GT, Q_GT] = indexes_evaluation(I_GT,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    MatrixResults(alg,:) = [Q_GT,Q_avg_GT,SAM_GT,ERGAS_GT,SCC_GT_GT];
    MatrixImage(:,:,:,alg) = I_GT;
end

%% EXP

if ismember('EXP',algorithms)
    alg = alg + 1;
    [Q_avg_EXP, SAM_EXP, ERGAS_EXP, SCC_GT_EXP, Q_EXP] = indexes_evaluation(I_MS,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    MatrixResults(alg,:) = [Q_EXP,Q_avg_EXP,SAM_EXP,ERGAS_EXP,SCC_GT_EXP];
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
    [Q_avg_BT_H, SAM_BT_H, ERGAS_BT_H, SCC_GT_BT_H, Q_BT_H] = indexes_evaluation(I_BT_H,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_BT_H,Q_avg_BT_H,SAM_BT_H,ERGAS_BT_H,SCC_GT_BT_H];
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

    [Q_avg_BDSD, SAM_BDSD, ERGAS_BDSD, SCC_GT_BDSD, Q_BDSD] = indexes_evaluation(I_BDSD,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_BDSD,Q_avg_BDSD,SAM_BDSD,ERGAS_BDSD,SCC_GT_BDSD];
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

    [Q_avg_C_BDSD, SAM_C_BDSD, ERGAS_C_BDSD, SCC_GT_C_BDSD, Q_C_BDSD] = indexes_evaluation(I_C_BDSD,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_C_BDSD,Q_avg_C_BDSD,SAM_C_BDSD,ERGAS_C_BDSD,SCC_GT_C_BDSD];
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

    [Q_avg_BDSD_PC, SAM_BDSD_PC, ERGAS_BDSD_PC, SCC_GT_BDSD_PC, Q_BDSD_PC] = indexes_evaluation(I_BDSD_PC,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_BDSD_PC,Q_avg_BDSD_PC,SAM_BDSD_PC,ERGAS_BDSD_PC,SCC_GT_BDSD_PC];
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

    [Q_avg_GS, SAM_GS, ERGAS_GS, SCC_GT_GS, Q_GS] = indexes_evaluation(I_GS,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_GS,Q_avg_GS,SAM_GS,ERGAS_GS,SCC_GT_GS];
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

    [Q_avg_GSA, SAM_GSA, ERGAS_GSA, SCC_GT_GSA, Q_GSA] = indexes_evaluation(I_GSA,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_GSA,Q_avg_GSA,SAM_GSA,ERGAS_GSA,SCC_GT_GSA];
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
    [Q_avg_C_GSA, SAM_C_GSA, ERGAS_C_GSA, SCC_GT_C_GSA, Q_C_GSA] = indexes_evaluation(I_C_GSA,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_C_GSA,Q_avg_C_GSA,SAM_C_GSA,ERGAS_C_GSA,SCC_GT_C_GSA];
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

    [Q_avg_PRACS, SAM_PRACS, ERGAS_PRACS, SCC_GT_PRACS, Q_PRACS] = indexes_evaluation(I_PRACS,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_PRACS,Q_avg_PRACS,SAM_PRACS,ERGAS_PRACS,SCC_GT_PRACS];
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

    [Q_avg_AWLP, SAM_AWLP, ERGAS_AWLP, SCC_GT_AWLP, Q_AWLP] = indexes_evaluation(I_AWLP,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_AWLP,Q_avg_AWLP,SAM_AWLP,ERGAS_AWLP,SCC_GT_AWLP];
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

    [Q_avg_MTF_GLP, SAM_MTF_GLP, ERGAS_MTF_GLP, SCC_GT_MTF_GLP, Q_MTF_GLP] = indexes_evaluation(I_MTF_GLP,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MTF_GLP,Q_avg_MTF_GLP,SAM_MTF_GLP,ERGAS_MTF_GLP,SCC_GT_MTF_GLP];
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

    [Q_avg_MTF_GLP_FS, SAM_MTF_GLP_FS, ERGAS_MTF_GLP_FS, SCC_GT_MTF_GLP_FS, Q_MTF_GLP_FS] = indexes_evaluation(I_MTF_GLP_FS,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MTF_GLP_FS,Q_avg_MTF_GLP_FS,SAM_MTF_GLP_FS,ERGAS_MTF_GLP_FS,SCC_GT_MTF_GLP_FS];
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

    [Q_avg_MTF_GLP_HPM, SAM_MTF_GLP_HPM, ERGAS_MTF_GLP_HPM, SCC_GT_MTF_GLP_HPM, Q_MTF_GLP_HPM] = indexes_evaluation(I_MTF_GLP_HPM,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MTF_GLP_HPM,Q_avg_MTF_GLP_HPM,SAM_MTF_GLP_HPM,ERGAS_MTF_GLP_HPM,SCC_GT_MTF_GLP_HPM];
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
    [Q_avg_MTF_GLP_HPM_H, SAM_MTF_GLP_HPM_H, ERGAS_MTF_GLP_HPM_H, SCC_GT_MTF_GLP_HPM_H, Q_MTF_GLP_HPM_H] = indexes_evaluation(I_MTF_GLP_HPM_H,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MTF_GLP_HPM_H,Q_avg_MTF_GLP_HPM_H,SAM_MTF_GLP_HPM_H,ERGAS_MTF_GLP_HPM_H,SCC_GT_MTF_GLP_HPM_H];
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
    [Q_avg_MTF_GLP_HPM_R, SAM_MTF_GLP_HPM_R, ERGAS_MTF_GLP_HPM_R, SCC_GT_MTF_GLP_HPM_R, Q_MTF_GLP_HPM_R] = indexes_evaluation(I_MTF_GLP_HPM_R,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MTF_GLP_HPM_R,Q_avg_MTF_GLP_HPM_R,SAM_MTF_GLP_HPM_R,ERGAS_MTF_GLP_HPM_R,SCC_GT_MTF_GLP_HPM_R];
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

    [Q_avg_MTF_GLP_CBD, SAM_MTF_GLP_CBD, ERGAS_MTF_GLP_CBD, SCC_GT_MTF_GLP_CBD, Q_MTF_GLP_CBD] = indexes_evaluation(I_MTF_GLP_CBD,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MTF_GLP_CBD,Q_avg_MTF_GLP_CBD,SAM_MTF_GLP_CBD,ERGAS_MTF_GLP_CBD,SCC_GT_MTF_GLP_CBD];
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
    [Q_avg_C_MTF_GLP_CBD, SAM_C_MTF_GLP_CBD, ERGAS_C_MTF_GLP_CBD, SCC_GT_C_MTF_GLP_CBD, Q_C_MTF_GLP_CBD] = indexes_evaluation(I_C_MTF_GLP_CBD,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_C_MTF_GLP_CBD,Q_avg_C_MTF_GLP_CBD,SAM_C_MTF_GLP_CBD,ERGAS_C_MTF_GLP_CBD,SCC_GT_C_MTF_GLP_CBD];
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

    [Q_avg_MF, SAM_MF, ERGAS_MF, SCC_GT_MF, Q_MF] = indexes_evaluation(I_MF,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_MF,Q_avg_MF,SAM_MF,ERGAS_MF,SCC_GT_MF];
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

    [Q_avg_FE_HPM, SAM_FE_HPM, ERGAS_FE_HPM, SCC_GT_FE_HPM, Q_FE_HPM] = indexes_evaluation(I_FE_HPM,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_FE_HPM,Q_avg_FE_HPM,SAM_FE_HPM,ERGAS_FE_HPM,SCC_GT_FE_HPM];
    MatrixImage(:,:,:,alg) = I_FE_HPM;
end

%% SR-D

if ismember('SR-D',algorithms)
    alg = alg + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TS = 7; % Tiling (dimensions of the patches are TS x TS)
    ol = 4; % Overlap (in pixels) between contiguous tile
    n_atoms = 10; % Max number of representation atoms (default value = 10)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd SR-D
    t2 = tic;
    I_SR_D = CS(I_MS,I_PAN,I_MS_LR,ratio,sensor,TS,ol,n_atoms);
    time_SR_D = toc(t2);
    fprintf('Elaboration time SR-D: %.2f [sec]\n',time_SR_D);
    cd ..

    [Q_avg_SR_D, SAM_SR_D, ERGAS_SR_D, SCC_GT_SR_D, Q_SR_D] = indexes_evaluation(I_SR_D,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_SR_D,Q_avg_SR_D,SAM_SR_D,ERGAS_SR_D,SCC_GT_SR_D];
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

    [Q_avg_PWMBF, SAM_PWMBF, ERGAS_PWMBF, SCC_GT_PWMBF, Q_PWMBF] = indexes_evaluation(I_PWMBF,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_PWMBF,Q_avg_PWMBF,SAM_PWMBF,ERGAS_PWMBF,SCC_GT_PWMBF];
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

    [Q_avg_TV, SAM_TV, ERGAS_TV, SCC_GT_TV, Q_TV] = indexes_evaluation(I_TV,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_TV,Q_avg_TV,SAM_TV,ERGAS_TV,SCC_GT_TV];
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

    [Q_avg_RR, SAM_RR, ERGAS_RR, SCC_GT_RR, Q_RR] = indexes_evaluation(I_RR,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);

    MatrixResults(alg,:) = [Q_RR,Q_avg_RR,SAM_RR,ERGAS_RR,SCC_GT_RR];
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
    
    [Q_avg_PNN, SAM_PNN, ERGAS_PNN, SCC_GT_PNN, Q_PNN] = indexes_evaluation(I_PNN,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    
    MatrixResults(alg,:) = [Q_PNN,Q_avg_PNN,SAM_PNN,ERGAS_PNN,SCC_GT_PNN];
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
    
    [Q_avg_PNN_IDX, SAM_PNN_IDX, ERGAS_PNN_IDX, SCC_GT_PNN_IDX, Q_PNN_IDX] = indexes_evaluation(I_PNN_IDX,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    
    MatrixResults(alg,:) = [Q_PNN_IDX,Q_avg_PNN_IDX,SAM_PNN_IDX,ERGAS_PNN_IDX,SCC_GT_PNN_IDX];
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
    
    [Q_avg_A_PNN, SAM_A_PNN, ERGAS_A_PNN, SCC_GT_A_PNN, Q_A_PNN] = indexes_evaluation(I_A_PNN,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    
    MatrixResults(alg,:) = [Q_A_PNN,Q_avg_A_PNN,SAM_A_PNN,ERGAS_A_PNN,SCC_GT_A_PNN];
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
    
    [Q_avg_A_PNN_FT, SAM_A_PNN_FT, ERGAS_A_PNN_FT, SCC_GT_A_PNN_FT, Q_A_PNN_FT] = indexes_evaluation(I_A_PNN_FT,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    
    MatrixResults(alg,:) = [Q_A_PNN_FT,Q_avg_A_PNN_FT,SAM_A_PNN_FT,ERGAS_A_PNN_FT,SCC_GT_A_PNN_FT];
    MatrixImage(:,:,:,alg) = I_A_PNN_FT;
end

%% Print in LATEX

matrix2latex(MatrixResults(:,[1,3,4]),'RR_Assessment.tex', 'rowLabels',algorithms,'columnLabels',[{'Q2n'},{'SAM'},{'ERGAS'}],'alignment','c','format', '%.4f');

%% View All

if size(I_GT,3) == 4
    vect_index_RGB = [3,2,1];
else
    vect_index_RGB = [5,3,2];
end

titleImages = algorithms;
  
figure, showImagesAll(MatrixImage,titleImages,vect_index_RGB,flag_cut_bounds,dim_cut,0);