import FFTW: fft, ifft, fftshift, ifftshift
using Statistics: mean
using ImageFiltering
using LinearAlgebra
using ImageTransformations
using ImageFiltering
using Statistics

function S2sharp(Yim,varargin)
#   Input:
#           Yim : 1x12 cell array containing the observed images for each of
#                 the nb bands of the Sentinel 2 sensor.
#         CDiter: Number of cyclic descent iterations.
#                 CDiter=10 is the default.
#              r: The subspace dimension, r=7 is the default.
#         lambda: The regularization parameter, lambda=0.005 is the
#                 default.
#          Xm_im: nl x nc x nb 3D matrix containing the true (10m resolution)
#              q: penalty weights, if r=7 then q= [1, 1.5, 4, 8, 15, 15, 20 ]'
#                 is the default otherwise the default is q=ones(p,1). Note
#                 that lam_i = lam*q_i
#             X0: Initial value for X = G * F'
#     Gstep_only: If Gstep_only=1 then perform the G-step (once). Assuming that F is fixed
#            GCV: If GCV=1 then the GCV value is computed.

    #declare the different inputs that aren't called in function
    CDiter=10
    r=7
    lambda=0.005
    Xm_im=""
    X0 = ""
    tolgradnorm = 0.1
    Gstep_only=0
    GCV = 0

    if r==7
        q = [1, 1.5, 4, 8, 15, 15, 20 ]'
    else
        println("NEED TO SUPPORT DIFFERENT r")
        q = (zeros(r,1).+1)'
    end

    output = Dict("SAMm" => [], "SAMm_2m" => [], "SRE" => [], "GVCscore" => [], "ERGAS_20m" => [], "ERGAS_60m" => [], "SSIM" => [], "aSSIM" => [], "RMSE" => [], "Time" => [] )

    for i in 1:2:(length(varargin))
        current = varargin[i]
        next = varargin[i+1]

        if current == "CDiter"
            CDiter = next
        elseif current == "r"
            r = next
        elseif current == "lambda"
            lambda = next
        elseif current == "Xm_im"
            Xm_im = next
        elseif current == "q"
            q = next
        elseif current == "X0"
            q = next
        elseif current == "tolgradnorm"
            tolgradnorm = next
        elseif current == "Gstep_only"
            Gstep_only = next
        elseif current == "GCV"
            GCV = next
        else
            println("Bad varargin: $current")

        end
    end

    #quick sanity checks here
    if length(q) != r
        println("The length of q has to match r")
    end

    #unsure if this will be needed
    q = q[:]

    ##NEED A DIFFERENT LENGTH FUNCTION THAN MATLAB HERE
    ##TODO check if this works
    #Yim = vec(Yim)
    #may need to convert data type here

    ##can't visualize what this does

    nl,nc = size(Yim[2]);
    L = length(Yim)
    n = nl*nc;


    Yim2, av = normaliseData(Yim)

    #subsampling factors (in pixels)
    d = [6 1 1 1 2 2 2 1 2 6 2 2]#'

    # convolution  operators (Gaussian convolution filters), taken from ref [5] from the Ulfarsson
    mtf = [ .32 .26 .28 .24 .38 .34 .34 .26 .33 .26 .22 .23]

    sdf = d.*sqrt.(-2*log.(mtf)/(π^2))

    #julia sucks bc i use python too much
    for i in 1:length(sdf)
        if d[i] ==1
            sdf[i]=0
        end
    end

    limsub=2
    dx=12
    dy=12
    FBM = createConvKernel(sdf,d,nl,nc,L,dx,dy)
    Y,M,F = initialization(Yim2,sdf,nl,nc,L,dx,dy,d,limsub,r)
    Mask = reshape(M,(n,L))'

    if X0 == ""
        Z = zeros(r,n)
    else
        #low rank approx
        X0, = normaliseData(X0)
        X0 = reshape(X0, (n,L))'
        U,Σ,V = svd(X0, full=false)
        U = U[:,1:r]
        Z = Σ[1:r,1:r]*V[:,1:r]'
    end

    FDH,FDV,FDHC,FDVC = createDiffkernels(nl,nc,r)
    sigmas = 1.0
    W = computeWeights(Y,d,sigmas,nl)
    Whalf = sqrt.(W)

    if GCV==1
        Gstep_only=1
    end

    if Gstep_only ≈ 1
        CDiter=1
    end

    # Konrad 4/14: stopped here

    for jCD in 1:CDiter
        Z, Jcost[jCD], options = Zstep(Y,FBM,F,lambda,nl,nc,Z,Mask,q,FDH,FDV,FDHC,FDVC,W,Whalf,tolgradnorm)

        if(Gstep_only==0)
           F1=Fstep(F,Z,Y,FBM,nl,nc,Mask)
           F=F1
        end

        if( GCV==1 )
        Ynoise = ( abs(Y) > 0 ) .* randn( size(Y) )

            Znoise = Zstep(Ynoise,FBM,F,lambda,nl,nc,Z,Mask,q,FDH,FDV,FDHC,FDVC,W,Whalf,tolgradnorm)

            HtHBXnoise = Mask.*ConvCM(F*Znoise,FBM,nl)

            Ynoise = Ynoise[[1,5,6,7,9,10,11,12],:]
            HtHBXnoise = HtHBXnoise[[1,5,6,7,9,10,11,12],:]


            den = tr(Ynoise*(Ynoise - HtHBXnoise)')
            HtHBX=Mask.*ConvCM(F*Z,FBM,nl)

            num = norm(AY[[1,5,6,7,9,10,11,12],:] - HtHBX[[1,5,6,7,9,10,11,12],:] , p=2).^2;

            output["GCVscore"] = [num / den]
        end
        #output["Time"] = [toc] unclear where they get this from in the Matlab

        if Xm_im != [] #idk if this is good syntax for isempty, but its the best i have
            Xhat_im = conv2im(F*Z,nl,nc,L)
            output["SAMm"], output["SAMm_2m"], output["SRE"], output["RMSE"], output["SSIM"], output["aSSIM"],output["ERGAS_20m"], output["ERGAS_60m"] = evaluate_performance(Xm_im,Xhat_im,nl,nc,L,limsub,d,av);
        end
    end
    Xhat_im = conv2im(F*Z,nl,nc,L)

    Xhat_im = Xhat_im[limsub+1:end-limsub,limsub+1:end-limsub,:]
    Xhat_im = unnormaliseData(Xhat_im,av)



#   Output:   output is a structure containing the following fields
#      Xhat_im: estimated image (3D) at high resolution (10m) for each
#               spectral channel
#         SAMm: mean SAM for the 60m and 20m bands (empty if Xm_im is not
#               available)
#      SAMm_2m: mean SAM for the 20m bands (empty if Xm_im is not available)
#          SRE: signal to reconstruction error for all the 12 bands
#               (empty if Xm_im is not available)
#          GCVscore: Contains the GCV score if it was GCV=1 otherwise GCV is
#          empty.
#      ERGAS_20m: ERGAS score for the 20 m bands
#      ERGAS_60m: ERGAS score for the 60 m bands
#          aSSIM: average Structural Similarity Index
#           RMSE: Root mean squared error
#           Time: computational time


    return Xhat_im, output
end
function initialization(Yim2::Array{Any}, sdf::Array{Float64},
    nl::Int,nc::Int,L::Int,
    dx::Int, dy::Int, d::Array{Int}, limsub::Int,r::Int)

    FBM2 = createConvKernelSubspace(sdf,nl,nc,L,dx,dy)

    Ylim = zeros(nl, nc, L)

    for i = 1:L
        Ylim[:,:,i] = imresize(Yim2[i], nl, nc)
    end

    Y2im = real.(ifft(fft(Ylim).*FBM2))
    Y2tr = Y2im[limsub+1:end-limsub,limsub+1:end-limsub,:]
    Y2n = reshape(Y2tr,(nl-4)*(nc-4),L)
    F, dd, P = svd(Y2n')#,"econ")
    F = F[:,1:r]
    M, Y = createSubsampling(Yim2, d, nl, nc, L)

    return Y, M, F
end

function createSubsampling(Yim::Array{Any}, d::Array{Int}, nl::Int, nc::Int, L::Int)
M = zeros(nl, nc, L)
indexes = Array{Any}(undef,L)
Y = zeros(L,nl*nc)
for i = 1:L
    im = ones(nl÷d[i],nc÷d[i])
    maux = zeros(d[i],d[i])
    maux[1,1] = 1
    M[:,:,i] = kron(im,maux)
    indexes[i] = findall(!iszero,vec(M[:,:,i]))
    Y[i,indexes[i]] = conv2mat(Yim[i],nl÷d[i],nc÷d[i],1)
end
return M, Y
end


function ConvCM(X::Array{Float64,2}, FKM::Array{ComplexF64,3}, nl::Int, nc::Int = 0, L::Int=0)
    if nc == 0 || L == 0
        L, n = size(X)
        nc = n ÷ nl
    end
    X = conv2mat(real.(ifft(fft(conv2im(X,nl,nc,L)).*FKM)))
    return X
end


function conv2mat(X,nl::Int=0,nc::Int=0,L::Int=0)
    if ndims(X) == 3
        nl,nc,L = size(X);
        X = reshape(X,nl*nc,L)';
    elseif ndims(squeeze(X)) == 2
        L = 1;
        nl,nc = size(X);
        X = reshape(X,nl*nc,L)';
    end
    return X
end


function grad_cost_G(Z::Array{Float64,2},
    F::Array{Float64,2},
    Y::Array{Float64,2},
    UBTMTy::Array{Float64,2},
    FBM::Array{ComplexF64,3},
    Mask::Array{Float64,2},
    nl::Int, nc::Int, r::Int,
    tau::Float64, q::Array{Float64},
    FDH::Array{ComplexF64,3},
    FDV::Array{ComplexF64,3},
    FDHC::Array{ComplexF64,3},
    FDVC::Array{ComplexF64,3},
    W::Array{Float64})

    X = F*Z
    BX=ConvCM(X,FBM,nl);
    HtHBX=Mask.*BX;
    ZH=ConvCM(Z,FDHC,nl);
    Zv=ConvCM(Z,FDVC,nl);
    ZHW=ZH.*W;
    ZVW=Zv.*W;
    grad_pen=ConvCM(ZHW,FDH,nl)+ConvCM(ZVW,FDV,nl);
    AtAg = F'*ConvCM(HtHBX,conj.(FBM),nl)+2*tau*(q*ones(1,nl*nc)).*grad_pen;
    gradJ=AtAg-UBTMTy;
    J = 1/2 * sum( Z .* AtAg ) - sum( Z.*UBTMTy );

    return J, gradJ, AtAg
end

function CG(Z::Array{Float64,2},
    F::Array{Float64,2},
    Y::Array{Float64,2},
    UBTMTy::Array{Float64,2},
    FBM::Array{ComplexF64,3},
    Mask::Array{Float64,2},
    nl::Int, nc::Int, r::Int,
    tau::Float64, q::Array{Float64},
    FDH::Array{ComplexF64,3},
    FDV::Array{ComplexF64,3},
    FDHC::Array{ComplexF64,3},
    FDVC::Array{ComplexF64,3},
    W::Array{Float64})

    maxiter = 1000;
    tolgradnorm = 0.1;  #%1e-6;   %0.1
    cost,grad = grad_cost_G(Z,F,Y,UBTMTy,FBM,Mask,nl,nc,r,tau,q,FDH,FDV,FDHC,FDVC,W);
    gradnorm = norm(grad[:]);
    iter = 0;
    res = -grad;
    while ( gradnorm > tolgradnorm & iter < maxiter )
        iter = iter + 1;
       # fprintf('%5d\t%+.16e\t%.8e\n', iter, cost, gradnorm);
        if( iter == 1 )
            desc_dir = res;
        else
            beta = ( res[:]' * res[:] ) / ( old_res[:]' * old_res[:] );
            desc_dir = res + beta * desc_dir;
        end
        _, _, AtAp = grad_cost_G(desc_dir,F,Y,UBTMTy,FBM,Mask,nl,nc,r,tau,q,FDH,FDV,FDHC,FDVC,W);
        alpha = ( res[:]' * res[:] ) / ( desc_dir[:]' * AtAp[:] );
        Z1 = Z + alpha * desc_dir;
        old_res = res;
        res = res - alpha* AtAp;
        gradnorm = norm( res[:] );
        # Transfer iterate info
        Z = Z1;
    end

end


function computeWeights(Y::Array{Float64,2}, d::Array{Int}, sigmas::Float64, nl::Int, method::String="")
    hr_bands = d==1;
    hr_bands = findall(hr_bands)';
    L = size(Y,1)
    nc = size(Y,2)÷nl
    grad = zeros(nl,nc,L)
    for i=hr_bands
        #TODO get intermediate to work
         Gy, Gx = imgradients(conv2im(Y[i,:],nl),KernelFactors.sobel);
         grad[:,:,i] = Gy.^2 .+ Gx.^2
        # switch method
        #     case 'original'
        #         grad[:,:,i] = imgradient(conv2im(Y[i,:],nl),'intermediate').^2;
        #     case 'sobel'
        #         grad[:,:,i] = (1/16)*(imgradient(conv2im(Y[i,:],nl),'sobel').^2);
        #     case 'prewitt'
        #         grad[:,:,i] = (1/9)*(imgradient(conv2im(Y[i,:],nl),'prewitt').^2);
        #
        # end
    end
    grad = sqrt.(maximum(grad,dims=3));
    grad = grad / quantile(grad[:],0.95);

    Wim = exp.(-grad.^2/2/sigmas^2);
    Wim[Wim.<0.5] .= 0.5;

    W = conv2mat(Wim,nl);

    return W
end


function Zstep(
    Y::Array{Float64,2},
    FBM::Array{ComplexF64,3},
    F::Array{Float64,2},
    tau::Float64, nl::Int, nc::Int,
    Z::Array{Float64,2},
    Mask::Array{Float64,2},
    q::Array{Float64},
    FDH::Array{ComplexF64,3},
    FDV::Array{ComplexF64,3},
    FDHC::Array{ComplexF64,3},
    FDVC::Array{ComplexF64,3},
    W::Array{Float64},
    Whalf::Array{Float64},
    tolgradnorm::Float64)


    r = size(F,2);
    n = nl*nc;
    UBTMTy=F'*ConvCM(Y,conj.(FBM),nl);
    Z = CG(Z,F,Y,UBTMTy,FBM,Mask,nl,nc,r,tau,q,FDH,FDV,FDHC,FDVC,W);
    xcost=1;
    options=[];

    return Z, xcost, options

end



function Fstep(F::Array{Float64,2},
        Z::Array{Float64,2},
        Y::Array{Float64,2},
        FBM::Array{ComplexF64,3},
        nl::Int,nc::Int,
        Mask::Array{Float64,2})
     F0=F;#%   U; % initialization
     BTXhat =  ConvCM(F0*Z,FBM,nl);
     MBTXhat=Mask.*BTXhat;
     L, r = size(F);
     for ii=1:L
        MBZT[:,:,ii] = repeat(Mask(ii,:),[r,1]).*
                        ConvCM(Z,repeat(FBM[:,:,ii],outer=[1,1,r]),outer=nl);
        A[:,:,ii]=MBZT[:,:,ii]*MBZT[:,:,ii]';
        ZBMTy[:,ii]=MBZT[:,:,ii]*Y[ii,:]';
     end
     ZBYT=ZBMTy';#   BTY*Z';

     manifold = Manopt.Stiefel(L,r)
     cost = F -> costF(F, MBZT, Y)
     grad = F -> egrad(F,A,ZBYT)

     F1 = Manopt.trust_regions(manifold, cost, grad,
            stopping_criterion = StopWhenGradientNormLess(1e-2))

    return F1

     #manifold = stiefelfactory(L,r,1);# %euclideanfactory(L,r);
     #problem.M = manifold;
     #problem.cost  = @(F) costF(F,MBZT,Y);
     #problem.egrad = @(F) egrad(F,A,ZBYT);
     #warning('off', 'manopt:getHessian:approx')
     #options.tolgradnorm = 1e-2;
     #options.verbosity=0;
     #[F1, xcost, info, options] = trustregions(problem,F0,options);
end



function costF(F,MBZT,Y)
    L=size(F,1);
    Ju=0;
    for i=1:L,
        fi=F[i,:]';
        yi=Y[i,:]';
        Ju=Ju+0.5*norm(MBZT[:,:,i]'*fi-yi,"fro")^2;
    end

    return Ju
end

function  egrad(F,A,ZBYT)
    p=size(A,3);
    Du=0*F;
    for ii=1:p
        Du[ii,:]=F[ii,:]*A[:,:,ii]'-ZBYT[ii,:];
    end

    return Du
end

function ERGAS(I1::Array{Float64,3},
    I2::Array{Float64,3},
    ratio::Int)

    #I1 = double(I1);
    #I2 = double(I2);

    Err=I1-I2;
    ERGAS_index=0.0;
    for iLR = 1:size(Err,3),
        ERGAS_index = ERGAS_index+mean(Err[:,:,iLR].^2)/(mean((I1[:,:,iLR])))^2;
    end

    ERGAS_index = (100/ratio) * sqrt((1/size(Err,3)) * ERGAS_index);

    return ERGAS_index
end

function  evaluate_performance(
        Xm_im::Array{Float64,3},
        Xhat_im::Array{Float64,3},
        nl::Int,nc::Int,L::Int,
        limsub::Int,
        d::Array{Int},
        av::Array{Float64})

    Xhat_im = Xhat_im[limsub+1:end-limsub,limsub+1:end-limsub,:];
    Xhat_im = unnormaliseData(Xhat_im,av);
    Xhat=reshape(Xhat_im,[(nl-4)*(nc-4),L]);
    #% Xm_im is the ground truth image
    Xm_im = Xm_im[limsub+1:end-limsub,limsub+1:end-limsub,:];
    if ( size(Xm_im,3) == 6 ) #% Reduced Resolution
        ind = findall( d==2 );
        SAMm=SAM(Xm_im,Xhat_im[:,:,ind]);
        SAMm_2m=SAMm;
        X = conv2mat(Xm_im);
        Xhat = conv2mat(Xhat_im);
        #% SRE - signal to reconstrution error
        for i=1:6
            SRE[i,1] = 10*log10(sum(X[i,:].^2)/ sum((Xhat(ind(i),:)-X[i,:]).^2));
            SSIM_index[i,1] = ssim(Xm_im[:,:,i],Xhat_im(:,:,ind(i)));
        end
        aSSIM=mean(SSIM_index);
        ERGAS_20m = ERGAS(Xm_im,Xhat_im[:,:,ind],2);
        ERGAS_60m = NaN;
        RMSE = norm(X - Xhat[ind,:],"fro") / size(X,2);
    else
        ind=findall(d .== 2 | d .== 6);
        SAMm=SAM(Xm_im[:,:,ind],Xhat_im[:,:,ind]);
        ind2=findall(d .== 2);
        SAMm_2m=SAM(Xm_im[:,:,ind2],Xhat_im[:,:,ind2]);
        ind6=findall(d .== 6);
        X = conv2mat(Xm_im);
        Xhat = conv2mat(Xhat_im);
        #% SRE - signal to reconstrution error
        for i=1:L
            SRE[i,1] = 10*log10(sum(X[i,:].^2)/ sum((Xhat[i,:]-X[i,:]).^2));
            SSIM_index[i,1] = ssim(Xm_im[:,:,i],Xhat_im[:,:,i]);
        end
        aSSIM=mean(SSIM_index(ind));
        ERGAS_20m = ERGAS(Xm_im[:,:,ind],Xhat_im[:,:,ind],2);
        ERGAS_60m = ERGAS(Xm_im[:,:,ind2],Xhat_im[:,:,ind2],6);
        RMSE = norm(X[ind,:] - Xhat[ind,:],"frp") / size(X,2);
    end


    return SAMm, SAMm_2m, SRE, RMSE, SSIM_index, aSSIM, ERGAS_20m, ERGAS_60m
end
function createDiffkernels(nl,nc,r)
    dh = zeros(nl,nc)
    dh[1,1]=1
    dh[1,nc]=-1

    dv = zeros(nl,nc)
    dv[1,1]=1
    dv[nl,1]=-1

    fft_dh = fft(dh)
    fft_dv = fft(dv)

    FDH = zeros(ComplexF64,nl,nc,r)
    FDV = zeros(ComplexF64,nl,nc,r)

    for i in 1:r
        FDH[:,:,i] = fft_dh
        FDV[:,:,i] = fft_dv
    end

    FDHC = conj(FDH)
    FDVC = conj(FDV)
    return FDH,FDV,FDHC,FDVC
end


function normaliseData(Yim::Array{Any})
    #Yim is the array that has all the different images in it

    nb = length(Yim)
    av = zeros(nb)

    #TODO see what this other case even means
    for i in 1:nb
        av[i] = mean(mean(Yim[i].^2))
        Yim[i] = sqrt.(Yim[i].^2/av[i])
    end

    return Yim, av
end

function createConvKernel(sdf,d,nl,nc,L,dx,dy)


    ##might not need this round here. i added this before seeing how round is done below
    middlel = nl÷2
    middlec = nc÷2

    B = zeros(nl,nc,L)
    FBM = zeros(ComplexF64,nl,nc,L)

    for i in 1:L
        if d[i] >1
            h = Kernel.gaussian(3) #this will produce a 13x13 kernel which is close the example 12x12 needed

            h = imresize(h,12,12)# HACK

            #this function won't work.. not adaptable enough
            #need to create a 12x12 filter here...

            B[(middlel-dy÷2+1:middlel+dy÷2) .- (d[i]÷2) .+ 1 , (middlec-dx÷2+1:middlec+dx÷2) .- (d[i]÷2) .+ 1, i] = h #not sure about this line here

            B[:,:,i]=fftshift( B[:,:,i] )/sum(sum(B[:,:,i]))
            FBM[:,:,i] = fft(B[:,:,i])
        else
            B[1,1,i]=1
            FBM[:,:,i]=fft(B[:,:,i])
        end
    end

    return FBM
end



function createConvKernelSubspace(sdf,nl,nc,L,dx,dy)
    middlel=(nl+1)÷2
    middlec=(nc+1)÷2

    dx = dx+1
    dy = dy+1

    # kernel filters expanded to size [nl,nc]
    B = zeros(nl,nc,L)
    # fft2 of kernels
    FBM2 = zeros(ComplexF64,nl,nc,L)
    s2 = maximum(sdf)

    ##this loop is nearly identical to whats happening above! this would be a good place to improve it
    for i in 1:L
        if sdf[i] < s2 # != would be a little better i think

            h = Kernel.gaussian(3) #this will produce a 13x13 kernel which is close the example 12x12 needed

            B[middlel-(dy-1)÷2:middlel+(dy-1)÷2,middlec-(dx-1)÷2:middlec+(dx-1)÷2,i] = h;
            B[:,:,i]=fftshift( B[:,:,i] )/sum(sum(B[:,:,i]))
            FBM2[:,:,i] = fft(B[:,:,i])
        else
            B[1,1,i]=1
            FBM2[:,:,i]=fft(B[:,:,i])
        end

    end

    return FBM2

end


function conv2mat(X,nl::Int,nc::Int,L::Int)
    # function X = conv2mat(X,nl,nc,L)
    #     if ndims(X) == 3
    #         [nl,nc,L] = size(X);
    #         X = reshape(X,nl*nc,L)';
    #     elseif ndims(squeeze(X)) == 2
    #         L = 1;
    #         [nl,nc] = size(X);
    #         X = reshape(X,nl*nc,L)';
    #     end
    # end

    if ndims(X) == 3
        nl,nc,L = size(X)
        X = reshape(X,nl*nc,L)'
        else #ndims(squeeze(X)) == 2 #unsure if squeeze() is in image transformations
        L = 1
        nl,nc = size(X)
        X = reshape(X,nl*nc,L)'
    end
    return X
end

function unnormaliseData(Yim::Array{Any},av::Array{Float64})
    # function [Yim] = unnormaliseData(Yim, av)
    #     if iscell(Yim)
    #         % mean squared power = 1
    #         nb = length(Yim);
    #         for i=1:nb
    #             Yim{i,1} = sqrt(Yim{i}.^2*av(i,1));
    #         end
    #     else
    #         nb = size(Yim,3);
    #         for i=1:nb
    #             Yim(:,:,i) = sqrt(Yim(:,:,i).^2*av(i,1));
    #         end
    #     end
    # end

    nb,=size(Yim)

    for i in 1:nb
        Yim[:,:,i] = sqrt.(Yim[i].^2 * av[i])
    end
    return Yim
end


function regularization(X1,X2,tau,mu,W,q)
    # function [Y1,Y2] = regularization(X1,X2,tau,mu,W,q)
    #     Wr = q*W;
    #     %Wr=ones(size(Wr));
    #     Y1 = (mu*X1)./(mu + tau*Wr);
    #     Y2 = (mu*X2)./(mu + tau*Wr);
    # end

    ##NOT USED IN THE EXAMPLE CODE
    ##WAS NOT ABLE TO CHECK THIS
    Wr = q*W;

    Y1 = (mu*X1)./(mu + tau*Wr);
    Y2 = (mu*X2)./(mu + tau*Wr);

    return Y1, Y2
end

function conv2im(X,nl,nc=0,L=0)

    #     if size(X,2)==1
    #         X = conv2mat(X,nl,nc,L);
    #     end
    #     if nargin == 2
    #         [L,n] = size(X);
    #         if n==1
    #             X = conv2mat(X,nl,nc,L);
    #         end
    #         nc = n/nl;
    #     end
    #     X = reshape(X',nl,nc,L);

    a,b=size(X)
    if b==1
        X = conv2mat(X,nl,nc,L)
    end

    if L==0 && nc==0 #this is the easiest sub for L,nc not being included
        nc = b/nl
    end

    return reshape(X',(nl,nc,L))
end
