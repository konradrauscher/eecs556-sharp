using FFTW
using Statistics: mean
using ImageFiltering
using Images
using LinearAlgebra
using ImageTransformations
using ImageFiltering
using Statistics
using OffsetArrays
using Printf
#using MIRT: diff_forw
using MATLAB # :(
MATLAB.__init__()#absolutely unbelievable that I have to put this!!!

function S2sharp(Yim;
    CDiter = 10,
    r=7,
    lambda=0.05,
    Xm_im=nothing,
    q=[1, 1.5, 4, 8, 15, 15, 20 ],
    X0=nothing,
    tolgradnorm=0.1,
    Gstep_only=false,
    GCV=false)
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

    vn() = Vector{Any}(nothing, CDiter)

    output = Dict("SAMm" => vn(), "SAMm_2m" => vn(), "SRE" => vn(), "GVCscore" => vn(),
        "ERGAS_20m" => vn(), "ERGAS_60m" => vn(), "SSIM" => vn(), "aSSIM" => vn(),
        "RMSE" => vn(), "Time" => vn() )
    FFTW.set_num_threads(6)

    tic = time()

    if length(q) != r
        println("The length of q has to match r")
    end

    q = q[:]

    nl,nc = size(Yim[2]);
    L = length(Yim)
    n = nl*nc;


    Yim2, av = normaliseData(Yim)

    #subsampling factors (in pixels)
    d = vec([6 1 1 1 2 2 2 1 2 6 2 2])#'

    # convolution  operators (Gaussian convolution filters), taken from ref [5] from the Ulfarsson
    mtf = vec([ .32 .26 .28 .24 .38 .34 .34 .26 .33 .26 .22 .23])

    sdf = d.*sqrt.(-2*log.(mtf)/(π^2))

    sdf[d .== 1] .= 0

    limsub=2
    dx=12
    dy=12
    FBM, fbm = createConvKernel(sdf,d,nl,nc,L,dx,dy)
    Y,M,F = initialization(Yim2,sdf,nl,nc,L,dx,dy,d,limsub,r)

    Mask = collect(reshape(M,(n,L))')

    if isnothing(X0)
        Z = zeros(r,n)
    else
        #low rank approx
        X0, = normaliseData(X0)
        X0 = reshape(X0, (n,L))'
        U,s,V = svd(X0, full=false)
        Σ = Diagonal(s)
        U = U[:,1:r]
        Z = Σ[1:r,1:r]*V[:,1:r]'
    end

    FDH,FDV,FDHC,FDVC = createDiffkernels(nl,nc,r)
    sigmas = 1.0
    W = collect(computeWeights(Y,d,sigmas,nl))
    Whalf = sqrt.(W)

    if GCV==1
        Gstep_only=1
    end

    if Gstep_only ≈ 1
        CDiter=1
    end

    init_matlab(Y,FBM,nl,nc,Mask)

    Jcost = zeros(CDiter)
    for jCD in 1:CDiter
        @show jCD

        Z, Jcost[jCD], options = Zstep(Y,FBM,fbm,F,lambda,nl,nc,Z,Mask,q,FDH,FDV,FDHC,FDVC,W,Whalf,tolgradnorm)

        if(Gstep_only==0)
           F1=Fstep_matlab(F,Z,Y,FBM,nl,nc,Mask)
           F=F1
        end

        if( GCV==1 )
            Ynoise = ( abs(Y) > 0 ) .* randn( size(Y) )

            Znoise = Zstep(Ynoise,FBM,fbm,F,lambda,nl,nc,Z,Mask,q,FDH,FDV,FDHC,FDVC,W,Whalf,tolgradnorm)

            HtHBXnoise = Mask.*ConvCM(F*Znoise,FBM,nl)

            Ynoise = Ynoise[[1,5,6,7,9,10,11,12],:]
            HtHBXnoise = HtHBXnoise[[1,5,6,7,9,10,11,12],:]


            den = tr(Ynoise*(Ynoise - HtHBXnoise)')
            HtHBX=Mask.*ConvCM(F*Z,FBM,nl)

            num = norm(AY[[1,5,6,7,9,10,11,12],:] - HtHBX[[1,5,6,7,9,10,11,12],:] , p=2).^2;

            output["GCVscore"] = [num / den]
        end

        output["Time"][jCD] = time() - tic

        if !isnothing(Xm_im)
            Xhat_im = collect(conv2im(F*Z,nl,nc,L))
            output["SAMm"][jCD], output["SAMm_2m"][jCD], output["SRE"][jCD],
                 output["RMSE"][jCD], output["SSIM"][jCD], output["aSSIM"][jCD],
                 output["ERGAS_20m"][jCD], output["ERGAS_60m"][jCD] =
                 evaluate_performance(Xm_im,Xhat_im,nl,nc,L,limsub,d,av);
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

function initialization(Yim2::Array{Any,1}, sdf::Array{Float64,1},
    nl::Int,nc::Int,L::Int,
    dx::Int, dy::Int, d::Array{Int,1}, limsub::Int,r::Int)

    FBM2 = createConvKernelSubspace(sdf,nl,nc,L,dx,dy)

    Ylim = zeros(nl, nc, L)

    for i = 1:L
        Ylim[:,:,i] = imresize(Yim2[i], nl, nc)
    end

    Y2im = real.(ifft2(fft2(Ylim).*FBM2))
    Y2tr = Y2im[limsub+1:end-limsub,limsub+1:end-limsub,:]
    Y2n = reshape(Y2tr,(nl-4)*(nc-4),L)
    F, dd, P = svd(Y2n')#,"econ")
    F = F[:,1:r]
    M, Y = createSubsampling(Yim2, d, nl, nc, L)

    return Y, M, F
end

function createSubsampling(Yim::Array{Any,1}, d::Array{Int,1}, nl::Int, nc::Int, L::Int)
    M = zeros(nl, nc, L)
    indexes = Array{Any,1}(undef,L)
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
    X = conv2mat(real.(ifft2(fft2(conv2im(X,nl,nc,L)).*FKM)))
    return X
end


function conv_sp(X::Array{Float64,2}, fkm::Vector{Any}, nl::Int, flip::Bool=true)

    Xim = conv2im(X,nl)
    L = size(Xim,3)
    Xfilt = Array{Float64,3}(undef, size(Xim))
    for ii = 1:L
        Xi = view(Xim,:,:,ii)
        Xfilti = view(Xfilt,:,:,ii)
        ker = length(fkm) == 1 ? fkm[1] : fkm[ii]
        if flip
            ker = reflect(ker)
        end
        imfilter!(Xfilti,Xi,ker,"circular")
    end

    return conv2mat(Xfilt)
end
function diffmat(X::Array{Float64,2}, dim::Int, nl::Int; flip::Bool = false)
        Xim = conv2im(X,nl)

        nl, nc, L = size(Xim)

        if flip
            Xim = reverse(Xim,dims=dim)
        end
        #pad for circular
        if dim == 1
            Xim = vcat(Xim,reshape(Xim[1,:,:],(1,nc,L)))
        else
            Xim = hcat(Xim,reshape(Xim[:,1,:],(nl,1,L)))
        end
        Xfilt = diff(Xim,dims=dim)
        if flip
            Xfilt = reverse(Xfilt,dims=dim)
        end
        return conv2mat(Xfilt)
end

function conv2mat(X,nl::Int=0,nc::Int=0,L::Int=0)
    if ndims(X) == 3
        nl,nc,L = size(X);
        X = reshape(X,nl*nc,L)';
    elseif ndims(X) == 2
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
    fbm::Vector{Any},
    Mask::Array{Float64,2},
    nl::Int, nc::Int, r::Int,
    tau::Float64, q::Array{Float64,1},
    FDH::Array{ComplexF64,3},
    FDV::Array{ComplexF64,3},
    FDHC::Array{ComplexF64,3},
    FDVC::Array{ComplexF64,3},
    W::Array{Float64,2})

    FFTCONV = false

    X = F*Z
    if FFTCONV
        BX = ConvCM(X,FBM,nl)
    else
        BX = conv_sp(X,fbm,nl);
    end
    HtHBX=Mask.*BX;
    if FFTCONV
        ZH=ConvCM(Z,FDHC,nl);
        Zv=ConvCM(Z,FDVC,nl);
    else
        ZH = diffmat(Z,2,nl)
        Zv = diffmat(Z,1,nl)
    end
    ZHW=ZH.*W;
    ZVW=Zv.*W;
    if FFTCONV
        grad_pen=ConvCM(ZHW,FDH,nl)+ConvCM(ZVW,FDV,nl);
        AtAg = F'*ConvCM(HtHBX,conj.(FBM),nl)+2*tau*(q*ones(1,nl*nc)).*grad_pen;
    else
        grad_pen = diffmat(ZHW,2,nl,flip=true) + diffmat(ZVW,1,nl,flip=true)
        AtAg = F'*conv_sp(HtHBX,fbm,nl,false)+2*tau*(q*ones(1,nl*nc)).*grad_pen;
    end
    gradJ=AtAg-UBTMTy;
    J = 1/2 * sum( Z .* AtAg ) - sum( Z.*UBTMTy );

    return J, gradJ, AtAg
end

function CG(Z::Array{Float64,2},
    F::Array{Float64,2},
    Y::Array{Float64,2},
    UBTMTy::Array{Float64,2},
    FBM::Array{ComplexF64,3},
    fbm::Vector{Any},
    Mask::Array{Float64,2},
    nl::Int, nc::Int, r::Int,
    tau::Float64, q::Array{Float64,1},
    FDH::Array{ComplexF64,3},
    FDV::Array{ComplexF64,3},
    FDHC::Array{ComplexF64,3},
    FDVC::Array{ComplexF64,3},
    W::Array{Float64,2})

    maxiter = 1000;
    tolgradnorm = 0.1;  #%1e-6;   %0.1
    cost,grad = grad_cost_G(Z,F,Y,UBTMTy,FBM,fbm,Mask,nl,nc,r,tau,q,FDH,FDV,FDHC,FDVC,W);
    gradnorm = norm(grad[:]);
    iter = 0;
    res = -grad;
    old_res = nothing
    desc_dir = nothing
    while ( gradnorm > tolgradnorm && iter < maxiter )
        iter = iter + 1;
        @printf "%5d\t%+.16e\t%.8e\n" iter cost gradnorm;
        if( iter == 1 )
            desc_dir = res;
        else
            beta = ( res[:]' * res[:] ) / ( old_res[:]' * old_res[:] );
            desc_dir = res + beta * desc_dir;
        end
        cost, _, AtAp = grad_cost_G(desc_dir,F,Y,UBTMTy,FBM,fbm,Mask,nl,nc,r,tau,q,FDH,FDV,FDHC,FDVC,W);
        alpha = ( res[:]' * res[:] ) / ( desc_dir[:]' * AtAp[:] );
        Z1 = Z + alpha * desc_dir;
        old_res = res;
        res = res - alpha* AtAp;
        gradnorm = norm( res[:] );
        # Transfer iterate info
        Z = Z1;
    end

    return Z
end


function computeWeights(Y::Array{Float64,2}, d::Array{Int,1}, sigmas::Float64, nl::Int, method::String="")
    hr_bands = d.==1;
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


function createConvKernel(sdf,d,nl,nc,L,dx,dy)

    middlel = nl÷2
    middlec = nc÷2

    B = zeros(nl,nc,L)
    FBM = zeros(ComplexF64,nl,nc,L)
    fbm = []

    for i in 1:L
        if d[i] >1
            h = _matlab_imagetoolbox_gaussian((dx,dy),sdf[i])

            h /= sum(h)
            #todo - this may be flipped and/or off by one
            push!(fbm,OffsetArray(h,
                    (-dy÷2+1:dy÷2) .- (d[i]÷2) .+ 0,
                    (-dx÷2+1:dx÷2) .- (d[i]÷2) .+ 0))


            B[(middlel-dy÷2+1:middlel+dy÷2) .- (d[i]÷2) .+ 1 , (middlec-dx÷2+1:middlec+dx÷2) .- (d[i]÷2) .+ 1, i] = h #not sure about this line here

            B[:,:,i]=fftshift( B[:,:,i] )#/sum(sum(B[:,:,i]))
            FBM[:,:,i] = fft2(B[:,:,i])
        else
            push!(fbm,centered(ones(1,1)))
            B[1,1,i]=1
            FBM[:,:,i]=fft2(B[:,:,i])
        end
    end


    return FBM, fbm
end



function createConvKernelSubspace(sdf,nl,nc,L,dx,dy)
    middlel=Int(round((nl+1)/2,RoundNearestTiesUp))
    middlec=Int(round((nc+1)/2,RoundNearestTiesUp))

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

            σ = sqrt(s2^2 - sdf[i]^2)
            h = Kernel.gaussian((σ,σ),(dx,dy))
            #h = Kernel.gaussian(3) #this will produce a 13x13 kernel which is close the example 12x12 needed

            B[middlel-(dy-1)÷2:middlel+(dy-1)÷2,middlec-(dx-1)÷2:middlec+(dx-1)÷2,i] = h;
            B[:,:,i]=fftshift( B[:,:,i] )/sum(sum(B[:,:,i]))
            FBM2[:,:,i] = fft2(B[:,:,i])
        else
            B[1,1,i]=1
            FBM2[:,:,i]=fft2(B[:,:,i])
        end

    end

    return FBM2

end



function Zstep(
    Y::Array{Float64,2},
    FBM::Array{ComplexF64,3},
    fbm::Vector{Any},
    F::Array{Float64,2},
    tau::Float64, nl::Int, nc::Int,
    Z::Array{Float64,2},
    Mask::Array{Float64,2},
    q::Array{Float64,1},
    FDH::Array{ComplexF64,3},
    FDV::Array{ComplexF64,3},
    FDHC::Array{ComplexF64,3},
    FDVC::Array{ComplexF64,3},
    W::Array{Float64,2},
    Whalf::Array{Float64,2},
    tolgradnorm::Float64)


    r = size(F,2);
    n = nl*nc;
    UBTMTy=F'*ConvCM(Y,conj.(FBM),nl);
    Z = CG(Z,F,Y,UBTMTy,FBM,fbm,Mask,nl,nc,r,tau,q,FDH,FDV,FDHC,FDVC,W);
    xcost=1;
    options=[];

    return Z, xcost, options

end


function init_matlab(Y,FBM,nl,nc,Mask)
    mat"cd('C:/Users/Konrad/Documents/EECS556/proj/eecs556-sharp/julia')"
    mxcall(:matlabInit,0)
    mat"Y = $Y;"
    mat"FBM = $FBM;"
    mat"Mask= $Mask;"
    mat"nl= $nl;"
    mat"nc= $nc;"
    println("Successfully initialized MATLAB environment")
end
function Fstep_matlab(F::Array{Float64,2},
    Z::Array{Float64,2},
    Y::Array{Float64,2},
    FBM::Array{ComplexF64,3},
    nl::Int,nc::Int,
    Mask::Array{Float64,2})
    mat"$F1 = Fstep($F,$Z,Y,FBM,nl,nc,Mask)"
    return F1
    #return mxcall(:Fstep,1,F,Z,Y,FBM,nl,nc,Mask)
end

# function Fstep(F::Array{Float64,2},
#         Z::Array{Float64,2},
#         Y::Array{Float64,2},
#         FBM::Array{ComplexF64,3},
#         nl::Int,nc::Int,
#         Mask::Array{Float64,2})
#      F0=F;#%   U; % initialization
#      BTXhat =  ConvCM(F0*Z,FBM,nl);
#      MBTXhat=Mask.*BTXhat;
#      L, r = size(F);
#      n = nl*nc
#      MBZT = zeros(r,n,L)
#      A = zeros(r,r,L)
#      ZBMTy = zeros(r,L)
#      for ii=1:L
#          #@show size(Mask)
#          #@show size(repeat(Mask[ii,:]',outer=(r,1)))
#          #@show size(FBM)
#          #@show size(Z)
#          #@show size(ConvCM(Z,repeat(FBM[:,:,ii],outer=(1,1,r)),nl))
#         MBZT[:,:,ii] = repeat(Mask[ii,:]',outer=(r,1)).*
#                         ConvCM(Z,repeat(FBM[:,:,ii],outer=(1,1,r)),nl);
#         A[:,:,ii]=MBZT[:,:,ii]*MBZT[:,:,ii]';
#         ZBMTy[:,ii]=MBZT[:,:,ii]*Y[ii,:];
#      end
#      ZBYT=ZBMTy';#   BTY*Z';
#
#      #You may need to install Manifolds package for this to work
#      manifold = Manopt.Stiefel(L,r)
#      cost(_, F) = costF(F, MBZT, Y)
#      grad(_, F) = egrad(F, A, ZBYT)
#      hessam = Manopt.ApproxHessianFiniteDifference(manifold, F0, grad)
#
#
#      F1 = Manopt.trust_regions(manifold, cost, grad, hessam, F0;
#             stopping_criterion = StopWhenGradientNormLess(1e-2))
#
#     return F1
#
#      #manifold = stiefelfactory(L,r,1);# %euclideanfactory(L,r);
#      #problem.M = manifold;
#      #problem.cost  = @(F) costF(F,MBZT,Y);
#      #problem.egrad = @(F) egrad(F,A,ZBYT);
#      #warning('off', 'manopt:getHessian:approx')
#      #options.tolgradnorm = 1e-2;
#      #options.verbosity=0;
#      #[F1, xcost, info, options] = trustregions(problem,F0,options);
# end

function ERGAS(I1::Array{Float64,3},
    I2::Array{Float64,3},
    ratio::Int)

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
        d::Array{Int,1},
        av::Array{Float64,1})

    Xhat_im = Xhat_im[limsub+1:end-limsub,limsub+1:end-limsub,:];
    Xhat_im = unnormaliseData(Xhat_im,av);
    Xhat=reshape(Xhat_im,((nl-4)*(nc-4),L));
    #% Xm_im is the ground truth image
    Xm_im = Xm_im[limsub+1:end-limsub,limsub+1:end-limsub,:];
    if ( size(Xm_im,3) == 6 ) #% Reduced Resolution
        ind = findall( d==2 );
        SAMm=SAM(Xm_im,Xhat_im[:,:,ind]);
        SAMm_2m=SAMm;
        X = conv2mat(Xm_im);
        Xhat = conv2mat(Xhat_im);
        SRE = zeros(6,1)
        SSIM_index = zeros(6,1)
        #% SRE - signal to reconstrution error
        for i=1:6
            SRE[i,1] = 10*log10(sum(X[i,:].^2)/ sum((Xhat(ind(i),:)-X[i,:]).^2));
            SSIM_index[i,1] = assess_ssim(Xm_im[:,:,i],Xhat_im(:,:,ind(i)));
        end
        aSSIM=mean(SSIM_index);
        ERGAS_20m = ERGAS(Xm_im,Xhat_im[:,:,ind],2);
        ERGAS_60m = NaN;
        RMSE = norm(X - Xhat[ind,:]) / size(X,2);
    else
        ind = findall((d .==2) .| (d .== 6))
        #@show ind
        #ind = [ii[2] for ii in ind]#?????
        SAMm=SAM(Xm_im[:,:,ind],Xhat_im[:,:,ind]);
        ind2=findall(d .== 2);
        SAMm_2m=SAM(Xm_im[:,:,ind2],Xhat_im[:,:,ind2]);
        ind6=findall(d .== 6);
        X = conv2mat(Xm_im);
        Xhat = conv2mat(Xhat_im);
        #% SRE - signal to reconstrution error
        SRE = zeros(L,1)
        SSIM_index = zeros(L,1)
        for i=1:L
            SRE[i,1] = 10*log10(sum(X[i,:].^2)/ sum((Xhat[i,:]-X[i,:]).^2));
            SSIM_index[i,1] = assess_ssim(Xm_im[:,:,i],Xhat_im[:,:,i]);
        end
        aSSIM=mean(SSIM_index);
        ERGAS_20m = ERGAS(Xm_im[:,:,ind],Xhat_im[:,:,ind],2);
        ERGAS_60m = ERGAS(Xm_im[:,:,ind2],Xhat_im[:,:,ind2],6);
        RMSE = norm(X[ind,:] .- Xhat[ind,:]) / size(X,2);
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

    fft_dh = fft2(dh)
    fft_dv = fft2(dv)

    FDH = zeros(ComplexF64,nl,nc,r)
    FDV = zeros(ComplexF64,nl,nc,r)

    for i in 1:r
        FDH[:,:,i] = fft_dh
        FDV[:,:,i] = fft_dv
    end

    FDHC = conj(FDH)
    FDVC = conj(FDV)

    return FDH,FDV,FDHC,FDVC, dh, dv
end


# julia ImageFiltering.kernel.Gaussian does not support even kernel lengths,
# which the MATLAB version of S2Sharp used. this function was adapted from
# fspecial in MATLAB Image Toolbox.
function _matlab_imagetoolbox_gaussian(size::NTuple{2,Int},σ::AbstractFloat)
    siz   = (size.-1)./2;
    std   = σ;

    x = -siz[2] : siz[2]
    y = transpose(-siz[1] : siz[1])

    arg   = -(x.*x .+ y.*y)./(2*std*std);

    h     = exp.(arg);
    h[h.<eps()*maximum(h[:])] .= 0;

    sumh = sum(h[:]);
    if sumh != 0
        h  = h/sumh;
    end

    return h
end

function normaliseData(Yim::Array{Any,1})
    #Yim is the array that has all the different images in it

    nb = length(Yim)
    av = zeros(nb)
    Yim = deepcopy(Yim)

    for i in 1:nb
        av[i] = mean(mean(Yim[i].^2))
        Yim[i] = sqrt.(Yim[i].^2/av[i])
    end

    return Yim, av
end


function unnormaliseData(Xim::Array{Float64,3},av::Array{Float64,1})

    Xim2 = zeros(size(Xim))
    nb = size(Xim,3)
    for i in 1:nb
        Xim2[:,:,i] = sqrt.(Xim[:,:,i].^2 * av[i])
    end
    return Xim2
end
function unnormaliseData(Yim::Array{Any, 1},av::Array{Float64,1})

    nb, =size(Yim)

    Yim2 = []

    for i in 1:nb
        append!(Yim2,sqrt.(Yim[i].^2 * av[i]))
    end
    return Yim2
end



function conv2im(X,nl,nc=0,L=0)

    if ndims(X) == 1
        X = transpose(X)
        #X = conv2mat(X,nl,nc,L)
    end

    if L==0 && nc==0 #this is the easiest sub for L,nc not being included
        L, n = size(X)
        if n == 1
            X = conv2mat(X,nl,nc,L)
        end
        nc = n÷nl
    end

    return reshape(X',(nl,nc,L))
end

function SAM(I1,I2)

    M, N, _ = size(I2);
    dot3(A,B) = sum(conj(A).*B,dims=3)
    prod_scal = dot3(I1,I2);
    norm_orig = dot3(I1,I1);
    norm_fusa = dot3(I2,I2);
    prod_norm = sqrt.(norm_orig.*norm_fusa);
    prod_map = prod_norm;
    prod_map[prod_map .== 0] .= eps();
    SAM_map = acos.(prod_scal./prod_map);
    prod_scal = reshape(prod_scal, M*N,1);
    prod_norm = reshape(prod_norm, M*N,1);
    nz = prod_norm .!= 0;
    prod_scal = prod_scal[nz]
    prod_norm = prod_norm[nz]
    angolo = sum(sum(acos.(prod_scal./prod_norm)))/(size(prod_norm,1));
    SAM_index = real.(angolo)*180/pi;

    return SAM_index#, SAM_map
end

#FFTs along first 2 dimensions - convenience functions to match matlab
fft2(X) = fft(X,(1,2))
ifft2(X) = ifft(X,(1,2))
