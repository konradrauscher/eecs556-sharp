using LinearAlgebra
using ImageTransformations
using FFTW
using ImageFiltering
using Statistics

function initialization(Yim2::Array{Any}, sdf::Array{Float64},
    nl::Int,nc::Int,L::Int,
    dx::Int, dy::Int, d::Array{Int}, limsub::Int,r::Int)

    FBM2 = Nothing
    #FBM2 = createConvKernelSubspace(sdf,nl,nc,L,dx,dy)

    Ylim = zeros(nl, nc, L)

    for i = 1:L
        Ylim[:,:,i] = imresize(Yim2[i],d[i])
    end

    Y2im = real.(ifft2(fft2(Ylim).*FBM2))
    Y2tr = Y2im[limsub+1:end-limsub,limsub+1:end-limsub,:]
    Y2n = reshape(Y2tr,[(nl-4)*(nc-4),L])
    F, d, P = svd(Y2n')#,"econ")
    F = F(:,1:r)
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
    indexes[i] = findall(M[:,:,i] == 1)
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


function conv2mat(X,nl::Int,nc::Int,L::Int)
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


function computeWeights(Y::Array{Float64,2}, d::Array{Int}, sigmas::Float64, nl::Int, method::String)
    hr_bands = d==1;
    hr_bands = findall(hr_bands)';
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
    grad = sqrt(max(grad,dims=3));
    grad = grad / quantile(grad(:),0.95);

    Wim = exp.(-grad.^2/2/sigmas^2);
    Wim[Wim.<0.5] = 0.5;

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
