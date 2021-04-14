import FFTW: fft,fftshift
using Statistics: mean
using ImageFiltering
using LinearAlgebra

function createDiffkernels(nl,nc,r)    
    dh = zeros(nl,nc)
    dh[1,1]=1
    dh[1,nc]=-1
    
    dv = zeros(nl,nc)
    dv[1,1]=1
    dv[nl,1]=-1
    
    fft_dh = fft(dh)
    fft_dv = fft(dv)
    
    FDH = zeros(nl,nc,r)
    FDV = zeros(nl,nc,r)

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
    
    nb,=size(Yim)
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
    middlel = round(nl/2)
    middlec = round(nc/2)
    
    B = zeros(nl,nc,L)
    FBM = zeros(nl,nc,L)
    
    for i in 1:L
        if d[i] >1
            f = Kernel.gaussian(3) #this will produce a 13x13 kernel which is close the example 12x12 needed
            
            #this function won't work.. not adaptable enough
            #need to create a 12x12 filter here...
            
            B[ (middlel-dy/2+1:middlel+dy/2)-d(i)/2+1, (middlec-dx/2+1:middlec+dx/2)-d(i)/2+1, i] = h #not sure about this line here
            
            B[:,:,i]=fftshift( B[:,:,i] )/sum(sum(B[:,:,i]))
            FBM[:,:,i] = fft(B[:,:,i])
        else
            B[1,1,i]=1
            FMB[:,:,i]=fft(B[:,:,i])
        end
    end
    
    return FBM
end



function createConvKernelSubspace(sdf,nl,nc,L,dx,dy):    
    middlel=round((nl+1)/2)
    middlec=round((nc+1)/2)

    dx = dx+1
    dy = dy+1

    # kernel filters expanded to size [nl,nc]
    B = zeros(nl,nc,L)
    # fft2 of kernels
    FBM2 = zeros(nl,nc,L)
    s2 = maximum(s2)
    
    ##this loop is nearly identical to whats happening above! this would be a good place to improve it
    for i in 1:L
        if sdf[i] < s2 # != would be a little better i think
            f = Kernel.gaussian(3) #TODO change this 
            
            B[ (middlel-dy/2+1:middlel+dy/2)-d(i)/2+1, (middlec-dx/2+1:middlec+dx/2)-d(i)/2+1, i] = h
            B[:,:,i]=fftshift( B[:,:,i] )/sum(sum(B[:,:,i]))
            FBM[:,:,i] = fft(B[:,:,i])
        else
            B[1,1,i]=1
            FMB[:,:,i]=fft(B[:,:,i])
            
  
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
