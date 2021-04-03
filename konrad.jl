using LinearAlgebra
using ImageTransformations
using FFTW

function initialization(Yim2::Array{Any}, sdf::Array{Float},
    nl::Int,nc::Int,L::Int,
    dx::Int, dy::Int, d::Array{Int}, limsub::Int,r::Int)

    FBM2 = Nothing
    #FBM2 = createConvKernelSubspace(sdf,nl,nc,L,dx,dy)

    Ylim = zeros(nl, nc, L)

    for i = 1:L
        Ylim[:,:,i] = imresize(Yim2[i],d[i])
    end

    Y2im = real.(ifft2(fft2(Ylim).*FBM2))
    Y2tr = Y2im(limsub+1:end-limsub,limsub+1:end-limsub,:)
    Y2n = reshape(Y2tr,[(nl-4)*(nc-4),L])
    F, d, P = svd(Y2n')#,"econ")
    F = F(:,1:r)
    M, Y = createSubsampling(Yim2, d, nl, nc, L)

end

function createSubsampling(Yim::Array{Any}, d::Array{Int}, nl::Int, nc::Int, L::Int)
M = zeros(nl, nc, L)
indexes = Array{Any}(undef,L)

for i = 1:L
    im = ones(nl÷d[i],nc÷d[i])
    maux = zeros(d[i],d[i])
    maux[1,1] = 1
    M[:,:,i] = kron(im,maux)
    indexes[i] = find(M[:,:,i] == 1)
    Y[i,indexes[i]] = conv2mat(Yim[i],nl÷d[i],nc÷d[i],1)
end

end

floori(x) = Int(floor(x))

function conv2mat(X,nl::Int,nc::Int,L::Int)
    if ndims(X) == 3
        nl,nc,L = size(X);
        X = reshape(X,nl*nc,L)';
    elseif ndims(squeeze(X)) == 2
        L = 1;
        nl,nc = size(X);
        X = reshape(X,nl*nc,L)';
    end
end
