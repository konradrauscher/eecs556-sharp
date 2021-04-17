function F1=Fstep(F,Z,Y,FBM,nl,nc,Mask)

     [L,r]=size(F);
     F0=F;%   U; % initialization
     BTXhat =  ConvCM(F0*Z,FBM,nl);
     MBTXhat=Mask.*BTXhat;
     for ii=1:L
       MBZT(:,:,ii)=repmat(Mask(ii,:),[r,1]).*ConvCM(Z,repmat(FBM(:,:,ii),[1,1,r]),nl);
       A(:,:,ii)=MBZT(:,:,ii)*MBZT(:,:,ii)';
       ZBMTy(:,ii)=MBZT(:,:,ii)*Y(ii,:)';
     end
     ZBYT=ZBMTy';%    BTY*Z';
     manifold = stiefelfactory(L,r,1); %euclideanfactory(L,r); 
     problem.M = manifold;
     problem.cost  = @(F) costF(F,MBZT,Y); 
     problem.egrad = @(F) egrad(F,A,ZBYT);  
     warning('off', 'manopt:getHessian:approx') 
     options.tolgradnorm = 1e-2;
     options.verbosity=0;
     [F1, xcost, info, options] = trustregions(problem,F0,options);
end

function [Ju]=costF(F,MBZT,Y)
    L=size(F,1);
    Ju=0;
    for i=1:L,
        fi=F(i,:)';
        yi=Y(i,:)';
        Ju=Ju+0.5*norm(MBZT(:,:,i)'*fi-yi,'fro')^2;
    end
end

function [Du]=egrad(F,A,ZBYT)
    p=size(A,3);
    Du=0*F;
    for ii=1:p
        Du(ii,:)=F(ii,:)*A(:,:,ii)'-ZBYT(ii,:);
    end
end


function X = ConvCM(X,FKM,nl,nc,L)

    if nargin == 3
        [L,n] = size(X);
        nc = n/nl;
    end
    X = conv2mat(real(ifft2(fft2(conv2im(X,nl,nc,L)).*FKM)));
end
function X = conv2im(X,nl,nc,L)

    if size(X,2)==1
        X = conv2mat(X,nl,nc,L);
    end
    if nargin == 2
        [L,n] = size(X);
        if n==1
            X = conv2mat(X,nl,nc,L);
        end
        nc = n/nl;
    end
    X = reshape(X',nl,nc,L);
end

function X = conv2mat(X,nl,nc,L)
    if ndims(X) == 3
        [nl,nc,L] = size(X);
        X = reshape(X,nl*nc,L)';
    elseif ndims(squeeze(X)) == 2
        L = 1;
        [nl,nc] = size(X);
        X = reshape(X,nl*nc,L)';
    end
end