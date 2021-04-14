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
        
    for i in 0:2:(length(varargin)-1)
        current = varargin[i]
        next = varargin[i+1]
        
        if current == "CDiter"
            CDiter = next
        elseif current == "r"
            r = next
        elseif current == "lamda"
            lamda = next
        elseif current == "XM_im"
            XM_im = next
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
            println("Bad varargin INPUT")
            
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
    Yim = vec(Yim)
    #may need to convert data type here
    
    ##can't visualize what this does
    
#     [nl,nc] = size(Yim{2});
#     n = nl*nc;
    
    
    Yim2, av = normalizeData(Yim)
    
    #subsampling factors (in pixels)
    d = [6 1 1 1 2 2 2 1 2 6 2 2]'
    
    # convolution  operators (Gaussian convolution filters), taken from ref [5] from the Ulfarsson
    mtf = [ .32 .26 .28 .24 .38 .34 .34 .26 .33 .26 .22 .23]
    
    sdf = d.*sqrt.(-2*log.(mtf)/(π^2))'
    
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
    Whalf = sqrt.(W)
    
    if GCV==1
        Gstep_only=1
    end
    
    if Gstep_only ≈ 1
        CDiter=1
    end
    
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
