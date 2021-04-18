using MAT
#TODO for anyone using this: change the paths here and in matlabInit.m to match
#your matlab installation and folders for this project
ENV["MATLAB_HOME"] = "C:\\Program Files\\MATLAB\\R2020b"
ENV["PATH"] = "C:\\Program Files\\MATLAB\\R2020b\\bin\\win64;" * ENV["PATH"]
include("S2Sharp.jl")

file = matopen("../matlab/S2Sharp/Data/Aviris_cell_3.mat")
Xm_im = read(file, "Xm_im")
Yim = vec(read(file, "Yim"))
close(file)

r = 8;# % subspace dimension / the rank
q  = [1, 0.3851, 6.9039, 19.9581, 47.8967, 27.5518, 2.7100, 34.8689];# % lam_i = lam*qi
ni = 5;#  % number of iterations / one can select fewer iterations (e.g. ni=2) which yields similar result
        #  % but is much quicker
lam = 1.8998e-04;


Xhat_im, output_S2 =S2sharp(Yim,["Xm_im",Xm_im,"r",r,"lambda",lam,"q",q,
            "CDiter",ni]);
##
S2sharp_SRE = output_S2["SRE"][end][[1,5,6,7,9,10,11,12]];
S2sharp_aSRE= mean(S2sharp_SRE);
S2sharp_SAM = output_S2["SAMm"][end];
S2sharp_RMSE = output_S2["RMSE"][end];
S2sharp_aSSIM = output_S2["aSSIM"][end];
#S2sharp_ERGAS_60m=output_S2.ERGAS_60m[end];
#S2sharp_ERGAS_20m=output_S2.ERGAS_20m[end];
S2sharp_time = output_S2["Time"][end];

@show S2sharp_SRE
@show S2sharp_SAM
@show S2sharp_aSRE
@show S2sharp_aSSIM
@show S2sharp_RMSE
@show S2sharp_time
println("B1    B5    B6    B7    B8a   B9    B11   B12")
@printf "%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n\n" S2sharp_SRE...
