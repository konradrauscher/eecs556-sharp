using MAT

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
