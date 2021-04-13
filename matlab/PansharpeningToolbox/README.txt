

Copyright (c) 2020.
All rights reserved. This work should only be used for nonprofit purposes.




MATLAB package for the following pansharpening algorithms:
	EXP: MS image interpolation using a polynomial kernel with 23 coefficients [43]
	BT-H [59]
	BDSD [33]
	C-BDSD [34]
	BDSD-PC [35]
	GS [24]
	GSA [30]
	C-GSA [32]
	PRACS [29]
	AWLP [61] with revised statistical matching between PAN and MS bands [31]
	MTF-GLP: GLP [43] with MTF-matched filter [47], unitary injection model, and revised statistical matching between PAN and MS bands [31]
	MTF-GLP-FS: GLP [43] with MTF-matched filter [47] with a new FS regression-based injection model [133]
	MTF-GLP-HPM: GLP with MTF-matched filter [47] with HPM injection model [143] and revised statistical matching between PAN and MS bands [31]
	MTF-GLP-HPM-H: GLP with MTF-matched filter [47] with HPM injection model [143] and haze correction [59]
	MTF-GLP-HPM-R: GLP with MTF-matched filter [47] and HPM injection model [143] with a new preliminary regression-based spectral matching phase [131]
	MTF-GLP-CBD: GLP [43] with MTF-matched filter [47] and regression-based injection model [132]
	C-MTF-GLP-CBD: context-based GLP [43] with MTF-matched filter [47] and regression-based injection model [132] with local parameter estimation exploiting clustering [32]
	MF: nonlinear decomposition scheme using MFs based on half gradient [55]
	FE-HPM: filter estimation based on a semiblind deconvolution framework and HPM injection model [50]
	SR-D: pansharpening based on sparse representation of injected details [88]
	PWMBF: model-based fusion using PCA and wavelets [79]
	TV: pansharpening based on TV [67]
	RR: model-based reduced-rank pansharpening [72], [134]
	PNN: proposed in [91]
	PNN-IDX: PNN augmenting the input by including several maps of nonlinear radiometric indexes as proposed in [91]
	A-PNN: proposed in [95]
	A-PNN-FT: proposed in [95]
For more information and to have a look at the relative references, please refer to the paper:
G. Vivone, M. Dalla Mura, A. Garzelli, R. Restaino, G. Scarpa, M.O. Ulfarsson, L. Alparone, and J. Chanussot, "A New Benchmark Based on Recent Advances in Multispectral Pansharpening: Revisiting pansharpening with classical and emerging pansharpening methods", IEEE Geoscience and Remote Sensing Magazine, doi: 10.1109/MGRS.2020.3019315.

Main functions to call:  
        "Fusion_Algorithms_Reduced_Resolution", "Fusion_Algorithms_Full_Resolution", and "Main_Reduced_Resolution"

How to use:
	For the full resolution assessment:
		1. Load the full resolution dataset (I_MS_LR is the low resolution MS image; I_PAN is the PAN image; I_MS is the upsampled MS image at PAN scale)
		2. Run Fusion_Algorithms_Full_Resolution.m
	For the reduced resolution assessment:
		1. Load the reduced resolution dataset (I_GT is the ground-truth image (high spatial resolution MS image); I_MS_LR is the low resolution MS image; I_PAN is the PAN image; I_MS is the upsampled MS image at PAN scale)
		NOTE: If the available dataset is at full resolution, you can obtain the reduced resolution version by properly modifying the file Main_Reduced_Resolution.m
		2. Run Fusion_Algorithms_Reduced_Resolution.m

System requirements:
        Matlab2019b or higher versions, with deep learning toolboxes.
        The running on previous Matlab versions is not guaranteed.

Execution Environment:
	The code of CS, MRA, and VO approaches run on CPU.
        The ML code runs on GPU when available or on CPU, otherwise.