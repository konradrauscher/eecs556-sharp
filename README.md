# eecs556-sharp

Here you can find the code used for the course project for EECS556 Image Processing (team members: Henry Do, John Gearig, Anusha Kikkeri and Konrad Rauscher).

To run the Julia code conversion, you'd need to download example.jl, S2sharp.jl, Fstep.m, and matlabinit.m. In order to do the Fstep, Julia will call a function that communicates to Matlab in order to use the Manopt library for Matlab. You will likely need to change file paths in matlabinit.m-- this was working on a Windows operating system and may encounter issues on Mac. 

To run the Matlab code, follow the path eecs556-sharp/matlab/S2Sharp/ and download example.m, manopt folder, data folder, and S2sharp.m. Running the example file will reproduce some of the pansharpened plots seen in our report. read_s2_data.m is a useful script to get any S2 jp2 image from our online sources into the Yim format that S2Sharp function requires to run.  

Using the same dataset, with manipulation of the input data, one can create pan-sharpened outputs using AWLP, MTF-GLP, and GS methods, located under the matlab folder. 
The base code is from: https://openremotesensing.net/knowledgebase/a-new-benchmark-based-on-recent-advances-in-multispectral-pansharpening-revisiting-pansharpening-with-classical-and-emerging-pansharpening-methods/

In order to create reference images when using real Sentinel-2 data, the downsampling.jl code can be used to filter and downsampling the original images to get lower-resolution images for testing.

Also in the eecs556-sharp/matlab folder is the pansharpening toolbox, which was used to compare our method against. 
