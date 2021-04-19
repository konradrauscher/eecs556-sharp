BASE = 'C:\Users\Konrad\Documents\EECS556\proj\S2A_MSIL2A_20200926T162041_N0214_R040_T16TGM_20200926T205347.SAFE\GRANULE\L2A_T16TGM_A027493_20200926T162748\IMG_DATA\';

Xm_im = [];

tl_fres = [901 6751];
br_fres = [1500 7350];
tl_20m = (tl_fres-1)/2+1;
br_20m = br_fres/2;
tl_60m = (tl_fres-1)/6+1;
br_60m = br_fres/6;

PR = 'PixelRegion';
pr10 = {[tl_fres(1) br_fres(1)], [tl_fres(2) br_fres(2)]};
pr20 = {[tl_20m(1) br_20m(1)], [tl_20m(2) br_20m(2)]};
pr60 = {[tl_60m(1) br_60m(1)], [tl_60m(2) br_60m(2)]};

Yim = {};
Yim{end+1} = imread([BASE '\R60m\T16TGM_20200926T162041_B01_60m.jp2'],PR,pr60);
Yim{end+1} = imread([BASE '\R10m\T16TGM_20200926T162041_B02_10m.jp2'],PR,pr10);
Yim{end+1} = imread([BASE '\R10m\T16TGM_20200926T162041_B03_10m.jp2'],PR,pr10);
Yim{end+1} = imread([BASE '\R10m\T16TGM_20200926T162041_B04_10m.jp2'],PR,pr10);
Yim{end+1} = imread([BASE '\R20m\T16TGM_20200926T162041_B05_20m.jp2'],PR,pr20);
Yim{end+1} = imread([BASE '\R20m\T16TGM_20200926T162041_B06_20m.jp2'],PR,pr20);
Yim{end+1} = imread([BASE '\R20m\T16TGM_20200926T162041_B07_20m.jp2'],PR,pr20);
Yim{end+1} = imread([BASE '\R10m\T16TGM_20200926T162041_B08_10m.jp2'],PR,pr10);
Yim{end+1} = imread([BASE '\R20m\T16TGM_20200926T162041_B8A_20m.jp2'],PR,pr20);
Yim{end+1} = imread([BASE '\R60m\T16TGM_20200926T162041_B09_60m.jp2'],PR,pr60);
%Yim{end+1} = imread([BASE '\R60m\T16TGM_20200926T162041_B10_60m.jp2'],PR,pr60);
Yim{end+1} = imread([BASE '\R20m\T16TGM_20200926T162041_B11_20m.jp2'],PR,pr20);
Yim{end+1} = imread([BASE '\R20m\T16TGM_20200926T162041_B12_20m.jp2'],PR,pr20);

RGB = cat(3,Yim{4},Yim{3},Yim{2});
figure; imshow(RGB*2^5);
