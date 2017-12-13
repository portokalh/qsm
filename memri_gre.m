%MEMRI GRE
%using Russel Recon

addpath(genpath('/Volumes/omega/Natalie/Russell_matlab/'));

cd '/Users/omega/Natalie/youngmice/MGE/workfiles/B03681.work'

myimg=RussReconADMask('bruker','center','fermi','save');

phimask=open_nii('maskME.nii.gz');
FreqC=open_nii('FreqC.nii.gz');
load('B0dir.mat');
B0 = 7.0;
load TE
load meanSNR


T2 = T2mapMaskedWeighted(myimg,TE,phimask(:,:,:,1),meanSNR);
XCFreqC = QSM_iLSQR(FreqC*2*pi/1e3,phimask(:,:,:,1),'B0',B0,'H',B0dir);

show3(XCFreqC);
show3(T2);

T22=abs(T2);
myT2nii=make_nii(T22,[0.1,0.1,0.1],[0, 0 ,0],16);
save_nii(myT2nii, ['B03681' 'T2.nii.gz']);

XCFreqC1=abs(XCFreqC);
myfreqnii=make_nii(XCFreqC1,[0.1,0.1,0.1],[0, 0, 0],16);
save_nii(myfreqnii, ['B03681' Chi.nii.gz']);