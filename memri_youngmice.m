%MEMRI GRE
%using Russel Recon

addpath(genpath('/Users/omega/Documents/Natalie/SCCAN_tutorial/Russell_matlab/'));

%runnos = {'B03618' 'B03681' 'B03705' 'B03710' 'B03715' 'B03720' 'B03725' 'B03730' 'B03735' 'B03740' 'B03819' 'B03824' 'B03829' 'B03835' 'B03853' 'B03859' 'B03865' 'B03871' 'B04007' 'B04012' 'B04017' 'B04252' 'B04256' 'B04262'};
runnos = {'B03618', 'B03681' 'B03705' 'B03710' 'B03715' 'B03720' 'B03725' 'B03730' 'B03735'  'B03740'  'B03819' 'B03824' 'B03829'  'B03835' 'B03853' 'B03859' 'B03865' 'B03871' 'B04007' 'B04012'  'B04017' 'B04252' 'B04256' 'B04262' } ;

flag_mask=2;
%:numel(runnos)

% for 
i=24;
runno=runnos{i};

mypath=['/Users/omega/Documents/Natalie/youngmice/MGE/workfiles/' char(runno) '.work'];
cd(mypath)

myimg=RussReconADMask_SE('bruker','fermi','save'); % Took off center 

% phimask=open_nii('maskME.nii.gz'); %auto generates

myGREnii=make_nii(myimg,[0.1,0.1,0.1],[0, 0 ,0],16);
save_nii(myGREnii, [char(runno) 'mGRE.nii.gz']);
%%
load('B0dir.mat');
B0 = 7.0;
load TE
load meanSNR
% XCFreqC = QSM_iLSQR(FreqC*2*pi/1e3,phimask1,'B0',B0,'H',B0dir);
Freq = open_nii('Freq.nii.gz');
star = zeros(size(Freq));
mask = open_nii('B04262_masked.nii.gz');
voxelsize=[.1 .1 .1];
for cpt = 1 %size(Freq,4)
star(:,:,:,cpt) = QSM_star(2*pi*Freq(:,:,:,cpt)*TE(1,cpt),single(mask), 'H',[1 0 0], 'voxelsize' ,voxelsize,'padsize',[12 12 12],'TE',TE(1,cpt)*1000,'B0',7);
end
% qsm = (mean(star(:,:,:,2:6),4));
star(star == 0) = NaN;
qsm = (nanmean(star(:,:,:,1),4));
[W,qsm2] = nanMEW(TE,.025,star);
%show3(qsm)
save_nii(make_nii(qsm,voxelsize,[0 0 0]), 'B04262_QSM.nii.gz');
%end
