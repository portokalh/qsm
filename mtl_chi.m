%MEMRI GRE
%using Russel Recon
addpath(genpath('/Users/omega/Natalie/Russell_matlab'));

% runnos = {'B03618' 'B03681' 'B03705' 'B03710' 'B03715' 'B03720' 'B03725' 'B03730' 'B03735' 'B03740' 'B03819' 'B03824' 'B03829' 'B03835' 'B03853' 'B03859' 'B03865' 'B03871' 'B04007' 'B04012' 'B04017' 'B04252' 'B04256' 'B04262'};
% % runnos={'B04027'
% % 'B04051'
% 'B04077'
% 'B04080'
% 'B04082'
% 'B04085'
% 'B04088'
% 'B04091'
% 'B04094'
% 'B04097'
% 'B04115'
% 'B04118'
% 'B04121'
% 'B04124'
% 'B04127'
% 'B04130'
% 'B04133'
% 'B04030'
% 'B04033'
% 'B04036'
% 'B04039'
% 'B04021'
% 'B04024'
% 'B04041'
% 'B04045'
% 'B04048'};

runnos = {'S67158'};
% flag_mask=1;
%:numel(runnos)

for i=1:length(runnos)
    runno=runnos{i};

mypath=['/Users/omega/Natalie/MTL/' char(runno) '.work'];
cd(mypath)
myimg=RussReconADMask('agilent','fermi','center','save');
phimask=open_nii('maskME.nii.gz'); %auto generates
%manmask_file=['/Users/omega/Natalie/Russell_matlab/LE_masks/' runno '_LEmask.nii'];%manual generated based on RARE aligned to MGE
%manmask=load_nii(manmask_file);
mask = uint8(ones(size(phimask(:,:,:,1))));
% TE=[0.0032 0.0087 0.0142 0.0197 0.0252  0.0307 0.036 0.0417];
TE = [4.24 10.6 16.96 23.32 ]*0.001; %3.44666666666667;
voxelsize=[.1 .1 .1];

% XCFreqC = QSM_iLSQR(FreqC*2*pi/1e3,phimask1,'B0',B0,'H',B0dir);
Freq = open_nii('Freq.nii.gz');
star = zeros(size(Freq));
for cpt = 1:4 %size(Freq,4)
star(:,:,:,cpt) = QSM_star(2*pi*Freq(:,:,:,cpt)*TE(1,cpt),single(mask), 'H',[1 0 0], 'voxelsize' ,voxelsize,'padsize',[12 12 12],'TE',TE(1,cpt)*1000,'B0',7);
end
% qsm = (mean(star(:,:,:,2:6),4));
star(star == 0) = NaN;
qsm = (nanmean(star(:,:,:,1:4),4));
[W,qsm2] = nanMEW(TE,.025,star);
%show3(qsm)
save_nii(make_nii(qsm,voxelsize,[0 0 0]),[char(runno) '_QSM.nii.gz']);
end
