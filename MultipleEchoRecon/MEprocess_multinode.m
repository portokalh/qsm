function MEprocess_multinode(workpath)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

cd(workpath);
raw = abs(open_nii([workpath '/img.nii']));
mask = single(open_nii([workpath '/mask.nii'])); 

% SNR and T2 calculations
% [SNR,NoiseParams,SignalParams,mask_new,meanSNR] = autoSNR(raw,mask);
% save([workpath '/meanSNR.mat'],'SNR','NoiseParams','SignalParams','mask_new','meanSNR');
autoSNR(raw,mask(:,:,:,1)); % since I've cd'd to the workpath, should save in the proper location, right?

load([workpath '/meanSNR.mat']);
load([workpath '/dataprops.mat']);
T2idx = 1:dims(4);
T2idx = T2idx(meanSNR > 3);%5
if length(T2idx) == 1
    T2idx = 1:2;
end
T2 = T2mapMaskedWeighted(raw(:,:,:,T2idx),TE(T2idx),mask(:,:,:,1),meanSNR(T2idx));
save_nii(make_nii(T2,vox,[0 0 0],16),[workpath '/T2.nii']);

end

