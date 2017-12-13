function [imgODD,imgEVEN,imgALL] = combineME(img_type,method)
%% function [imgODD,imgEVEN,imgALL] = combineME(img_type,dims)
%   combineME.m: Combine multi-echo data to create enhanced images
%   Creates combined images of odd, even, and all echoes and displays the
%   results.  The combined image volumes are also saved with the extension
%   img_type.  This script uses readMEraw, which requires the input 
%   dimensions.  Dims can be a scalar or vector up to length 3.  If a
%   scalar, m, volume is assumed to be m x m x m.  If dims is 1x2, volume
%   is m x n x n.  If dims is 1x3, volume is m x n x p.
    
if ischar(img_type) == 1
    img = readMEraw(img_type);
    if ((strcmp(img_type,'frq') || strcmp(img_type,'frq.nii')) && ~exist('method','var'))
        method = 'avg';
    end
else
    img = img_type;
    img_type = inputname(1);
end

echoes = size(img,4);
z = ceil(size(img,3)/2);

if ~exist('method','var')
    method = 'fic';
end

switch method
    case 'avg'
        % Averaging
        imgODD = mean(img(:,:,:,1:2:echoes),4);
        imgEVEN = mean(img(:,:,:,2:2:echoes),4);
        imgALL = mean(img,4);
    case 'sum'
        % Summing
        imgODD = sum(img(:,:,:,1:2:echoes),4);
        imgEVEN = sum(img(:,:,:,2:2:echoes),4);
        imgALL = sum(img,4);
    case 'rms'
        % RMS - can be very slow
        imgODD = rms(img(:,:,:,1:2:echoes),4);
        imgEVEN = rms(img(:,:,:,2:2:echoes),4);
        imgALL = rms(img,4);
    case 'fic'
        % FIC image
        imgODD = ficME(img(:,:,:,1:2:echoes));
        imgEVEN = ficME(img(:,:,:,2:2:echoes));
        imgALL = ficME(img);
end

z = cslice(img);

subplot(241); show(img(:,:,:,1),z); title([img_type ' - 1st']);
subplot(242); show(imgODD,z); title([img_type ' - Odd']);
subplot(243); show(imgEVEN,z); title([img_type ' - Even']);
subplot(244); show(imgALL,z); title([img_type ' - All']);

subplot(2,4,5:8);
show([img(:,:,z,1) imgODD(:,:,z) imgEVEN(:,:,z) imgALL(:,:,z)]); colorbar;
title([img_type ' - 1st Echo, Odd Echoes, Even Echoes, All Echoes']); 

ref_nii = dir('mag.odd.nii');
if ~isempty(ref_nii)
    ref_nii = load_nii(ref_nii(1).name);
    
    niiODD = ref_nii; niiEVEN = ref_nii; niiALL = ref_nii;
    niiODD.img = imgODD;
    niiEVEN.img = imgEVEN;
    niiALL.img = imgALL;

    save_nii(niiODD,[img_type '.odd.nii']);
    save_nii(niiEVEN,[img_type '.evn.nii']);
    save_nii(niiALL,[img_type '.all.nii']);
elseif exist('vox.mat','file')
    load vox;
    niiODD = make_nii(imgODD,vox,[0 0 0],16);
    niiEVEN = make_nii(imgEVEN,vox,[0 0 0],16);
    niiALL = make_nii(imgALL,vox,[0 0 0],16);
    save_nii(niiODD,[img_type '.odd.nii']);
    save_nii(niiEVEN,[img_type '.evn.nii']);
    save_nii(niiALL,[img_type '.all.nii']);
else
    filename = ['odd.' img_type];
    writeRaw(imgODD,filename);

    filename = ['evn.' img_type];
    writeRaw(imgEVEN,filename);

    filename = ['all.' img_type];
    writeRaw(imgALL,filename);
end

end

