function msk = brainext(img,threshold_zero,threshold_mult)



load vox;
dd = ceil(size(img,1)/256); %256 is max dim for readout (incr speed)
% mag = mean(abs(img),4);
if size(img,4) > 1
    mag = mefic(abs(img));
else 
    mag = abs(img);
end
niiname = 'mag.tst.nii';
save_nii(make_nii(mag,vox,[0 0 0],16),niiname);

% Calculate skull-stripping mask
ref_nii = 'msk.nii';

if exist([pwd '/' ref_nii],'file') == 2
    msk = load_nii(ref_nii);
    msk = msk.img;
else
    fprintf('   Generating mask...\n');
    if nargin < 2
        msk = stripMask(niiname,dd);
    elseif nargin < 3
        msk = stripMask(niiname,dd,threshold_zero);
    else
        msk = stripMask(niiname,dd,threshold_zero,threshold_mult);
    end
    nii = make_nii(single(msk),vox,[0 0 0],2);
    save_nii(nii,ref_nii);
end

if dd <= 2
    mag = single(msk).*single(mag); 
    nii = make_nii(mag,vox,[0 0 0],16);
    save_nii(nii,'mag.nii');
    figure(1); dock(1); show(mag);
end
    
end

