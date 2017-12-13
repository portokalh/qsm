function mskM = brainextM2(img,threshold_high)



load vox;
dd = ceil(size(img,1)/256); %256 is max dim for readout (incr speed)
% mag = mean(abs(img),4);


% Calculate skull-stripping mask
ref_nii = 'mskM.nii';
nechoes = size(img,4);

if exist([pwd '/' ref_nii],'file') == 2
    mskM = open_nii(ref_nii);
else
    mskM = zeros(size(img));

for k = 1:nechoes
    
    if nechoes > 1
        
        if k > 1
%             mskM = open_nii('mskM.nii');
%             mag = abs(img(:,:,:,k)).*single(mskM(:,:,:,k-1));
            mag = abs(img(:,:,:,k));
        else
            load TE
            mag = mean(abs(img(:,:,:,TE<.03)),4);
        end
    end
    niiname = 'mag.tst.nii';
    save_nii(make_nii(mag,vox,[0 0 0],16),niiname);


        fprintf('   Generating mask...\n');
        if nargin < 2
%             mskM(:,:,:,k) = stripMaskM(niiname,dd);
            mskM(:,:,:,k) = stripMaskM2(niiname,dd);
        elseif nargin < 3
%             mskM(:,:,:,k) = stripMaskM(niiname,dd,threshold_high);
            mskM(:,:,:,k) = stripMaskM2(niiname,dd,threshold_high(k));
        end
        
%         save_nii(make_nii(single(mskM),vox,[0 0 0],2),ref_nii);
% %         load stripMaskVars
%         num_morphM(k) = num_morph;
%         radiusM(k) = radius;
%         threshold_highM(k) = threshold_high(k);
%         delete stripMaskVars.mat
        fprintf(['   Showing echo #' num2str(k) '\n']);
        show3(single(mskM(:,:,:,k)).*abs(img(:,:,:,k)));
 
   
 
end 

% save(['stripMaskVarsM' num2str(k) '.mat'],'num_morphM','radiusM','threshold_highM');
nii = make_nii(single(mskM),vox,[0 0 0],2);
save_nii(nii,ref_nii);
end

if dd <= 2
    magM = single(mskM).*single(abs(img)); 
    nii = make_nii(magM,vox,[0 0 0],16);
    save_nii(nii,'magM.nii');
    figure(1); dock(1); show3(magM);
end


end

