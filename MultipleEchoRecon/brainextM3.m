function [mskM,vmask] = brainextM3(img,threshold_high)



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
    

        if k == 1
            load TE
            niiname = 'mag.tst.nii';
            save_nii(make_nii(mean(abs(img(:,:,:,TE<.03)),4),vox,[0 0 0],16),niiname);
            mskM(:,:,:,k) = stripMaskM(niiname,dd,threshold_high(k)*2^15);
            msk1 = mskM(:,:,:,k);
            mat2nii(msk1);     
%             vmag = mean(abs(img(:,:,:,3:6)),4).*single(msk1);
%             vmask = (0 < vmag) & (vmag < .15);
           
          
        end      
        msk1 = open_nii('msk1.nii');    
        mag = abs(img(:,:,:,k)).*single(msk1);
        if k > 1
            [tri,hys] = hysteresis3d(mag,threshold_high(k),threshold_high(1),6);
    %         mskM(:,:,:,k) = ((mag > threshold_high(k)/2^15) & msk1);
            mskM(:,:,:,k) = hys&msk1; 
        else
            mskM(:,:,:,k) = msk1;
        end
        radius = 1;
        [xx,yy,zz] = ndgrid(-radius:radius);
        nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= radius;
        radius0 = 2; 
        nhood0=strel('ball', radius0, radius0, 0);
%         mskM(:,:,:,k) = imclose(mskM(:,:,:,k), nhood0);

        
        mskM(:,:,:,k) = imopen(mskM(:,:,:,k), nhood0);
        mskM(:,:,:,k) = imfill(mskM(:,:,:,k),'holes');
        % or imfill here?
        mskM(:,:,:,k) = imerode(logical(mskM(:,:,:,k)),nhood);
        mskM(:,:,:,k) = bwareaopen(mskM(:,:,:,k),1e5,6);
%         msktemp = mskM(:,:,:,k);
%         msktemp(vmask) = 0;
%         mskM(:,:,:,k) = msktemp;

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

