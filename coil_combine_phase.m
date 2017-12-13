function img_out = coil_combine_phase(img)
% img is a 4D image volume where the 4th dim is coils

padsize = [0 0 6];
load vox
[m,n,p,ncoil] = size(img);

phiuwp = zeros(m,n,p,ncoil);

% need to unwrap phases so that the combination image doesn't combine
% combine conflicting phase wraps from multiple coils
for icoil = 1:ncoil
    phiuwp(:,:,:,icoil) = unwrap_phase(angle(img(:,:,:,icoil)),vox,padsize);          
end

% Weighted sum of the phase images
magweights = abs(img)./repmat(sum(abs(img),4),[1 1 1 ncoil]);
img_out = abs(sum(img,4)).*exp(1i*sum(magweights.*phiuwp,4)); 


end

