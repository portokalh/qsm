function img_out = coil_combine_complex(img,mask,padsize)
% img is a 4D image volume where the 4th dim is coils

padsize = [ 0 0 0];
load vox
[m,n,p,ncoil] = size(img);

phi = zeros(m,n,p,ncoil,'single');

% need to unwrap phases so that the combination image doesn't combine
% combine conflicting phase wraps from multiple coils
for icoil = 1:ncoil
% phi(:,:,:,icoil) = iHARPERELLA(angle(img(:,:,:,icoil)),mask,'voxelsize',vox);
        phi(:,:,:,icoil) = unwrap_phase(angle(img(:,:,:,icoil)),vox,padsize);
    if nargin > 1
        phi(:,:,:,icoil) = V_SHARP(phi(:,:,:,icoil),single(mask),'smvsize',16);%larger sphere make the results more accurate and take longer
    end
end

% % Weighted sum of the phase images
if nargin < 2
    magweights = abs(img)./repmat(sum(abs(img),4),[1 1 1 ncoil]);%fractional weights of the total magnitude
    img_out = sum(abs(img),4).*exp(1i*sum(magweights.*phi,4)); 
else
    % Weighted sum of complex images, may produce incorrect magnitude images
%     img_out = mean(abs(img).*exp(1j*phi),4); % why should this be different from above?
    img_out = sum(abs(img),4).*exp(1j*angle(mean(abs(img).*exp(1j*phi),4))); % why should this be different from above?
end

end

