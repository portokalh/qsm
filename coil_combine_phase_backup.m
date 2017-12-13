function img_out = coil_combine_phase(img)
% img is a 4D image volume where the 4th dim is coils

padsize = [0 0 6];
load vox
[m,n,p,ncoil] = size(img);

phiuwp = zeros(m,n,p,ncoil);

%%% Chunlei - weight unwrapped phase, no filter until after weight/sum?
for icoil = 1:ncoil
    phiuwp(:,:,:,icoil) = unwrap_phase(angle(img(:,:,:,icoil)),vox,padsize);
%     philow(:,:,:,icoil) = fermi(phiuwp(:,:,:,icoil)-mean(phiuwp,4), ...
%                                 filter_ratio,vox,padsize);            
end

% img_out=sum(magweights.*exp(1i*(phiuwp)),4);
% img_out=sum(abs(img).*exp(1i*(phiuwp)),4)./sum(abs(img),4); %weighting, division by total signal from both channels
% img_out=sum(abs(img).*exp(1i*(phiuwp)),4); %This is what works

magweights = abs(img)./repmat(sum(abs(img),4),[1 1 1 ncoil]);
img_out = abs(sum(img,4)).*exp(1i*sum(magweights.*phiuwp,4)); 
%%%



%%% best method so far is direct summing?  Why?

% %%% Wei's implementation - bad
% for icoil = 1:ncoil
%     phiuwp(:,:,:,icoil) = unwrap_phase(angle(img(:,:,:,icoil)),vox,padsize);
%     philow(:,:,:,icoil) = fermi(phiuwp(:,:,:,icoil)-mean(phiuwp,4), ...
%                                 filter_ratio,vox,padsize);            
% end
% img_out=sum(abs(img).*exp(1i*(phiuwp-philow)),4);
% %%%

%%% Wei's implementation, do not subtract mean phase  -- really bad
% for icoil = 1:ncoil
%     phiuwp(:,:,:,icoil) = unwrap_phase(angle(img(:,:,:,icoil)),vox,padsize);
%     philow(:,:,:,icoil) = fermi(phiuwp(:,:,:,icoil), ...
%                                 filter_ratio,vox,padsize);            
% end
% img_out=sum(abs(img).*exp(1i*(phiuwp-philow)),4);
%%%




end

