function img_out = coil_combine_kspace2(raw_k,fermir,fermiw)

[m,n,p,ncoil] = size(raw_k);

if ~exist('fermir','var')
    fermir = 8/512;
end
if ~exist('fermiw','var')
    fermiw = 2/512; %32? was 8
end
fermir=fermir*m;
fermiw=fermiw*m;
    
%%% Fermi Window
[x,y,z] = ndgrid(-m/2:m/2-1,-n/2:n/2-1,-p/2:p/2-1);
krad = sqrt(x.^2+y.^2+z.^2);
FW = 1./(1+exp((krad-fermir)/fermiw));

img_out = zeros(m,n,p,ncoil);
for icoil = 1:ncoil
    img = ifftnc(raw_k(:,:,:,icoil));
    ipha_low = angle(ifftnc(raw_k(:,:,:,icoil).*FW));
%     ipha_low = ifftnc(fftnc(angle(img)).*FW);
%     ipha_low = unwrap(angle(ifftnc(raw_k(:,:,:,icoil).*FW)));
    img_out(:,:,:,icoil) = abs(img).*exp(1i*(angle(img) - ipha_low));
end
% img_out = sum(img_out,4)./sum(abs(img_out),4);
img_out = sum(img_out,4);

end

