function img_out = coil_combine_kspace(raw_k,filter_ratio)

if ~exist('filter_ratio','var')
    filter_ratio = .4;
end
padsize = [0 0 6];
load vox
[m,n,p,ncoil] = size(raw_k);

img_out = zeros(m,n,p,ncoil);
img = img_out;
phiuwp = img;
philow = img;

for icoil = 1:ncoil
    img(:,:,:,icoil) = ifftnc(raw_k(:,:,:,icoil));
    phiuwp(:,:,:,icoil) = unwrap_phase(angle(img(:,:,:,icoil)),vox,padsize);
    philow(:,:,:,icoil) = fermi(phiuwp(:,:,:,icoil)-mean(phiuwp,4), ...
                                filter_ratio,vox,padsize);            
end

img_out=sum(abs(img).*exp(1i*(phiuwp-philow)),4);


end



% 
%         disp(ke)
%         kspace=readPfileFun(Pfile,ke);
%         magecho=zeros(size(kspace));
%         phaseecho=magecho;
%         for kc=1:nc
%             [magecho(:,:,:,kc), phaseecho(:,:,:,kc)]=CalMagPhase_v1(kspace,kc);
%         end
%         disp('image reconstructed....')
%         phaseuwp=zeros(size(phaseecho));
%         for kc=1:nc
%             phaseuwp(:,:,:,kc)=MRPhaseUnwrap_v1(phaseecho(:,:,:,kc),SpatialRes,padsize);
%         end
%         clear phaseecho
%         phaseuwpmean=mean(phaseuwp,4);
%         disp('calc diff ...')
%         phi0_Low=phaseuwp;
%         for ii=1:nc
%             diff=phaseuwp(:,:,:,ii)-phaseuwpmean;
%             phi0_Low(:,:,:,ii) = Fermi_v1(diff,0.4,SpatialRes,padsize);
%         end
%         disp(' Fermi Filtering done....')
%         clear diff
%         phasecorr=phaseuwp-phi0_Low;
%         clear phaseuwp phi0_Low 
%         newcomplex=magecho.*exp(1i*phasecorr);
%         magmean=mean(magecho,4);
%         if ke==5
%             ezsave_nii_v1(magmean, 'magecho5.nii',[],SpatialRes);
%         end
%         newcomplex=sum(newcomplex,4);
%         newphase=angle(newcomplex);
%         disp(' newphase done....')
%         phaseuwp=MRPhaseUnwrap_v1(newphase,SpatialRes,padsize);
%         save(['magmean' num2str(ke)], 'magmean','-v7.3')
%         save(['phaseuwp' num2str(ke)], 'phaseuwp','-v7.3')
%         magall(:,:,:,ke)=magmean;
%         phaseall=phaseall+phaseuwp;
%         disp('ke')
%         clear magmean newcomplex magecho

