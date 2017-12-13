%%function [Freq] = MapFreqME(pfile)
% Before calling the function, put the pfiles of GRE sequences in a
% separate folder. Call the function in the created folder. Otherwise, you
% have to provide the full path for the pfile.
% enter MapFreqME() 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ORIG Input:
% function [Freq] = MapFreqME()
% Input:     The function will read in all pfiles with a name in 
%            the format of 'split_me_trains_gaps.*' from the current
%            directory.
% Output:
% Freq  -    Frequency maps in Hz, stored in binary float format 
%            nrows x ncols x nslices with the file name of echo#.freq
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Chunlei Liu, Duke University, 06/2009
% slg - 100715 V2 removed Dim as param MapFreq, use
%        dimensions from pfile header (as stored by civm 3D cartesian scans)
%        Changed raw data reader to use header dimensions and bl skip.
%        V3 added fft2c and fftnc and inverses, provided by Chunlei.
%        Wrapped GE header reader to notice missing Pfile.
% slg - 100716 V4 Now able to discern raw format and read short and int/EDR Pfiles.
% rmd - 110310 V_Multi_Echo - Only works for CIVM multi echo files reading
%        in split_me_trains_gaps.* files created by radish reconstruction
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%function [Freq] = MapFreqME()

% 
% options.dim = '3D';
% options.ncoil = 1;
% options.fermir = 8/512;
% options.fermiw = 8/512;

function img_out = CalcPhase2(img)

if (nargin < 2)
    ncoil = 1;
end

options.dim = '3D';
options.ncoil = 1;
options.fermir = 8/512;
options.fermiw = 8/512;

% if ( nargin == 2 && strcmp(options.dim,'3D'))
%     fprintf('3D...\n');
    opxres = size(img,1);
    opyres = size(img,2);
    opzres = size(img,3);
    ncoil = size(img,4);
    
    %%% Fermi Window
    [y,x,z] = meshgrid(-opyres/2:opyres/2-1,-opxres/2:opxres/2-1,-opzres/2:opzres/2-1);
    if (~isfield(options,'fermir'))
        fermir=8*opxres/512;%100,12,8,32; 4 for 3T
    else
        fermir=options.fermir*opxres;
    end
    
    if (~isfield(options,'fermiw'))
        fermiw=8*opxres/512;%12,8, 32;
    else
        fermiw=options.fermiw*opxres;
    end
    krad = sqrt(x.^2+y.^2+z.^2);
    FW = 1./(1+exp((krad-fermir)/fermiw));
    
    img_out = zeros(opxres,opyres,opzres,ncoil);
    %%% compute phase
    for icoil = 1:ncoil
        ipha_low = ifftnc(fftnc(img(:,:,:,icoil)).*FW);
        img_out(:,:,:,icoil) = abs(img(:,:,:,icoil)).*exp(1i*(angle(img(:,:,:,icoil)) - angle(ipha_low)));
    end
    img_out = angle(sum(img_out,4));
% end



end

function im = fftnc(d)
% im = fftnc(d)
%
% fftnc perfomrs a centered fftn
%
im = fftshift(fftn(fftshift(d)));
end

function im = ifftnc(d)
% im = ifftnc(d)
%
% ifftnc perfomrs a centered ifftn
%
im = fftshift(ifftn(fftshift(d)));
end

