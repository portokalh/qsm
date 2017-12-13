function [img] = niis2stack(nii_name,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

files = dir([nii_name '*']); 

echoes = length(files);

if nargin > 1
    option = varargin{1};
    if (strcmp(option,'X') || strcmp(option,'Freq'))
        load TE;
        B0 = input('What is B0?: ');
        coeffs = Xcoeff(B0,TE);
    end
else
    option = '';
end

if exist('dims.mat','file') == 2
    load dims;
    img = zeros([dims(1:3) echoes]);
end 

back = 0; time = 0;
for k = 1:echoes
    [back,time] = progress(k,echoes,'Stacking echo',back,time);
    filename = files(k).name;
    img(:,:,:,k) = open_nii(filename);
    if strcmp(option,'X')
        img(:,:,:,k) = img(:,:,:,k)*coeffs.X(k);
    end
    if strcmp(option,'Freq')
        img(:,:,:,k) = img(:,:,:,k)*coeffs.Freq(k);
    end
    clear(sprintf('%s',filename(1:end-4)));
end

% montageME(img);


end

