function [r,b] = nanLinearFit4D(Volume4D,TE,mask,SNR_for_weighting)
%%  T2mapMaskedWeighted.m: Creates a T2 (T2*) map using a 4D data set.  
%   Russell Dibb, Duke University, 03/27/2013
%
%   T2 = T2mapfastW(Volume4D,mask,TE,mask,SNR_for_weighting)
%   T2 - T2 map for the image volue
%   Volume4D - Time series of image volumes
%   TE - vector of TE values (time points)
%   mask - mask determines index of voxels to calculate T2 (opt.)
%   SNR_for_weighting - vector of relative SNR of the image volumes (opt.)

[m,n,p,echoes] = size(Volume4D);
Volume4D = reshape(Volume4D,[m*n*p echoes]);

if ~exist('SNR_for_weighting','var')
    SNR_for_weighting = ones(echoes,1);
end

if ~exist('mask','var')
%     i = ones(m*n*p,1); % calculates every voxel
    i = all(~isinf(Volume4D),2); % should work faster for premasked data
else 
    i = mask(:);
end

r = NaN(size(i));
b = r;

data = Volume4D(i==1,:)';
idata = ~isnan(data);
ndata = sum(idata);
i2 = ndata > 1;
data2 = data(:,i2);
idata2 = ~isnan(data2);
ndata2 = sum(idata2);
v1 = zeros(2,size(data2,2));


parfor k = 1:size(data2,2);
    v1(:,k) = lscov([ones(ndata2(k),1) TE(idata2(:,k))'],data2(idata2(:,k),k),SNR_for_weighting(idata2(:,k)));
end

r2 = zeros(size(data,2),1);
r2(i2) = v1(2,:);
dd = data(:,~i2);
dd(isnan(dd)) = 0;
r2(~i2) = nanmean(TE*dd);
r(i==1) = r2;

% r(i==1) = v1(2,:);
% b(i==1) = v1(1,:);
r = reshape(r,[m n p]);
% b = reshape(b,[m n p]);
% m(m < 0) = 0;
% T2 = 1./m; % you will get NaNs

end

