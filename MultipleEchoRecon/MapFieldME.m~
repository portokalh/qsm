function [fld] = MapFieldME(phs,echo_spacing)
%% function [fld] = MapFieldME(phs,echo_spacing,dims)
%   MapFieldME.m: Summary of this function goes here
%   Detailed explanation goes here

if (ischar(phs) == 1)
    phs = readMEraw(phs);
end

gamma = 267.51336e6;
fld = zeros(size(phs));
fld(:,:,:,end) = [];

for k = 1:size(phs,4)-1
    fld(:,:,:,k) = (phs(:,:,:,k+1)-phs(:,:,:,k))/(gamma*echo_spacing);
end

fld = mean(fld,4);

figure(2);
subplot(235); show(fld,ceil(size(fld,3)/2)); title('Field Map');

% writeRaw(fld(:,:,ceil(size(fld,3)/2),:),'fieldmap.raw');



end

