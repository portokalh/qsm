function img = orientME(img,orientation)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fprintf('   Reorienting image to standard orientation...\n');
if orientation(1) == -1             % frequency dir flip
    img = img(end:-1:1,:,:,:);
end
if orientation(2) == -1             % phase 1 dir flip
    img = img(:,end:-1:1,:,:);
end
if orientation(3) == -1             % phase 2 dir flip
    img = img(:,:,end:-1:1,:);
end


end

