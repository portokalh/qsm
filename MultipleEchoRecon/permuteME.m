function [img] = permuteME(img,permute_order)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

figure(85)
subplot(211); z = size(img,3);
montageME(img,z,'1D'); title('Orientation 1');
img = permute(img,permute_order);
img = imrotate(img(end:-1:1,:,:,:),-90);
z = ceil(size(img,3)/2);
subplot(212); montageME(abs(img),z,'1D'); title('Orientation 2');

end

