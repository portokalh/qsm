function [img] = Prestore(img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if exist('shift_vector.mat','file') == 2
    load shift_vector
    for k = 1:size(img,4)
        img(:,:,:,k) = circshift(img(:,:,:,k),-shift_vector);
    end
end

if exist('permute_vector.mat','file') == 2
    load permute_vector

    unpermute(1) = find((permute_vector == 1));
    unpermute(2) = find((permute_vector == 2));
    unpermute(3) = find((permute_vector == 3));
    
    img = permute(img,[unpermute 4]);
end



end

