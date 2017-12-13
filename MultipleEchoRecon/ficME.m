function [ficimage] = ficME(img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if size(img,4) > 1
    ficimage=fft(img,[],4);  
    ficimage=abs(ficimage);  
    ficimage=max(ficimage,[],4);
else
    ficimage = img;
end

end

