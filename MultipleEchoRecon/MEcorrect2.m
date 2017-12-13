% MEcorrect2.m

close all;

load dims
m = dims(1);
n = dims(2);
p = dims(3);
x = floor(-m/2):ceil(m/2)-1;

xwin = floor(3*m/8)+1:floor(5*m/8);
yrow = floor(n/2);
zrow = floor(p/2);
% 
load raw;
img(:,:,:,1) = ifftnc(raw(:,:,:,1));
img(:,:,:,2) = ifftnc(raw(:,:,:,2));

img = abs(img); 
% For some reason, this shift correction requires me to 
% take the absolute value of the image.  This is probably
% bad math.

% ftx = zeros(size(raw));
% tic
% for k = 1:dims(4)
%     for mm = 1:m
%         ftx(mm,:,:,k) = ifft2c(raw(mm,:,:,k));
%     end
% end
% toc
                
                
ftx = fftshift(fft(fftshift(img,1),[],1),1); % (1)  1D FT(x)
n = 0;
w = 2;
p = 2*w+1;


tic
Hi = angle(conj(ftx(:,yrow-w:yrow+w,zrow-w:zrow+w,2)) ...
                .*ftx(:,yrow-w:yrow+w,zrow-w:zrow+w,1));  % (2)  Get 2a
Hi = unwrap(Hi);
toc

slopev = zeros(p^2,1);
yintv = slopev;


for z = 1:p
    for y = 1:p
        n = n+1;
        curve = fit((xwin-128)',Hi(xwin,y,z),'a*x+b');
        slopev(n) = curve.a;
        yintv(n) = curve.b;
    end
end

slope = median(slopev);
yint = yintv(floor(n/2)+1);

bvec = exp(1i*slope*x');

ftxc = zeros(dims(1:3));

for z = 1:dims(3)
    for y = 1:dims(2)
        ftxc(:,y,z) = ftx(:,y,z,2).*bvec;  % (3)  correct phase    
    end
end

imgb = abs(fftshift(ifft(fftshift(ftxc,1),[],1),1)); % (4)  1D FT(x)
k2b = ftxc(:,yrow,zrow);

figure(1);
subplot(2,1,1); plot(x,abs(abs(img(:,yrow,zrow,1))),x,abs(img(:,yrow,zrow,2)),x,abs(imgb(:,yrow,zrow))); 
legend('1st Echo','2nd Echo','2nd Echo Corrected'); title('Readout profile');
subplot(2,1,2); plot(x,Hi(:,w+1,w+1),x,x*slope+yint,'--k',x,Hi(:,1,1),x,Hi(:,p,p));
% subplot(2,1,2); plot(x,Hi(:,w+1,w+1),x,x*slope+yint,'--k');

figure(2);
show(abs(img(:,:,:,1))-abs(img(:,:,:,2)), ...
      abs(img(:,:,:,1))-abs(imgb));


figure(4); dock;
show([abs(img(:,:,:,1))+abs(img(:,:,:,1)) ...
      abs(img(:,:,:,1))+abs(img(:,:,:,2)) ... 
      abs(img(:,:,:,1))+imgb]);

% zslice = 64;
% show3(cat(3,imgb(:,:,zslice),img(:,:,zslice,1),img(:,:,zslice,2)));

