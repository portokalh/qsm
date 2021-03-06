% MEcorrect.m

close all;

load dims
m = dims(1);
x = floor(-m/2):ceil(m/2)-1;

yrow = 64;
zrow = 64;

iftkx = fftshift(ifft(fftshift(raw(:,:,:,1:2),1),[],1),1); % (1)  1D IFT(k_x)

img1 = abs(ifftnc(raw(:,:,:,1)));
img2 = abs(ifftnc(raw(:,:,:,2)));

k1 = raw(:,yrow,zrow,1);
k2 = raw(:,yrow,zrow,2);
i1 = iftkx(:,yrow,zrow,1);
i2 = iftkx(:,yrow,zrow,2);

% Hi = angle(conj(i2).*i1);
Hi = angle(conj(iftkx(:,:,:,2)).*iftkx(:,:,:,1));  % (2)  Get 2a
Hi = unwrap(Hi);

[curve,goodness] = fit(x(65:192)',Hi((65:192),yrow,zrow),'a*x+b');
bvec = exp(1i*curve.a*x');

iftkxc = zeros(dims(1:3));
% rawb = iftkxc;

for z = 1:dims(3)
    for y = 1:dims(2)
    iftkxc(:,y,z) = iftkx(:,y,z,2).*bvec;  % (3)  correct phase    
%     imgb(mm,:,:) = abs(fftshift(ifft2(fftshift(squeeze(iftkxc(mm,:,:))))));
    end
end



rawb = fftshift(fft(fftshift(iftkxc,1),[],1),1); % (4)  1D FT(x)
k2b = rawb(:,yrow,zrow);


imgb = abs(ifftnc(rawb));                       % (5)  3D IFT(k_x)

figure(1);
subplot(3,1,1); plot(x,abs(img1(:,yrow,zrow)),x,abs(img2(:,yrow,zrow)),x,abs(imgb(:,yrow,zrow))); 
legend('1st Echo','2nd Echo','2nd Echo Corrected'); title('Readout profile');
% subplot(3,1,2); plot(x,Hi(:,yrow,zrow),'.',x,x*slope+yint);
subplot(3,1,2); plot(x,Hi(:,yrow-4,zrow-4),x,Hi(:,yrow+4,zrow+4),x,Hi(:,yrow,zrow),x,x*curve.a+curve.b,'--k');
subplot(3,1,3); plot(x,abs(k1),x,abs(k2),x,abs(k2b)); legend('k1','k2','k2b');

% figure(3);
% show(img1,img2,imgb);
% figure(4);
% show(img1-img2,img1-imgb,img2-imgb);
% 
% zslice = 64;
% show3(cat(3,imgb(:,:,zslice),img1(:,:,zslice),img2(:,:,zslice)));
