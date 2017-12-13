function raw = phaseME2(raw)
% phaseME2.m is a 2D implementation of 
% Only phase corrects in readout (x) direction

for d = 1:3
    


dims = size(raw);
csize = 128;
dim_divisor = dims(1)/csize;
M = (dims(1)-csize)/2+1:(dims(1)+csize)/2;
N = (dims(2)-csize*dims(2)/dims(1))/2+1:(dims(2)+csize*dims(2)/dims(1))/2;
P = (dims(3)-csize*dims(3)/dims(1))/2+1:(dims(3)+csize*dims(3)/dims(1))/2;

fprintf('      Performing transforms...\n');
decraw = raw(M,N,P,:);
iftkx1 = fftshift(ifft(fftshift(decraw(:,:,:,1),1),[],1),1); % (1)  1D IFT(k_x)
iftkx2 = fftshift(ifft(fftshift(decraw(:,:,:,2),1),[],1),1); % (1)  1D IFT(k_x)
raw(:,:,:,2:2:end) = fftshift(ifft(fftshift(raw(:,:,:,2:2:end),1),[],1),1);

img1 = angle(ifftnc(decraw(:,:,:,1)));
img2 = angle(ifftnc(decraw(:,:,:,2)));

% i1 = iftkx(:,yrow,zrow,1);
% i2 = iftkx(:,yrow,zrow,2);

fprintf('      Calculating shift vector...\n');
[m,n,p,~] = size(decraw);
x = floor(-m/2):ceil(m/2)-1;
xwin = floor(3*m/8)+1:floor(5*m/8);
% xwin = floor(m/4)+1:floor(3*m/4);
yrow = floor(n/2);
zrow = floor(p/2);
w = 5;
s = 1;
zvec = zrow-w*s:s:zrow+w*s;
yvec = yrow-w*s:s:yrow+w*s;
Hi = angle(conj(iftkx2(:,yvec,zvec)).*iftkx1(:,yvec,zvec));
Hi = unwrap(Hi);

v = (2*w+1);
h = ceil(v/2);
slopev = zeros(v^2,1);
yintv = slopev;
q = 0;

for z = 1:v
    for y = 1:v
        q = q+1;
        curve = fit((xwin-ceil(m/2))',Hi(xwin,y,z),'a*x+b');
        slopev(q) = curve.a;
        yintv(q) = curve.b;
    end
end

slopev = removeoutliers(slopev);
slope = median(slopev);
slope2 = mean(slopev);
sig = std(slopev);

% Get ride of data that didn't unwrap properly
while (2*abs(slope-slope2)/abs(slope+slope2) > .025)
    figure(1); dock; subplot(3,1,3); hist(slopev,w^2);
    title('Phase slope estimate distribution');
    fprintf('         median=%.3f, mean=%.3f, sigma=%.4f\n',slope,slope2,sig);
    slopev = slopev(slopev < slope2+2*sig);
    slopev = slopev(slopev < slope2+2*sig);
    slope = median(slopev);
    slope2 = mean(slopev);
    sig = std(slopev); 
end

yint = yintv(floor(q/2)+1);

bvec = exp(1i*slope/dim_divisor*(1:dims(1))');
% bvec = exp(1i*slope*(1:dims(1))');

fprintf('         median=%.3f, mean=%.3f, sigma=%.4f\n',slope,slope2,sig);
fprintf('      Correcting kspace data...\n');
for q = 2:2:dims(4)
    
    for z = 1:dims(3)
        for y = 1:dims(2)
            raw(:,y,z,q) = raw(:,y,z,q).*bvec;  % (3)  correct phase    
        end
    end
    
end

raw(:,:,:,2:2:end) = fftshift(fft(fftshift(raw(:,:,:,2:2:end),1),[],1),1); % (4)  1D FT(x)

% for z = 1:dims(3)
%     for y = 1:dims(2)
%     iftkxc(:,y,z) = iftkx(:,y,z,2).*bvec;  % (3)  correct phase    
% %     imgb(mm,:,:) = abs(fftshift(ifft2(fftshift(squeeze(iftkxc(mm,:,:))))));
%     end
% end

rawb = raw(M,N,P,2); 
k1 = decraw(:,yrow,zrow,1);
k2 = decraw(:,yrow,zrow,2);
k2b = rawb(:,yrow,zrow);


imgb = angle(ifftnc(rawb));                       % (5)  3D IFT(k_x)

figure(1);
subplot(3,1,1); plot(x,abs(abs(k1)),x,abs(k2),x,abs(k2b)); 
legend('1st Echo','2nd Echo','2nd Echo Corrected'); title('Readout profile');
title(['Readout profile - shift correction = ' ... 
        num2str(-slope*dims(1)/2/pi/dim_divisor,2) ' pixels']);
% subplot(3,1,2); plot(x,Hi(:,yrow,zrow),'.',x,x*slope+yint);
subplot(3,1,2); plot(x,Hi(:,h,h),x,x*slope+yint,'--k',x,Hi(:,1,1),x,Hi(:,v,v));
title(['Image phase, readout slope = ' num2str(slope)]);
% subplot(3,1,1); plot(x,abs(k1),x,abs(k2),x,abs(k2b)); legend('k1','k2','k2b');
% subplot(3,1,3); hist(slopev,w^2);
% title('Phase slope estimate distribution');

figure(3);
show(img1,img2,imgb);
figure(4);
show(img1-img2,img1-imgb,img2-imgb);
% 
% zslice = 64;
% show3(cat(3,imgb(:,:,zslice),img1(:,:,zslice),img2(:,:,zslice)));

raw = permute(raw,[2 1 3 4]);

end

end
