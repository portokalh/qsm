% MEcorrect3.m has backup MEcorrect2.m

clear;
close all;

load raw;
load dims;
m = dims(1);
n = dims(2);
p = dims(3);
dim_divisor = m/256;

decraw = raw(1:dim_divisor:end,1:dim_divisor:end,1:dim_divisor:end,1:2);
img1 = abs(ifftnc(decraw(:,:,:,1)));
img2 = abs(ifftnc(decraw(:,:,:,2)));

ftx1 = fftshift(fft(fftshift(img1,1),[],1),1); % (1)  1D FT(x)
ftx2 = fftshift(fft(fftshift(img2,1),[],1),1); % (1)  1D FT(x)

Hi = angle(conj(ftx2).*ftx1);
Hi = unwrap(Hi);

[m,n,p,~] = size(decraw);
x = floor(-m/2):ceil(m/2)-1;
% xwin = floor(3*m/8)+1:floor(5*m/8);
xwin = floor(3*m/8)+1:floor(5*m/8);
yrow = floor(n/2);
zrow = floor(p/2);

w = 2;
vlength = (2*w+1)^2;
slopev = zeros(vlength,1);
yintv = slopev;
q = 0;
zvec = zrow-w:zrow+w;
yvec = yrow-w:yrow+w;

for z = zvec
    for y = yvec
        q = q+1;
        curve = fit((xwin)',Hi(xwin,y,z),'a*x+b');
        slopev(q) = curve.a;
        yintv(q) = curve.b;
    end
end

slope = median(slopev);
yint = yintv(floor(q/2)+1);

bvec = exp(1i*slope/dim_divisor*(1:dims(1))');

for q = 1:dims(4)
    for z = 1:dims(3)
        for y = 1:dims(2)
            raw(:,y,z,q) = raw(:,y,z,q).*bvec;  % (3)  correct phase    
        end
    end
end

imgb = ifftnc(raw(:,:,:,2)); % (4)  1D FT(x)
decimgb = abs(imgb(1:dim_divisor:end,1:dim_divisor:end,1:dim_divisor:end));

figure(1);
subplot(2,1,1); plot(x,abs(abs(img1(:,yrow,zrow))),x,abs(img2(:,yrow,zrow)),x,abs(decimgb(:,yrow,zrow))); 
legend('1st Echo','2nd Echo','2nd Echo Corrected'); title('Readout profile');
subplot(2,1,2); plot(x,Hi(:,yrow,zrow),x,x*slope+yint,'--k',x,Hi(:,yvec(1),zvec(1)),x,Hi(:,yvec(end),zvec(end)));
% subplot(2,1,2); plot(x,Hi(:,w+1,w+1),x,x*slope+yint,'--k');

figure(2);
show(abs(img1)-abs(img2), ...
      abs(img1)-abs(decimgb));


figure(4); dock;
show([abs(img1)+abs(img1) ...
      abs(img1)+abs(img2) ... 
      abs(img1)+decimgb]);

zslice = 64;
show3(cat(3,decimgb(:,:,zslice),img1(:,:,zslice),img2(:,:,zslice)));

