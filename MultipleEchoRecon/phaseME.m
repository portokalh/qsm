function raw = phaseME(raw)
% phaseME.m has backup MEcorrect.m
% Only phase corrects in readout (x) direction

dims = size(raw);
if dims(1) >= 256
    csize = 128;
elseif ((mod(dims(1),64) == 0) && (mod(dims(2),64) ==0))
    csize = 64;
end
if csize > dims(1)
    csize = dims(1);
end
dim_divisor = dims(1)/csize;
if ((dim_divisor == 3) || (dim_divisor == 6))
    csize = dim_divisor/2*csize;
    dim_divisor = dims(1)/csize;
end
M = (dims(1)-csize)/2+1:(dims(1)+csize)/2;
N = (dims(2)-csize*dims(2)/dims(1))/2+1:(dims(2)+csize*dims(2)/dims(1))/2;
if dims(3) >= csize
    P = (dims(3)-csize*dims(3)/dims(1))/2+1:(dims(3)+csize*dims(3)/dims(1))/2;
else
    P = 1;
end

fprintf('      Performing transforms...\n');
decraw = raw(M,N,P,:);
iftkx1 = fftshift(ifft(fftshift(decraw(:,:,:,1),1),[],1),1); % (1)  1D IFT(k_x)
iftkx2 = fftshift(ifft(fftshift(decraw(:,:,:,2),1),[],1),1); % (1)  1D IFT(k_x)
raw(:,:,:,2:2:end) = fftshift(ifft(fftshift(raw(:,:,:,2:2:end),1),[],1),1);

img1 = angle(ifftnc(decraw(:,:,:,1)));
img2 = angle(ifftnc(decraw(:,:,:,2)));

fprintf('      Calculating shift vector...\n');
[m,n,p,~] = size(decraw);
x = floor(-m/2):ceil(m/2)-1;
xwin = floor(3*m/8)+1:floor(5*m/8);
yrow = floor(n/2);
zrow = ceil(p/2);
w = 5;
s = 1;
if dims(3) >= csize
    zvec = zrow-w*s:s:zrow+w*s;
else
    zvec = 1;
end
yvec = yrow-w*s:s:yrow+w*s;
Hi = angle(conj(iftkx2(:,yvec,zvec)).*iftkx1(:,yvec,zvec));
Hi = unwrap(Hi);

v = (2*w+1);
h = ceil(v/2);
slopev = zeros(v*length(zvec),1);
yintv = slopev;
q = 0;

for z = 1:length(zvec)
    for y = 1:v
        q = q+1;
        curve = fit((xwin-ceil(m/2))',double(Hi(xwin,y,z)), ...
                                      'a*x+b','Startpoint',[0 0]);
        slopev(q) = curve.a;
        yintv(q) = curve.b;
    end
end

slopev = removeoutliers(slopev);
slope = median(slopev);
slope2 = mean(slopev);
sig = std(slopev);

% Get rid of data that didn't unwrap properly
nruns = 0;
while ((2*abs(slope-slope2)/abs(slope+slope2) > .01) && (nruns < 10))
    nruns = nruns+1;
    figure(81); dock; subplot(3,1,3); hist(slopev,w^2);
    title('Modified phase slope estimate distribution');
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
raw(:,:,:,2:2:end) = raw(:,:,:,2:2:end).*repmat(bvec,[1 dims(2:3) floor(dims(4)/2)]);

raw(:,:,:,2:2:end) = fftshift(fft(fftshift(raw(:,:,:,2:2:end),1),[],1),1); % (4)  1D FT(x)

rawb = raw(M,N,P,2); 
k1 = decraw(:,yrow,zrow,1);
k2 = decraw(:,yrow,zrow,2);
k2b = rawb(:,yrow,zrow);

imgb = angle(ifftnc(rawb));                       % (5)  3D IFT(k_x)

figure(81);
subplot(3,1,1); plot(x,abs(abs(k1)),x,abs(k2),x,abs(k2b)); 
legend('1st Echo','2nd Echo','2nd Echo Corrected'); title('Readout profile');
title(['Readout profile - shift correction = ' ... 
        num2str(-slope*dims(1)/2/pi/dim_divisor,2) ' pixels']);

if z == 1
    zv = z;
    zh = z;
else
    zv = v;
    zh = h;
end
    
subplot(3,1,2); plot(x,Hi(:,h,zh),x,x*slope+yint,'--k',x,Hi(:,1,1),x,Hi(:,v,zv), ...
    x,x*slope+yintv(1),'--k',x,x*slope+yintv(end),'--k');
title(['Image phase, readout slope = ' num2str(slope)]);
subplot(3,1,3); hist(slopev,w^2);
title('Phase slope estimate distribution');

figure(3);
show(img1,img2,imgb);
end
