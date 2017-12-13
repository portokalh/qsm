function bvec = shiftME(raw)
% shiftME.m has backup MEcorrect3.m

% close all;

dims = size(raw);
csize = 128;
dim_divisor = dims(1)/csize;
M = (dims(1)-csize)/2+1:(dims(1)+csize)/2;
N = floor((dims(2)-csize*dims(2)/dims(1))/2)+1:floor((dims(2)+csize*dims(2)/dims(1))/2);
P = floor((dims(3)-csize*dims(3)/dims(1))/2)+1:floor((dims(3)+csize*dims(3)/dims(1))/2);

fprintf('      Performing transforms...\n');
% decraw = raw(M,N,P,1:2);
% img1 = abs(ifftnc(decraw(:,:,:,1)));
% img2 = abs(ifftnc(decraw(:,:,:,2)));

decraw = raw(M,N,P,:);
img1 = max(abs(ifftnc(decraw(:,:,:,1:2:end))),[],4);
img2 = max(abs(ifftnc(decraw(:,:,:,2:2:end))),[],4);

ftx1 = fftshift(fft(fftshift(img1,1),[],1),1); % (1)  1D FT(x)
ftx2 = fftshift(fft(fftshift(img2,1),[],1),1); % (1)  1D FT(x)

img1 = abs(ifftnc(decraw(:,:,:,1)));
img2 = abs(ifftnc(decraw(:,:,:,2)));

fprintf('      Calculating shift vector...\n');
[m,n,p,~] = size(decraw);
x = floor(-m/2):ceil(m/2)-1;
xwin = floor(3*m/8)+1:floor(5*m/8);
% xwin = floor(m/4)+1:floor(3*m/4);
yrow = floor(n/2);
zrow = floor(p/2);
w = 5;
s = 2;
zvec = zrow-w*s:s:zrow+w*s;
yvec = yrow-w*s:s:yrow+w*s;
Hi = angle(conj(ftx2(:,yvec,zvec)).*ftx1(:,yvec,zvec));
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
while ((2*abs(slope-slope2)/abs(slope+slope2) > .025) && abs(slope) > .025) %abs(slope)?
    figure(1); dock; subplot(3,1,3); hist(slopev,w^2);
    title('Phase slope estimate distribution');
    fprintf('         median=%.3f, mean=%.3f, sigma=%.4f\n',slope,slope2,sig);
    slopev = slopev(slopev < slope2+2*sig);
    slopev = slopev(slopev < slope2+2*sig);
    slope = median(slopev);
    slope2 = mean(slopev);
    sig = std(slopev); 
end

yint = yintv([1 floor(q/2)+1 q]);

decbvec = exp(1i*slope/dim_divisor*(1:dims(1))');
bvec = exp(1i*slope*(1:dims(1))');


fprintf('         mean=%.3f, median=%.3f, sigma=%.4f\n',slope,slope2,sig);
fprintf('      Correcting kspace data...\n');

% if slope > .025;
%     for q = 2:2:dims(4)
%         for z = 1:dims(3)
%             for y = 1:dims(2)
%                 raw(:,y,z,q) = raw(:,y,z,q).*bvec;  % (3)  correct phase    
%             end
%         end
%     end  
% end

fprintf('      Generating image results...\n');

decrawb = zeros(m,n,p);
for z = 1:p
    for y = 1:n
        decrawb(:,y,z) = decraw(:,y,z).*decbvec(M);  % (3)  correct phase    
    end
end

decimgb = abs(ifftnc(decrawb)); % (4)  1D FT(x)

figure(1); dock;
subplot(3,1,1); plot(x,abs(abs(img1(:,yrow,zrow))),x,abs(img2(:,yrow,zrow)),x,abs(decimgb(:,yrow,zrow))); 
legend('1st Echo','2nd Echo','2nd Echo Corrected'); title('Readout profile');
title(['Readout profile - shift correction = ' num2str(-slope*dims(1)/2/pi,2) ' pixels']);
subplot(3,1,2); plot(x,Hi(:,h,h),x,x*slope+yint(1),'--k', ...
                     x,x*slope+yint(2),'--k',x,x*slope+yint(3),'--k', ...
                     x,Hi(:,1,1),x,Hi(:,end,end));
title(['Image phase, readout slope = ' num2str(slope)]);
% subplot(2,1,2); plot(x,Hi(:,w+1,w+1),x,x*slope+yint,'--k');
subplot(3,1,3); hist(removeoutliers(slopev),w^2);
title('Phase slope estimate distribution');

figure(2); dock; 
show([abs(img1)-abs(img2) abs(img1)-abs(decimgb)]); 
title('Difference between echoes 1 and 2 before and after correction');

% figure(4);
% show([abs(img1)+abs(img1) ...
%       abs(img1)+abs(img2) ... 
%       abs(img1)+decimgb]);

% zslice = ceil(p/2);
% show3(cat(3,decimgb(:,:,zslice),img1(:,:,zslice),img2(:,:,zslice)));

fprintf('      Even echo shift correction complete!\n');

end

