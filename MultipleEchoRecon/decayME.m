% decayME.m
clc;
clear;

data = load('raw.mat');
load TE;
load vox;
load fov;
load dims;

necho = 10;
T2 = 14e-3;
dwell = 4e-6;
T2v = exp(-TE/T2);
noise_mean = 1;
noise_std = .05;

raw = zeros([size(data.raw) necho]);
mag = raw;
raw(:,:,:,1) = data.raw;
img(:,:,:,1) = ifftnc(data.raw);
% data.raw = data.raw/T2v(1);
% noiseROI1 = img(11:30,11:30,61:70,1);
rraw = real(data.raw/T2v(1));
iraw = imag(data.raw/T2v(1));


for k = 2:necho
    k
    T2v = repmat(TE-dims(1)/2*dwell:dwell:TE+(dims(1)-1)*dwell,[1 dims(2) dims(3)]);
    T2v = exp(-TE(k)./T2v);
%     rnoise = random('normal',noise_mean,noise_std(k-1),size(data.raw));
%     inoise = random('normal',noise_mean,noise_std(k-1),size(data.raw))*1i;
    rnoise = random('normal',noise_mean,noise_std,size(data.raw));
    inoise = random('normal',noise_mean,noise_std,size(data.raw));
    raw(:,:,:,k) = T2v(k)*complex(rraw.*rnoise,iraw.*inoise);
    img(:,:,:,k) = ifftnc(raw(:,:,:,k));
%     noiseROIk = img(11:30,11:30,61:70,k);
%     scalar = mean(noiseROI1(:))/mean(noiseROIk(:));
%     img(:,:,:,k) = img(:,:,:,k)*scalar;

end

figure(5);
subplot(3,1,1);
montageME(log(1+abs(raw)),64,'1D'); colorbar;
subplot(3,1,2);
montageME(abs(img),64,'1D'); colorbar;
subplot(3,1,3);
montageME(angle(img),64,'1D'); colorbar;
