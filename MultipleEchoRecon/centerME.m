function img = centerME(img)
% centerME.m: Helps you center an image volume (up to 5D) with GUI.

dims = size(img);

mc = ceil(dims(1)/2);
nc = ceil(dims(2)/2);
zc = ceil(dims(3)/2);
cmask1 = zeros(dims(1),1);
cmask1(mc-floor(mc/4):floor(mc/8):mc+floor(mc/4)) = 1;
cmask2 = zeros(1,dims(2));
cmask2(nc-floor(nc/4):floor(nc/8):nc+floor(nc/4)) = 1;
cmask3 = zeros(1,1,dims(3));
cmask3(zc-floor(zc/4):floor(zc/8):zc+floor(zc/4)) = 1;
cmask = repmat(cmask1,[1 dims(2:3)])+repmat(cmask2,[dims(1) 1 dims(3)])+ ...
               repmat(cmask3,[dims(1:2) 1]);

figure(90); dock;
subplot(311); show(squeeze(abs(img(mc,:,:,1,1)))'); 
subplot(312); show(squeeze(abs(img(:,nc,:,1,1)))');
subplot(313); show(squeeze(abs(img(:,:,zc,1,1)))); 
cstring = 'n';

if exist('shift_vector.mat','file')
    load shift_vector
else
    
    dd = floor(dims(2)/128);
    dimg = abs(img(1:dd:end,1:dd:end,1:dd:end,1));
    cmask = cmask(1:dd:end,1:dd:end,1:dd:end,1);
    imax = max(dimg(:))/2;

    dims = size(dimg);

    mc = ceil(dims(1)/2);
    nc = ceil(dims(2)/2);
    zc = ceil(dims(3)/2);
    mtotal = 0;
    ntotal = 0;
    ztotal = 0;
        
while strcmp(cstring,'y') ~= 1

    figure(91)
    show(dimg(:,:,zc)+cmask(:,:,zc)*imax);
    [n,m] = ginput(1);
    n = floor(n);
    m = floor(m);
    dimg = circshift(dimg,[mc-m nc-n 0]);

    mtotal = mtotal+mc-m;
    ntotal = ntotal+nc-n;

figure(90);
    subplot(311); show(squeeze(dimg(mc,:,:))'); 
    subplot(312); show(squeeze(dimg(:,nc,:))');
    subplot(313); show(squeeze(dimg(:,:,zc)));
    
    
    figure(91)
    show(squeeze(dimg(:,nc,:)+cmask(:,nc,:)*imax));
    [z,m] = ginput(1);
    z = floor(z);
    m = floor(m);
    dimg = circshift(dimg,[mc-m 0 zc-z]);

    mtotal = mtotal+mc-m;
    ztotal = ztotal+zc-z;

figure(90);
    subplot(311); show(squeeze(dimg(mc,:,:))'); 
    subplot(312); show(squeeze(dimg(:,nc,:))');
    subplot(313); show(squeeze(dimg(:,:,zc)));
    
    figure(91)
    show(squeeze(dimg(mc,:,:)+cmask(mc,:,:)*imax)');
    [n,z] = ginput(1);
    n = floor(n);
    z = floor(z);
    dimg = circshift(dimg,[0 nc-n zc-z]);
        
    ntotal = ntotal+nc-n;
    ztotal = ztotal+zc-z;

figure(90);
    subplot(311); show(squeeze(dimg(mc,:,:))'); 
    subplot(312); show(squeeze(dimg(:,nc,:))');
    subplot(313); show(squeeze(dimg(:,:,zc)));



    cstring = waitinput('   Is this alignment okay? (y/n): ',15,'s');
    if isnan(cstring)
        fprintf('\n   No response, I''ll assume you meant yes.\n');
        cstring = 'y';
    end
    

mtotal = mtotal+mc-m;
ntotal = ntotal+nc-n;
ztotal = ztotal+zc-z;
    
end

shift_vector = dd*[mtotal ntotal ztotal];
save shift_vector shift_vector
end

fprintf('   Shifting image center...\n');
% for k = 1:nechoes
%     img(:,:,:,k) = circshift(img(:,:,:,k),shift_vector);
% end

img = circshift(img,[shift_vector 0 0]);

figure(90);
subplot(311); show(squeeze(abs(img(mc,:,:,1,1)))'); 
subplot(312); show(squeeze(abs(img(:,nc,:,1,1)))');
subplot(313); show(squeeze(abs(img(:,:,zc,1,1))));


end

