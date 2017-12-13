function img = centerME(img)
% centerME.m: Helps you center an image volume (3D or 4D) with GUI.

if ~exist('dims','var')
    load dims;
end

mc = ceil(dims(1)/2);
nc = ceil(dims(2)/2);
zc = ceil(dims(3)/2);
nechoes = dims(4);
mtotal = 0;
ntotal = 0;
ztotal = 0;

dd = floor(dims(2)/128);
dimg = abs(img(1:dd:end,1:dd:end,1:dd:end,1));

figure(90);
subplot(311); show(squeeze(abs(img(mc,:,:,1)))'); 
subplot(312); show(squeeze(abs(img(:,nc,:,1)))');
subplot(313); show(squeeze(abs(img(:,:,zc,1))));
cstring = 'n';

if exist('shift_vector.mat','file')
    fprintf('   Shifting image center...\n');
    load shift_vector
    for k = 1:nechoes
        img(:,:,:,k) = circshift(img(:,:,:,k),shift_vector);
    end
    figure(90);
    subplot(311); show(squeeze(abs(img(mc,:,:,1)))'); 
    subplot(312); show(squeeze(abs(img(:,nc,:,1)))');
    subplot(313); show(squeeze(abs(img(:,:,zc,1))));
else
while strcmp(cstring,'y') ~= 1

    figure(91)
    show(img(:,:,zc,1));
    [n,m] = ginput(1);
    n = floor(n);
    m = floor(m);
    for k = 1:nechoes
        img(:,:,:,k) = circshift(img(:,:,:,k),[mc-m nc-n 0]);
    end
    mtotal = mtotal+mc-m;
    ntotal = ntotal+nc-n;

figure(90);
subplot(311); show(squeeze(abs(img(mc,:,:,1)))'); 
subplot(312); show(squeeze(abs(img(:,nc,:,1)))');
subplot(313); show(squeeze(abs(img(:,:,zc,1))));
    
    
    figure(91)
    show(squeeze(img(:,nc,:,1)));
    [z,m] = ginput(1);
    z = floor(z);
    m = floor(m);
    for k = 1:nechoes
        img(:,:,:,k) = circshift(img(:,:,:,k),[mc-m 0 zc-z]);
    end
    mtotal = mtotal+mc-m;
    ztotal = ztotal+zc-z;

figure(90);
subplot(311); show(squeeze(abs(img(mc,:,:,1)))'); 
subplot(312); show(squeeze(abs(img(:,nc,:,1)))');
subplot(313); show(squeeze(abs(img(:,:,zc,1))));
    
    figure(91)
    show(squeeze(img(mc,:,:,1))');
    [n,z] = ginput(1);
    n = floor(n);
    z = floor(z);
    for k = 1:nechoes
        img(:,:,:,k) = circshift(img(:,:,:,k),[0 nc-n zc-z]);
    end
    ntotal = ntotal+nc-n;
    ztotal = ztotal+zc-z;

figure(90);
subplot(311); show(squeeze(abs(img(mc,:,:,1)))'); 
subplot(312); show(squeeze(abs(img(:,nc,:,1)))');
subplot(313); show(squeeze(abs(img(:,:,zc,1))));



    cstring = waitinput('   Is this alignment okay? (y/n): ',15,'s');
    if isnan(cstring)
        fprintf('\n   No response, I''ll assume you meant yes.\n');
        cstring = 'y';
    end
    

mtotal = mtotal+mc-m;
ntotal = ntotal+nc-n;
ztotal = ztotal+zc-z;
    
end

shift_vector = [mtotal ntotal ztotal];
save shift_vector shift_vector
end




end

