function [mag,frq] = BrukerReconME(bruker)
%%  Original file by Evan Calabrese, 6/8/11
%   rmd - Adapted for multi-echo read, 6/10/11
%
%   Script requires use of Evan's other m-files: readBrukerDirectory.m,
%   readBrukerHeader.m, and readBrukerFID.m.  Required functions from
%   Russell Dibb's ME Toolkit should all be included at the end of the main
%   function.

%% Select Directories

if isempty(bruker)
   display('   Select directory containing fid file to recon.')
   [bruker,dirpath] = readBrukerDirectory();
   fprintf('   Input: %s\n',dirpath);
end

display('   Select destination directory for raw multi-echo file output.')
save_path=uigetdir('/~', 'select save path');
fprintf('   Output: %s\n',[save_path '/']);

% Choose if you want to wait for it to calculate frequency maps (ME only)
if nargout < 2
    freqoption = input('   Calculate frequency maps as well? (y/n): ','s');
else
    freqoption = 'y';
end

%% Read matrix size

m = bruker.method.PVM_EncMatrix(1);             % Frequency encode direction
n = bruker.method.PVM_EncMatrix(2);             % Phase encode 1 direction
if length(bruker.method.PVM_EncMatrix) < 3
    p = 1;
else
    p = bruker.method.PVM_EncMatrix(3);         % Phase encode 2 direction
end
echoes = bruker.method.PVM_NEchoImages;         % Number of echoes
nrecs = bruker.method.PVM_EncNReceivers;        % Number of receivers

if (strcmp(freqoption,'y') == 1 || strcmp(freqoption,'Y'))
    options.ncoil = nrecs;
    options.fermir = 8/512;
    options.fermiw = 8/512; % Change this variable to adjust unwrapping
end

%% Calculate magnitude and frequency map images

if echoes > 1    
    data_mat=reshape(bruker.fid,m,nrecs,echoes,n,p);
    data_mat = permute(data_mat,[1 4 5 3 2]);
    back = 0; time = 0;
    for k = 1:echoes
        [back,time] = progress(k,echoes,'Generating magnitude image',back,time);
        for q = 1:nrecs
            data_mat(:,:,:,k,q) = fftshift(ifftn(fftshift(data_mat(:,:,:,k,q))));
        end
        if (strcmp(bruker.method.EchoAcqMode,'allEchoes') == 1 && mod(k,2) == 0)
            data_mat(:,:,:,k,:) = data_mat(end:-1:1,:,:,k,:);
        end
    end
    
%     figure(85); % Plot data for each receiver
%     for k = 1:nrecs
%         subplot(nrecs*2,1,k); montageME(abs(data_mat(:,:,:,:,k)),ceil(p/2),'1D');
%         title(['Magnitude images for coil ' num2str(k)]);
%         subplot(nrecs*2,1,k+nrecs); montageME(angle(data_mat(:,:,:,:,k)),ceil(p/2),'1D');
%         title(['Angle map for coil ' num2str(k)]);
%     end

    mag = abs(sum(data_mat,5));
    figure(1); % Display
    montageME(mag,ceil(p/2)); title('Magnitude image echoes');

    
    if (strcmp(freqoption,'y') == 1 || strcmp(freqoption,'Y')) % Frequency maps
        frq = zeros(size(mag));
        back = 0; time = 0;
        for k = 1:echoes
            [back,time] = progress(k,echoes,'Generating frequency map',back,time);
            frq(:,:,:,k) = CalcPhaseBruker(squeeze(data_mat(:,:,:,k,:)),options)/pi;
        end
        figure(1); % Display echo montages
        if echoes <= 10
            dispstr = '1D';
        else
            dispstr = '';
        end
        subplot(211); montageME(mag,ceil(p/2),dispstr); 
        title('Magnitude image echoes');
        subplot(212); montageME(frq,ceil(p/2),dispstr); 
        title('Frequency map echoes');
     
    end
else % Unproven code for magnitude and frequency images for single echo data
%     m = 256; n = 360; p = 256;
    data_mat=squeeze(reshape(bruker.fid,m,nrecs,n,p));
    data_mat = permute(data_mat,[1 3 4 2]);
    for q = 1:nrecs
        data_mat(:,:,:,q) = fftshift(ifftn(fftshift(data_mat(:,:,:,q))));
    end
    mag = abs(sum(data_mat,4));
    m=max(mag(:));
    mag=mag/m*(2^15-1); %for converting to 16 bit  
%     nii=make_nii(mag, bruker.method.PVM_SpatResol, zeros(1,length(bruker.method.PVM_Matrix)), 512);
%     out_file=strcat(bruker.subject.SUBJECT_name_string,'_',bruker.subject.SUBJECT_date,'.nii');
%     save_nii(nii,fullfile(save_path,out_file));
    if (strcmp(freqoption,'y') == 1 || strcmp(freqoption,'Y'))
        frq = CalcPhaseBruker(data_mat,options)/2/pi/bruker.method.PVM_EchoTime1;
    end
end

% writeME(mag,[save_path '/']);
if exist('frq','var')
%     writeME(frq,[save_path '/']);
end
fprintf('   Process complete!\n');

end

function [img_out] = CalcPhaseBruker(img,options)
% Compute the phase map at 7T (initially used for GE Signa 7T)
% The phase map is computed by removing the coil phase. The coil phase
% is estimated using a low resolution image. Assume original images are 
% stored in the same directory. Odd-numbered images are magnitude image;
% even-numbered images are phase image (x1000). The last two images for 
% each slice are the sum-of-squares images.
% options.dim "3D", "2D"
% options.ncoil
% options.fermir: number of voxel per 512, e,g. 8/512
% options.fermiw: number of voxel per 512, e.g. 8/512
% Chunlei Liu, Stanford University, 08/18/2007,
%              Duke University, 2009
% rmd - Adapted for Russell Dibb's Multi-echo Toolkit, 05/31/2011
%              Duke Universtiy, 2011
% rmd - Adapted for use with Bruker raw data, 06/10/2011

[opxres,opyres,opzres,ncoil] = size(img);

%%% Fermi Window
[y,x,z] = meshgrid(-opyres/2:opyres/2-1,-opxres/2:opxres/2-1,-opzres/2:opzres/2-1);
if (~isfield(options,'fermir'))
    fermir=8*opxres/512;%100,12,8,32; 4 for 3T
else
    fermir=options.fermir*opxres;
end

if (~isfield(options,'fermiw'))
    fermiw=8*opxres/512;%12,8, 32;
else
    fermiw=options.fermiw*opxres;
end
krad = sqrt(x.^2+y.^2+z.^2);
FW = 1./(1+exp((krad-fermir)/fermiw));

img_out = zeros(opxres,opyres,opzres,ncoil);
%%% compute phase
for icoil = 1:ncoil
    ipha_low = ifftnc(fftnc(img(:,:,:,icoil)).*FW);
    img_out(:,:,:,icoil) = abs(img(:,:,:,icoil)).*exp(1i*(angle(img(:,:,:,icoil)) - angle(ipha_low)));
end
img_out = angle(sum(img_out,4));

end

function im = fftnc(d)
% im = fftnc(d)
%
% fftnc perfomrs a centered fftn
%
im = fftshift(fftn(fftshift(d)));
end

function im = ifftnc(d)
% im = ifftnc(d)
%
% ifftnc perfomrs a centered ifftn
%
im = fftshift(ifftn(fftshift(d)));
end

function [f] = montageME(img,slice_num,option,dims)
%montageME: Automates the generation of a montage of echoes.
%   Detailed explanation goes here
%   img: either a 4-D array (m,n,p,echoes) or a 3-letter, ME string code
%   slice_num: the slice that you want
%   option 'write': writes the echoes of the slice to an image stack
%   option '1D': displays the montage in a 1 x echoes array
%   dims: 1x1 or 1x2 if not square.  Only required if "img" is a string.

if ~exist('option','var')
    option = '';
end

if ischar(img) == 1 % Execute these lines if reading raw files
    
    % Establish variable names and dimensions
    img_type = img;
    m = dims(1);
    img = dir(['echo*' img_type]);
    if length(dims) > 1
        n = dims(2);
    else
        n = m;
    end
    
    % Read raw files
    if isempty(img)
        fprintf('   There are no .%s files in the folder!\n', img_type);
    else
        echoes = length(img);
        f = zeros(m,n,1,echoes);
        for k = 1:echoes
            fid = fopen(img(k).name,'r','l');
            fread(fid,m*n*(slice_num-1),'float');   % Skip until desired slice
            f(:,:,1,k) = fread(fid,[dims(1) dims(2)],'float');
            fclose(fid);
        end               
    end       
else    % Execute these lines if input data is 4D matrix
    f = img(:,:,slice_num,:);
    img_type = inputname(1);
end

if exist('f','var') == 1
    % Prep data for montage
    f = double(f);
    fmax = max(f(:));
    fmin = min(f(:));

    % Display montage, 1D format if designated in options
    if strcmp(option,'1D') == 1
        montage(f,'DisplayRange',[fmin fmax],'Size',[1 size(f,4)]);
    elseif strcmp(option,'write') == 1
        if ~isempty(findall(0,'Type','Figure')) == 1
            h = gcf;
            close(h);
            figure(h);
        else
            figure(1)
            h = 1;
        end
        montage(f,'DisplayRange',[fmin fmax]);
    else
        montage(f,'DisplayRange',[fmin fmax]);
    end
    f = squeeze(f);
    
    % Write image if designated as an option
    if strcmp(option,'write') == 1
        filename = ['slice' num2str(slice_num) '.' img_type];
        fid = fopen(filename,'w');
        fwrite(fid,f,'float');
        fclose(fid);
        fprintf('   Wrote 32-bit real, %s\n',filename);        
        saveas(h,['slice' num2str(slice_num) '.tif']);
    end
end

end

function [] = writeME(img,save_path)
% writeME.m: Writes a series of raw data files from multi-echo data using
% script writeRaw.m
%
%   Make sure to use the proper image type name (e.g. mag, frq, phs, scl, sus) 
%   because that will become the raw data file extension.  Function works
%   for echo trains up to 99 echoes long.
%   
%   writeME(img,index)
%   img = 4-D array (m,n,p,echoes)
%   index = echo index used for writing

img_type = inputname(1);
if ~exist('save_path','var')
    save_path = '';
end

for k = 1:size(img,4)
    if k < 10
        imgName = ['echo0' num2str(k) '.' img_type];
    else
        imgName = ['echo' num2str(k) '.' img_type];
    end
    writeRaw(img(:,:,:,k),imgName,save_path);    
end

end

function [] = writeRaw(img,name,save_path)
% writeRaw.m: Writes a 32-bit real (float), raw data file from provided 
%             data and name.  Save path is optional.

if ~exist('save_path','var')
    save_path = '';
end

fid = fopen([save_path name],'w');
fwrite(fid,img,'float');
fclose(fid);
fprintf('   Wrote 32-bit real, %s\n',name);
    
end

function [new_backspaces,time] = progress(index,total,title,backspaces,time)
%% [new_backspaces,time] = progress(index,total,title,backspaces,time)

if index > 1   
    time(index) = toc;
%     tdisp = ceil((total-index+1)*time(index));
    if index > 12
        tdisp = ceil((total-index+1)*median(time(index-10:index)));
    else
        tdisp = ceil((total-index+1)*median(time(3:index)));
    end
    if tdisp == 1
        str3 = 'second';
    else
        str3 = 'seconds';
    end
    for k = 1:backspaces
        fprintf('\b');
    end
    
    if tdisp > 59
        tm = floor(tdisp/60);
        if tm == 1
            str4 = 'minute';
        else
            str4 = 'minutes';
        end
        tdisp = mod(tdisp,60);
        str4 = [num2str(tm) ' ' str4 ', '];
    else
        str4 = '';
    end   
else
    time = zeros(1,total);
end

tic
str1 = ['   ' title ': ' num2str(index) '/' num2str(total)];
fprintf('%s',str1);

if (index == total && index > 2)
    str2 = [' - ' str4 num2str(tdisp) ' ' str3 ' remaining...'];
    fprintf('%s',str2);
    fprintf('\n');
elseif index > 2
    str2 = [' - about ' str4 num2str(tdisp) ' ' str3 ' remaining...'];
    fprintf('%s',str2);
elseif total > index
    str2 = ' - estimating time to completion...';
    fprintf('%s',str2);
else
    str2 = '';
    fprintf('\n');
end

new_backspaces = length(str1)+length(str2);

end


