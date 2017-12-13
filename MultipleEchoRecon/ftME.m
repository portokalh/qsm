function raw = ftME(raw)
% Converts raw k-space into raw image for 4D datasets

dims = size(raw);

back = 0; time = 0;
if length(dims) == 3
    for k = 1:dims(3)
        raw(:,:,k) = fftnc(raw(:,:,k));
    end    
    figure(90); montageME(abs(permute(raw,[1 2 4 3])),1,'1D',dims(1:2));
else
    for k = 1:dims(4)
        [back,time] = progress(k,dims(4),'Performing fft on echo',back,time);
        raw(:,:,:,k) = fftnc(raw(:,:,:,k));
    end
    figure(90); montageME(abs(raw),ceil(dims(3)/2),'1D',dims(1:3));
end

end

function montageME(img,slice_num,option,dims)
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
if ~exist('dims','var')
    if ~ischar(img)
        dims = size(img);
    else
        load dims;
    end
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
    % Prep data for montage-ing
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