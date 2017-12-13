function [img5D] = niiME5D(img_type,dir_name,option)
%%  niiME5D.m: Generates a 5D .nii file that can be read into ImageJ
%   niiME5D(img_type,dir_name)
%   img_type:   Multi-echo image type ('mag','frq',etc.) corresponding to
%               Russell Dibb's ME Toolbox naming convention
%   dir_name:   Base directory name holding image files from each channel
%
%   Example:    niiME5D('mag','ProHance');

dirs = dir([dir_name '*']);
channels = length(dirs);
echo_files = dir([dirs(1).name '/echo*' img_type '*']);  
echoes = length(echo_files);
ref_nii = load_nii(fullfile(dirs(1).name,echo_files(1).name));
pixdim = ref_nii.hdr.dime.pixdim(2:4);
img5D = zeros([size(ref_nii.img) echoes channels]);
back = 0; time = 0;

for k = 1:channels
    [back,time] = progress(k,length(dirs),'Reading channel',back,time);
    cdir = dirs(k).name;   
    echo_files = dir([cdir '/echo*' img_type '*']);  
    echoes = length(echo_files);
    for q = 1:echoes
        nii = load_nii(fullfile(cdir,echo_files(q).name));
        img5D(:,:,:,q,k) = nii.img;
    end
end

filename = ['5D_' dir_name '.' img_type '.nii'];

if ~exist('option','var')
    option = '';
end
if strcmp(option,'avg')
    filename = ['cumulative_' filename];
    back = 0; time = 0;
    for k = 1:channels
        [back,time] = progress(k,length(dirs), ...
                      'Generating cumulative image for channel',back,time);
        for q = 1:echoes
            img5D(:,:,:,q,k) = mean(img5D(:,:,:,1:q,k),4);
        end
    end
end

fprintf('   Converting file type to facilitate visualization...\n');
img5D = img5D/max(abs(img5D(:)))*(2^15-1);
fprintf('   Saving signed, 16-bit, 5D hyperstack as %s...\n',filename);
nii5D = make_nii(img5D,pixdim,[0 0 0],4);
save_nii(nii5D,filename);

end



