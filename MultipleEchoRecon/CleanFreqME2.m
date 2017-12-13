% CleanFreqME2.m
filetype = 'mag';
directory = pwd;
files = {'echo01.mag.nii','echo02.mag.nii'};


for q = 1:1%length(files);
    filename = strcat(fullfile(directory,files{q}));
    msk = strip_mask_T2w(filename,1,3);
    for k = q:2:echoes
        if k < 10
            zstr = '0';
        else
            zstr = '';
        end
        echofile = strcat(filename(1:end-9),zstr,num2str(k),'.',filetype,'.nii');
        nii = load_nii(echofile);
        nii.img = msk.*nii.img;
        save_nii(nii,fullfile(directory,echofile))
    end
end
