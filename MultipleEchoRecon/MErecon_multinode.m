
function MErecon_multinode(workpath,scanner,runno,study,series,auto_shift_flag)
pause(20);
if exist([workpath '/' runno 'recon.mat'],'file') ~= 2
    agilent2glusterspace_multinode(scanner,runno,study,series,workpath);
end

if exist([workpath '/' runno 'recon.mat'],'file') ~= 2
    error('missing data')
end

cd(workpath);

%% Handle some variables
load([workpath '/' runno 'recon.mat'], '-regexp', '^(?!workpath)\w');

vox = res; 
if isfield(procpar,'TE')
    TE = procpar.TE/1000; 
else
    TE = procpar.te;   
end
dims = [voldims nechoes];
if dims(3) == 0
    dims(3) = 1;
end

save([workpath '/TE.mat'],'TE');
save([workpath '/vox.mat'],'vox');
save([workpath '/fov.mat'],'fov');
save([workpath '/dims.mat'],'dims');
save([workpath '/dataprops.mat'],'TE','vox','fov','dims');

%%
if ((exist([workpath '/img.nii'],'file') ~= 2) && (exist([workpath '/img.nii.gz'],'file') ~= 2))
    
%%% Sort multiecho data and save
fidvec = zeros(1,nechoes);
for k = 1:nechoes
    fidvec(k) = fopen([workpath '/k_single_echo_' num2str(k) '.raw'],'w+');
end
sort_fid(fidpath,max_blocks*nechoes,ntraces,npoints,bitdepth,1,fidvec,nechoes);
for k = 1:nechoes
    fclose(fidvec(k));
end

raw = readMEraw('k',workpath,'single',1:nechoes);
% raw = phase_ramp_remove1(raw);
raw = iftME(raw);
if auto_shift_flag == 1
    if exist([workpath '/shift_vector.mat'],'file') == 2
        raw = circshift(raw,[shift_vector 0 0]);
    else
        [raw,shift_vector] = autoC(raw); 
        save([workpath '/shift_vector.mat'],'shift_vector');
    end
end
save_nii(make_nii(raw,vox,[0 0 0],32),[workpath '/img.nii']);
% complex_nii(raw,vox,[workpath '/img.nii']);
delete([workpath '/k_single_echo_*.raw']);
end

save_nii(make_nii(abs(open_nii([workpath '/img.nii'],1)),vox,[0 0 0],16),[workpath '/mag1.nii']);

display('Complex-valued images reconstructed and saved as img.nii');
end


