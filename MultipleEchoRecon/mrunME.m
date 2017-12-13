% mrunME.m

runfolder = 'Run';
folders = dir([runfolder '*']);
iBIAC = 221;
orientation = [1 -1 -1];

for r = 1:length(folders)
    cd(folders(r).name)
    Pconvert('');
    if r == 1
        load dims
        load hdr
        load TE
        load vox
        load fov
        save('../dims.mat','dims');
        save('../hdr.mat','hdr');
        save('../TE.mat','TE');
        save('../vox.mat','vox');
        save('../fov.mat','fov');
    end
    cd ..
end

save orientation orientation
raw = sumMEraw('k','Run','short',1:dims(4),orientation);
writeMEraw(raw,'kmean','','single');
% raw = Pseries5('kmean','single','img','shift','mask',['BIAC' num2str(iBIAC)],'matsave');
raw = Pseries5('kmean','single','img','matsave');
show3(abs(raw));