function procpar2reconvars(workpath,scanner,runno,study,series)
%% Converts procpar file to recon.mat in same directory, populated with recon variables.
% The only required input is workpath. All the others are included for
% reference only, and are not necessary for reconstruction.
%% read headers collect metadata and set up work directories

%consider adding the metadata GUI here
%first assemble the path to the data
fidpath=dir([workpath '/*fid']);
fidpath = [workpath '/' fidpath.name];
procpar_path= dir([workpath '/*procpar']);
procpar_path = [workpath '/' procpar_path.name];
procpar = readprocparSTI(procpar_path);
[~,runno] = fileparts(procpar_path);

%read the fid file header
[npoints,nblocks,ntraces,bitdepth] = load_fid_hdr(fidpath);
display(['fid file has ' num2str(npoints) ' points; ' num2str(ntraces) ' traces; ' num2str(nblocks) ' blocks; and bitdepth ' bitdepth]);

%% memory checks
%check system memory
display('Checking system memory and purging, this will take a few seconds');
% system('purge');
meminfo=imaqmem; %check available memory
% while meminfo.TotalPhys/2>meminfo.AvailPhys
%     input('Less than half of physical memory is currently available, close some programs and then press enter >> ','s');
% end
divisor=3; %this is the expected multiplication in memory, reduce for efficient memory use
max_total_points=meminfo.AvailPhys/(4*divisor); %calculate number of single percision points we can put in memory
% handle ignore memory limit options
% if ignore_memory_boolean==1;
%     display('you have chosen to ignore this machine''s memory limits, this machine may crash');
%     max_total_points=npoints*ntraces*nblocks;
% end
save([workpath '/procpar.mat'],'procpar');
voldims=[procpar.np/2 procpar.nv procpar.nv2];
if voldims(3) == 0;
    voldims(3) = 1;
end
nvols=(npoints/2*ntraces*nblocks)/prod(voldims);
blocks_per_vol=nblocks/nvols;
fov=[procpar.lro procpar.lpe procpar.lpe2].*10; %fov in mm this may not be right for multislice data
res=fov./voldims;

if isfield(procpar,'TE')
    nechoes = numel(procpar.TE);
else
    nechoes = 1;
end

%check to see if we need to do this in chunks or not
if nvols>1 %recon one volume at a time
    numchunks=nvols;
    max_blocks=blocks_per_vol;
%     if max_blocks*npoints*ntraces>max_total_points
%         error('volume size is too large, consider closing programs and restarting')
%     end
else %if its just one volume, see if we can do it all at once or need to do chunks
    max_blocks=floor(max_total_points/(ntraces*npoints)); %number of blocks we can work on at a time
    numchunks=ceil(nblocks/max_blocks);
end
%Now we have all the information we need to decide what type of recon to do
save([workpath '/' runno 'recon.mat']);

display('Agilent files pulled successfully');

