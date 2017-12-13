function  raw = Pseries5(varargin)
%%function raw = Pseries(varargin)
% Before calling the function, put the pfiles (single- or multi-echo) in a
% separate folder. Call the function in the created folder. 
%
%   
% varargin is a list of options in 'single quotes' separated by commas
%   (no input)  - output is complex, raw, frequency space data
%   'img'       - output is complex data in image space
%   'echofill'  - replaces faulty data in last echo when phase1*phase2 is
%                 not divisible by baseline_spacing (usually 2048)
%   'shift'     - linearly corrects even echo shift to match odd echoes
%   [-1 -1 -1]  - flips readout, phase1, and phase2 array orientations
%                 (good for getting data into WHS if you scanned it weird)
%   'backphase' - removes background phase, only corrects readout direction
%                 (more coming soon)
%   'mask'      - creates a binary, skull-stripped mask
%   'BIAC#'     - saves data in format for susceptibility computation by
%                 the BIAC cluster
%                   
%
% Example: imagedata = Pseries4('img',[-1,1,-1]);
% imagedata will be an image space volume of complex data with the readout
% and phase2 directions reversed.
%
% img = Pseries4('img','shift','mask');
%
%
% Based on the function m-file, MapFreq.m, created by Chunlei Liu
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ORIG Input:
% function [Freq] = MapFreq()
% Input:     The function will read in all pfiles from the current
%            directory.
% Output:
% Freq  -    Frequency maps in Hz, stored in binary float format 
%            nrows x ncols x nslices with the file name of echo#.freq
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Chunlei Liu, Duke University, 06/2009
%
% slg - 100715 V2 removed Dim as param MapFreq, use
%       dimensions from pfile header (as stored by civm 3D cartesian scans)
%       Changed raw data reader to use header dimensions and bl skip.
%       V3 added fft2c and fftnc and inverses, provided by Chunlei.
%       Wrapped GE header reader to notice missing Pfile.
% slg - 100716 V4 Now able to discern raw format and read short and int/EDR 
%       Pfiles.
% rmd - 111025 V5 Now mimics the effect of radish reconstruction to pull
%       data from multiple, multi-echo pfiles and sort it into echoes
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% User-defined parameters
SHORT_SIZE = 2; % in bytes
for narg = 1:nargin
    if ischar(varargin{narg})
        if strcmp(varargin{narg},'theworks')
            imgswitch = 1;
            shiftswitch = 1;
            phaseswitch = 1;
            maskswitch = 1;
        elseif strcmp(varargin{narg},'nomask')
            imgswitch = 1;
            shiftswitch = 1;
            phaseswitch = 1;
        elseif strcmp(varargin{narg},'img')
            imgswitch = 1;
        elseif strcmp(varargin{narg},'rp')
            rp_flag = 1;
            rp_file = dir('*.rp');
        elseif strcmp(varargin{narg},'echofill') 
            fillswitch = 1;
        elseif strcmp(varargin{narg},'shift')
            shiftswitch = 1;
        elseif strcmp(varargin{narg},'phaseshift')
            phaseswitch = 1;
        elseif strcmp(varargin{narg},'mask')
            maskswitch = 1;
        elseif strcmp(varargin{narg},'naxos')
            ram = 32e9; 
        elseif strcmp(varargin{narg},'andros')
            ram = 64e9;
        elseif strcmp(varargin{narg},'center')
            centerswitch = 1;
            imgswitch = 1;
        elseif strcmp(varargin{narg},'single')
            pfile_type = 'single';
        elseif strcmp(varargin{narg},'int')
            pfile_type = 'int';
        elseif strcmp(varargin{narg},'binsave')
            binsave = 1;
        elseif strcmp(varargin{narg},'matsave')
            matsave = 1;
        elseif strcmp(varargin{narg},'short')
            pfile_type = 'short';    
        elseif strcmp(varargin{narg},'convert')
            SHORT_SIZE = 2; % in bytes
            convertflag = 1;
        elseif length(varargin{narg}) < 5
            scan_name = varargin{narg};
        elseif ~strcmp(varargin{narg}(1:4),'BIAC')
            scan_name = varargin{narg};
        elseif strcmp(varargin{narg}(1:4),'BIAC')
            BIACswitch = 1;
            maskswitch = 1;
            iBIAC = str2num(varargin{narg}(5:end));
        end
    else
        if length(varargin{narg}) == 3
            
            if length(varargin{narg}(varargin{narg} > 1)) > 1
                permuteswitch = 1;
                permute_vector = varargin{narg};
            else
                orientswitch = 1;
                orientation = varargin{narg};
            end
        end
    end   
end


%% Convert p-file series to echo files

if ~exist('scan_name','var')
    scan_name = '';
end

if exist('convertflag','var')
    Pconvert(scan_name);
    scan_name = [scan_name 'k'];
%     pfile_type = 'short';
elseif exist('rp_flag','var')
    [hdr,header_bytes] = WrapGEheader(rp_file(1).name);
    bp_half_element = hdr.rdb.point_size;
    baseline_spacing = hdr.rdb.nframes;
    SHORT_SIZE = 2;  
    flag='l';  % default byte swap
    if bp_half_element == SHORT_SIZE
        pfile_type = 'short';
    else
        if bp_half_element ~= 4
            error('\nPfile rdb.point_size unknown.\n');
        end
        pfile_type = 'int';
    end
    % Save scan paramters
    dims = [hdr.rdb.da_xres, hdr.rdb.user7, ...
        hdr.rdb.user8, hdr.rdb.nechoes]; save dims dims
    fov = [hdr.rdb.user16, hdr.rdb.user17, hdr.rdb.user18]; save fov fov
    vox = fov./dims(1:3); save vox vox   % voxel size, mm
    if  ((hdr.rdb.user33 == 0) && dims(4) > 1)
        hdr.rdb.user33 = waitinput('What is the echo spacing (in ms)? ',10);
        if isnan(hdr.rdb.user33)
            hdr.rdb.user33 = 2.896;
        end
    end
    echo_spacing = hdr.rdb.user33/1e3; % Only valid in the PSD Russell uses
    TE = (0:dims(4)-1)*echo_spacing+hdr.rdb.te/1e6; save TE TE
    save hdr hdr header_bytes

    raw=zeros(dims(1),dims(2),dims(3));
    back = 0; time = 0;
    fid = fopen(rp_file(1).name,'r','l');
    skip=fread(fid,[1 header_bytes/SHORT_SIZE],'short');  % header
    for z=1:dims(3)
        [back,time] = progress(z,dims(3),'Reading slice',back,time);
            for y=1:dims(2)    
                dump=fread(fid,[2 dims(1)], pfile_type);
                raw(:,y,z)=squeeze(dump(2,1:end)+1i*dump(1,:));
            end
    end
    fclose(fid);   
end


load dims;
load hdr hdr;


%% Read in raw data

% % This automatically converts the data to single if it is too large for
% % particular computer to handle.f
%
% databytes = dims(1)*dims(2)*dims(3)*dims(4)*2;
% doubleram = databytes*4;
% 
% if ~exist('threshbytes','var')
%     ram = 6e9; % Jeeves/default
% end
% 
% if doubleram > ram
%     singleflag = 1;    
% end

if dims(4) > length(dir([scan_name '*']))
    fprintf(['   Number of echo files does not match dims.mat.\n' ...
             '   Reading only the number of available echo files.\n']);
end

if ~exist('rp_flag','var')
    % SHORT_SIZE = 2;  % in bytes
    if ~exist('pfile_type','var')
        if hdr.rdb.point_size == SHORT_SIZE
            pfile_type = 'short';
        elseif SHORT_SIZE == 4
            pfile_type = 'single';
        else
            if hdr.rdb.point_size ~= 4
                error('\nPfile rdb.point_size unknown.\n');
            end
            pfile_type = 'int';
        end
    end
    raw = readMEraw(scan_name,'',pfile_type, 1:length(dir([scan_name '*'])));
end

if exist('singleflag','var')
    fprintf('   Converting to single...\n');
    raw = single(raw);
end

if exist('orientswitch','var')
    raw = orientME(raw,orientation);
    save orientation orientation
end

if exist('permuteswitch','var')
    raw = permute(raw,permute_vector);
    save permute_vector permute_vector
end

% Only phase corrects readout direction
if exist('phaseswitch','var')
    fprintf('  Performing odd/even phase correction in readout direction...\n');    
    raw = phaseME(raw);
end

if exist('shiftswitch','var')
    fprintf('   Correcting even echo shift...\n');
    bvec = shiftME(raw);
    for q = 2:2:dims(4)
        for z = 1:dims(3)
            for y = 1:dims(2)
                raw(:,y,z,q) = raw(:,y,z,q).*bvec;   
            end
        end
    end 
end

% Perform IFFT
if exist('imgswitch','var')
    back = 0; time = 0;
    for k = 1:dims(4)
        [back,time] = progress(k,dims(4),'Performing ifft on echo',back,time);
        raw(:,:,:,k) = ifftnc(raw(:,:,:,k));
    end
    figure(71); 
    subplot(211); montageME(abs(raw),ceil(dims(3)/2),'1D',dims(1:3));
    subplot(212); montageME(angle(raw),ceil(dims(3)/2),'1D',dims(1:3));
    data_str = 'img';
    if exist('centerswitch','var')
        raw = centerME(raw);
        figure(1); montageME(abs(raw),ceil(dims(3)/2),'1D',dims(1:3));
    end
else
    data_str = 'raw';
end

if exist('maskswitch','var')
    [~] = brainext(raw);
    if exist('BIACswitch','var')
        fprintf('   Saving mask and phase information for BIAC computations...\n');
        BIACprep(raw,iBIAC)
    end
end
if exist('binsave','var')
    fprintf('   Saving binary echo files...\n');
    writeMEraw(raw,data_str,'','single');
end

if exist('matsave','var')
    fprintf('   Saving data as img.mat...\n');
    save -V7.3 img raw
end

fprintf('   P-file read complete!\n');
end
