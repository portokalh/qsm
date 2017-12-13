function  [mag,ang,TE] = MapAngleME(echo_spacing,extbrain)
%%function [Freq] = MapFreqME(pfile)
% Before calling the function, put the pfiles of GRE sequences in a
% separate folder. Call the function in the created folder. Otherwise, you
% have to provide the full path for the pfile.
% enter MapFreqME() 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ORIG Input:
% function [Freq] = MapFreqME()
% Input:     The function will read in all pfiles with a name in 
%            the format of 'split_me_trains_gaps.*' from the current
%            directory.
% Output:
% Freq  -    Frequency maps in Hz, stored in binary float format 
%            nrows x ncols x nslices with the file name of echo#.freq
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Chunlei Liu, Duke University, 06/2009
% slg - 100715 V2 removed Dim as param MapFreq, use
%        dimensions from pfile header (as stored by civm 3D cartesian scans)
%        Changed raw data reader to use header dimensions and bl skip.
%        V3 added fft2c and fftnc and inverses, provided by Chunlei.
%        Wrapped GE header reader to notice missing Pfile.
% slg - 100716 V4 Now able to discern raw format and read short and int/EDR Pfiles.
% rmd - 110310 V_Multi_Echo - Only works for CIVM multi echo files reading
%        in split_me_trains_gaps.* files created by radish reconstruction
% rmd - 110503 V_Multi_Echo - Sets background values to zero in freq images
%        using CleanFreq.m (later became obsolete).
% rmd - 110614 V_ME - Several additiional tweaks to work with multi-echo
%        toolbox including saving dims.mat for easy readability later.
%        Also, asks for echo_spacing and view window for changing
%        filtering, though this is likely to change in the future.
% rmd - 110630 V_ME - Several more tweaks including mask generation using
%        stripMask.m to extract brain tissue only, saving data as .nii
%        files, and monitoring the appropriateness of the view window when
%        doing a phase unwrap.
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% parameters

if ~exist('extbrain','var')
    extbrain = 0;
end
if (~exist('pfile','var'))
    pfile = dir('split_me_trains_gaps.*');
    if isempty(pfile)
        pfile = dir('PN*');
        echo_spacing = 0;
    end
    if isempty(pfile)
        fprintf('   There are no pfiles in the folder!\n');
    end
    num_pfile = length(pfile);
else
    num_pfile = size(pfile,1);
    printf('   Total number of pfiles to read: %d\n',num_pfile);
    for ipfile = 1:num_pfile
        temp_pfile(ipfile).name = pfile(ipfile,:);
    end
    pfile = temp_pfile;
end
fprintf('   Total number of pfiles to read %d:\n',num_pfile);

% Properly sorts the echo train files since their numbering messes things up
if num_pfile >= 10
    for k = 1:num_pfile
        pfile(k).name = ['split_me_trains_gaps.' num2str(k)];
    end
end

%% read in raw data
options.dim = '3D';
options.ncoil = 1;
options.fermir = 8/512;
view_window = 1;
% options.fermiw_first_echo = 8/512*view_window;

for ipfile = 1:num_pfile
    fprintf('   %d/%d: %s',ipfile, num_pfile, pfile(ipfile).name);
    [hdr, header_bytes] = WrapGEheader(pfile(ipfile).name);
    
    % Initialize variables after wrapping GE header
    if ipfile == 1
        TE = zeros(1,num_pfile);
        dims = [hdr.rdb.da_xres, hdr.rdb.user7, hdr.rdb.user8]; % civm PSD 3D size info
        fov = [hdr.rdb.user16, hdr.rdb.user17, hdr.rdb.user18]; % FOV, mm
        vox = fov./dims;                                        % voxel size, mm
        mag = zeros([dims num_pfile]);
        ang = mag;
        z = ceil(size(mag,3)/2);
    end
    
    % Read and display parameters, calculate fermi window size based on TE
    TE(ipfile)  = hdr.rdb.te*1.0e-6+echo_spacing*(ipfile-1);
    fprintf(', %s %d %d %d, %s %g %s\n', 'Dim:', dims(:), 'TE:', TE(ipfile)*1000, 'ms');
    
    % Read pfile and align echoes properly
    img = ifftnc(readepicCIVM(pfile(ipfile).name));
    if (mod(ipfile,2) == 0)
        img = img(end:-1:1,:,:);
    end
    
    % Extract brain if extbrain is properly set
    if extbrain > 0
        if extbrain == 2    % Make two masks for odd and even echoes
            if ipfile == 1
                % Calculate skull-stripping mask for odd echoes
                ref_nii = 'mask_odd.nii';
                if exist(ref_nii,'file') == 2
                    msk1 = load_nii(ref_nii);
                    msk1 = double(msk1.img);
                else
                    fprintf('   Generating mask for odd echoes.\n');
                    nii = make_nii(abs(img),vox,[0 0 0],16);
                    save_nii(nii,ref_nii);
                    msk1 = stripMask(ref_nii,1);
                    save mask_odd.mat msk1
                end
            elseif ipfile == 2
                % Calculate skull-stripping mask for even echoes
                ref_nii = 'mask_even.nii';
                if exist(ref_nii,'file') == 2
                    msk2 = load_nii(ref_nii);
                    msk2 = double(msk2.img);
                else            
                    fprintf('   Generating mask for even echoes.\n');
                    nii = make_nii(abs(img),vox,[0 0 0],16);
                    save_nii(nii,ref_nii);
                    msk2 = stripMask(ref_nii,1);
                    save mask_even.mat msk2
                end
            end
            if mod(ipfile,2) == 1
                msk = msk1;
            else
                msk = msk2;
            end
        elseif (ipfile == 1 && extbrain == 1)
            % Calculate skull-stripping mask for all echoes
            ref_nii = 'mask_all.nii';
            if exist(ref_nii,'file') == 2
                msk = load_nii(ref_nii);
                msk = double(msk.img);
            else
                fprintf('   Generating mask for all echoes.\n');
                nii = make_nii(abs(img),vox,[0 0 0],16);
                save_nii(nii,ref_nii);
                msk = stripMask(ref_nii,1); 
                save mask_all.mat msk
            end
        end
        img = msk.*img;
    end
    
    % Calculate and write magnitude and phase image files
    ang(:,:,:,ipfile) = angle(img)/2/pi; %/TE(ipfile)
    mag(:,:,:,ipfile) = abs(img);    
    
    figure(1);
    subplot(211); montageME(mag,z,'1D'); title('Magnitude');
    subplot(212); montageME(ang,z,'1D'); title('Raw Angle Map');
%     if ipfile < 10
%         imgName = ['echo0' num2str(ipfile)];
%     else
%         imgName = ['echo' num2str(ipfile)];
%     end
%     mag_nii = make_nii(mag(:,:,:,ipfile),vox,[0 0 0],16);
%     save_nii(mag_nii,[imgName '.mag.nii']);
%     frq_nii = make_nii(ang(:,:,:,ipfile),vox,[0 0 0],16);
%     save_nii(frq_nii,[imgName '.frq.nii']);
    
end
if num_pfile > 1
    figure(111);
    combineME(mag);
    figure (222);
    combineME(ang,'sum');
end
save rawphase.mat ang
fprintf('   Raw angle map complete!\n');
end

%%
function raw=readepicCIVM(filename)
%function raw=readepic7T_3(filename)
% function raw=readepic7T_3(filename)
% function raw=readepic7T_2(dim1,dim2,dim3,filename,header_bytes,flag)
% Matlab function to read Signa EPIC raw data file into a complex matrix.
% Uses headfile dimension information set by CIVM PSD's.
% Discerns and handles Pfile regular and EDR formats.
% Biggest difference between this script and the regular MapFreq script is
% that echo files already have their baselines removed, so MapFreqME
% doesn't not try to remove them a second time.

[hdr, header_bytes] = WrapGEheader(filename);
dim1 = hdr.rdb.da_xres; % get raw dimensions from pfile header
dim2 = hdr.rdb.user7;
dim3 = hdr.rdb.user8;
bp_half_element = hdr.rdb.point_size;
baseline_spacing = hdr.rdb.nframes;

SHORT_SIZE = 2;  % in bytes
if bp_half_element == SHORT_SIZE
    pfile_type = 'short';
else
    if bp_half_element ~= 4
        error('\nPfile rdb.point_size unknown.\n');
    end
    pfile_type = 'int';
end
fprintf('         2x%d byte %s Pfile, nframes %d\n', bp_half_element, pfile_type, baseline_spacing);

flag='l';  % default byte swap
fd=fopen(filename,'r',flag);
skip=fread(fd,[1 header_bytes/SHORT_SIZE],'short');  % header
%skip=fread(fd,[2 dim1],'int');                     
%skip=fread(fd,[2 dim1], pfile_type);                 % first baseline

raw=zeros(dim1,dim2,dim3);
for z=1:dim3
    for y=1:dim2
        
       %dump=fread(fd,[2 dim1],'int');      
       dump=fread(fd,[2 dim1], pfile_type);
       raw(:,y,z)=squeeze(dump(2,1:end)+1i*dump(1,:));
%         raw(:,y,z)=squeeze([0,dump(2,1:end-1)]+1i*dump(1,:));
       if ( mod(((z-1)*dim2+y),baseline_spacing)==0)
            %skip=fread(fd,[2 dim1],'int');
            %skip=fread(fd,[2 dim1], pfile_type);  % skip other baselines
       end
    end
end
fclose(fd);
end

%%%%% function [hdrp, header_bytes]=WrapGEheader(Pfilename)
%%%%% Wrap the GE header reader, adding pre-check for existance of Pfile
%%%%% slg
function [hdrp, header_bytes]=WrapGEheader(Pfilename)

[fd, msg]=fopen(Pfilename,'r'); 
holdfd = fd;
fclose(fd);
if holdfd < 0 
    mymsg = sprintf ('\nHeader reader unable to open Pfile %s: %s', Pfilename, msg)
    error(mymsg);
end
[hdrp, header_bytes] = read_gehdr_local(Pfilename);  % now we know Pfile is there
end

%%
%%%%% function [img] = CalcPhase(img,options)
%%%%% compute the phase map at 7T
%%%%% The phase map is computed by removing the coil phase. The coil phase
%%%%% is estimated using a low resolution image. Assume original images are stored
%%%%% in the same directory. Odd-numbered images are magnitude image;
%%%%% even-numbered images are phase image (x1000). The last two
%%%%% images for each slice are the sum-of-squares images.
%%%%% options.dim "3D", "2D"
%%%%% options.ncoil
%%%%% options.fermir: number of voxel per 512, e,g. 8/512
%%%%% options.fermiw: number of voxel per 512, e.g. 8/512
%%%%% Chunlei Liu, Stanford University, 08/18/2007,
%%%%%              Duke University, 2009
function [img_out] = CalcPhase(img,options)

if (nargin == 1)
    ncoil = 1;
end

if ( nargin == 2 && strcmp(options.dim,'3D'))
%     fprintf('3D...\n');
    opxres = size(img,1);
    opyres = size(img,2);
    opzres = size(img,3);
    ncoil = size(img,4);
    
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

end


function im = fft2c(d)
% im = fft2c(d)
%
% fft2c perfomrs a centered fft2
%
im = fftshift(fft2(fftshift(d)));
end

function im = ifft2c(d)
% im = ifft2c(d)
%
% ifft2c perfomrs a centered ifft2
%
im = fftshift(ifft2(fftshift(d)));
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

%% ++++++++++++++++++++++++++++++++++++++++++++
% sub routines to get header info
% +++++++++++++++++++++++++++++++++++++++++++++
%%
function [h, header_size] = read_gehdr_local(pfile)
%READ_GEHDR Read the header data in a GE LX raw data P-file
%
%  [H, BYTES] = READ_GEHDR(PFILE) reads header data from a GE LX raw
%  data P-file. PFILE can be a string containing the P-file name or the
%  integer run number for the series containing the P-file. The
%  returned structure, H, is similar to that found in GE include files
%  rdbm.h and imagedb.h and contains the following fields: H.RDB,
%  H.EXAM, H.SERIES, H.IMAGE, and H.DAQTAB.
%
%  The returned value, BYTES, is the number of bytes in the the header,
%  i.e. total header length.  this value can be used to open the data
%  file and skip past the header with FSEEK to go directly to the start
%  of the data section.
%
%  If PFILE does not correspond to a "proper" GE LX P-file, H will be
%  empty and a warning message will be displayed.  A "proper" P-file
%  must have a recognized version stored in the first 4 bytes and it
%  must have a minimum size of BYTES.
%
%  This function reads GE LX header versions 7.0 through 11.0.  Note
%  that these versions DO NOT necessarily correspond to the system
%  software version.
%
%  Most of the code in this function was generated automatically by the
%  Perl script gehdr2mat.pl (by DB Clayton) which parses the GE include
%  files, rdbm.h and imagedb.h.  See the header in gehdr2mat.pl for
%  more details on how to upgrade this script when new header versions
%  are released.
%
%  DB Clayton - v11.2 2005/10/7
%
%  Jian Zhang - 12/17/2008 - Add support for 20.0
%                            Missing prescan (starting 14.3)
%                            swiftcoilinfo may be wrong?
%
%  JZ         - 02/19/2009 - Full support for 14.3 & 20.0


% initialize outputs
h = [];  % header struct
header_size = 0;  % total header size
msgid = 'GE:readHeaderFailed';  % warning message identifier

% change numeric arg into pfile name
if (isnumeric(pfile)),
    pfile = sprintf('P%05d.7', pfile);
end

% get file version
id = fopen(pfile, 'r', 'l');  % try little-endian
fseek(id, 0, 'eof');
s = ftell(id);  % file size
if (s < 4)
    warning(msgid, 'file too small for GE raw data: file=%s, size=%d', pfile, s);
    fclose(id);
    return;
end
frewind(id);
ver = fread(id, 1, 'float');
if (~isfinite(ver) || ver < 7.0 || ver > 21.0)
    fclose(id);
    id = fopen(pfile, 'r', 'b');  % try big-endian
    ver = fread(id, 1, 'float');
    if (~isfinite(ver) || ver < 7.0 || ver > 21.0)
        warning(msgid, 'unrecognized GE header version: file=%s', pfile);
        fclose(id);
        return;
    end
end

% header parameters
prescan_header_length = 1500;
if (ver >= 7.0 && ver < 8.0)
    daqtab_size = 512;
    daqtab_offset = 10240;
    exam_offset = 36872;
    header_size = 39984;
elseif (ver >= 8.0 && ver < 9.0)
    daqtab_size = 1024;
    daqtab_offset = 10240;
    exam_offset = 57352;
    header_size = 60464;
elseif (ver >= 9.0 && ver < 11.0)
    daqtab_size = 1024;
    daqtab_offset = 10240;
    exam_offset = 57352;
    header_size = 61464;
elseif (ver == 11.0)
    daqtab_size = 1024;
    daqtab_offset = 10240;
    exam_offset = 61448;
    header_size = 66072;
elseif (ver > 12 && ver <= 14)     % 14X
    daqtab_size = 2048;
    daqtab_offset = 34816;
    exam_offset = 139272;
    header_size = 144408;
elseif (ver >= 14.2 && ver <= 14.5)    % 14.3
    daqtab_size = 2048;
    daqtab_offset = 34816;
    exam_offset = 140772;
    header_size = 145908;
elseif (ver >= 15.0 && ver <= 15.1)    % 15.0
    daqtab_size = 2048;
    daqtab_offset = 34816;
    exam_offset = 140772;
    header_size = 145908;    
elseif (ver >= 20)                  % 20X
    daqtab_size = 2048;
    daqtab_offset = 36864;
    exam_offset = 142820;
    header_size = 149788;
end

% check for minimum file size
fseek(id, 0, 'eof');
s = ftell(id);
if (s < header_size)
    warning(msgid, 'file too small for GE raw data: file=%s, size=%d', pfile, s);
    fclose(id);
    return;
end

% read rdb section
frewind(id);
if ver >= 20
    h.rdb = read_rdb_20(id,ver);
else
    h.rdb = read_rdb(id, ver);
end

% read daqtab section
fseek(id, daqtab_offset, 'bof');
h.daqtab = read_daqtab(id, ver, daqtab_size);

% read prescan section if ver>=14.3
if ver >= 14.3
    fseek(id,exam_offset-prescan_header_length,'bof');
    h.prescan = read_prescan(id);
end
    
% read exam, series, and image sections
fseek(id, exam_offset, 'bof');
h.exam = read_exam(id, ver);
h.series = read_series(id, ver);
h.image = read_image(id, ver);

% done
fclose(id);

end

%
%    subroutines
%


function rdb = read_rdb(id, ver)
rdb.rdbm_rev = fread(id, 1, 'float');  %
rdb.run_int = fread(id, 1, 'long');  % Rdy pkt Run Number
rdb.scan_seq = fread(id, 1, 'short');  % Rdy pkt Sequence Number
rdb.run_char = modchar(fread(id, 6, '*char'));  % Rdy pkt Run no in char
rdb.scan_date = modchar(fread(id, 10, '*char'));  %
rdb.scan_time = modchar(fread(id, 8, '*char'));  %
rdb.logo = modchar(fread(id, 10, '*char'));  % rdbm used to verify file
rdb.file_contents = fread(id, 1, 'short');  % Data type 0=emp 1=nrec 2=rw 0, 1, 2
rdb.lock_mode = fread(id, 1, 'short');  % unused
rdb.dacq_ctrl = fread(id, 1, 'short');  % rhdacqctrl bit mask 15 bits
rdb.recon_ctrl = fread(id, 1, 'short');  % rhrcctrl bit mask 15 bits
rdb.exec_ctrl = fread(id, 1, 'short');  % rhexecctrl bit mask 15 bits
rdb.scan_type = fread(id, 1, 'short');  % bit mask 15 bits
rdb.data_collect_type = fread(id, 1, 'short');  % rhtype bit mask 15 bits
rdb.data_format = fread(id, 1, 'short');  % rhformat bit mask 15 bits
rdb.recon = fread(id, 1, 'short');  % rhrecon proc-a-son recon 0 - 100
rdb.datacq = fread(id, 1, 'short');  % rhdatacq proc-a-son dacq
rdb.npasses = fread(id, 1, 'short');  % rhnpasses passes for a scan 0 - 256
rdb.npomp = fread(id, 1, 'short');  % rhnpomp pomp group slices 1,2
rdb.nslices = fread(id, 1, 'short');  % rhnslices slices in a pass 0 - 256
rdb.nechoes = fread(id, 1, 'short');  % rhnecho echoes of a slice 1 - 32
rdb.navs = fread(id, 1, 'short');  % rhnavs num of excitiations 1 - 32727
rdb.nframes = fread(id, 1, 'short');  % rhnframes yres 0 - 1024
rdb.baseline_views = fread(id, 1, 'short');  % rhbline baselines 0 - 1028
rdb.hnover = fread(id, 1, 'short');  % rhhnover overscans 0 - 1024
rdb.frame_size = fread(id, 1, 'ushort');  % rhfrsize xres 0 - 32768
rdb.point_size = fread(id, 1, 'short');  % rhptsize 2 - 4
rdb.vquant = fread(id, 1, 'short');  % rhvquant 3d volumes 1
rdb.cheart = fread(id, 1, 'short');  % RX Cine heart phases 1 - 32
rdb.ctr = fread(id, 1, 'float');  % RX Cine TR in sec 0 - 3.40282e38
rdb.ctrr = fread(id, 1, 'float');  % RX Cine RR in sec 0 - 30.0
rdb.initpass = fread(id, 1, 'short');  % rhinitpass allocate passes 0 - 32767
rdb.incrpass = fread(id, 1, 'short');  % rhincrpass tps autopauses 0 - 32767
rdb.method_ctrl = fread(id, 1, 'short');  % rhmethod 0=recon, 1=psd 0, 1
rdb.da_xres = fread(id, 1, 'ushort');  % rhdaxres 0 - 32768
rdb.da_yres = fread(id, 1, 'short');  % rhdayres 0 - 2049
rdb.rc_xres = fread(id, 1, 'short');  % rhrcxres 0 - 1024
rdb.rc_yres = fread(id, 1, 'short');  % rhrcyres 0 - 1024
rdb.im_size = fread(id, 1, 'short');  % rhimsize 0 - 512
rdb.rc_zres = fread(id, 1, 'long');  % power of 2 > rhnslices 0 - 128
rdb.raw_pass_size = fread(id, 1, 'ulong');  % rhrawsize 0 - 2147483647
rdb.sspsave = fread(id, 1, 'ulong');  % rhsspsave 0 - 2147483647
rdb.udasave = fread(id, 1, 'ulong');  % rhudasave 0 - 2147483647
rdb.fermi_radius = fread(id, 1, 'float');  % rhfermr fermi radius 0 - 3.40282e38
rdb.fermi_width = fread(id, 1, 'float');  % rhfermw fermi width 0 - 3.40282e38
rdb.fermi_ecc = fread(id, 1, 'float');  % rhferme fermi excentiricty 0 - 3.40282e38
rdb.clip_min = fread(id, 1, 'float');  % rhclipmin 4x IP limit +-16383
rdb.clip_max = fread(id, 1, 'float');  % rhclipmax 4x IP limit +-16383
rdb.default_offset = fread(id, 1, 'float');  % rhdoffset default offset = 0 +-3.40282e38
rdb.xoff = fread(id, 1, 'float');  % rhxoff scroll img in x +-256
rdb.yoff = fread(id, 1, 'float');  % rhyoff scroll img in y +-256
rdb.nwin = fread(id, 1, 'float');  % rhnwin hecho window width 0 - 256
rdb.ntran = fread(id, 1, 'float');  % rhntran hecho trans width 0 - 256
rdb.scalei = fread(id, 1, 'float');  % PS rhscalei +-3.40282e38
rdb.scaleq = fread(id, 1, 'float');  % PS rhscaleq def = 0 +-3.40282e38
rdb.rotation = fread(id, 1, 'short');  % RX 0 90 180 270 deg 0 - 3
rdb.transpose = fread(id, 1, 'short');  % RX 0, 1 n / y transpose 0 - 1
rdb.kissoff_views = fread(id, 1, 'short');  % rhblank zero image views 0 - 512
rdb.slblank = fread(id, 1, 'short');  % rhslblank slice blank 3d 0 - 128
rdb.gradcoil = fread(id, 1, 'short');  % RX 0=off 1=Schnk 2=Rmr 0 - 2
rdb.ddaover = fread(id, 1, 'short');  % rhddaover unused
rdb.sarr = fread(id, 1, 'short');  % SARR bit mask 15 bits
rdb.fd_tr = fread(id, 1, 'short');  % SARR feeder timing info
rdb.fd_te = fread(id, 1, 'short');  % SARR feeder timing info
rdb.fd_ctrl = fread(id, 1, 'short');  % SARR control of feeder
rdb.algor_num = fread(id, 1, 'short');  % SARR df decimation ratio
rdb.fd_df_dec = fread(id, 1, 'short');  % SARR which feeder algor
rdb = read_multi_rcv(id, rdb); % kluge for RDB_MULTI_RCV_TYPE
rdb.user0 = fread(id, 1, 'float');  % rhuser0 +-3.40282e38
rdb.user1 = fread(id, 1, 'float');  % rhuser1 +-3.40282e38
rdb.user2 = fread(id, 1, 'float');  % rhuser2 +-3.40282e38
rdb.user3 = fread(id, 1, 'float');  % rhuser3 +-3.40282e38
rdb.user4 = fread(id, 1, 'float');  % rhuser4 +-3.40282e38
rdb.user5 = fread(id, 1, 'float');  % rhuser5 +-3.40282e38
rdb.user6 = fread(id, 1, 'float');  % rhuser6 +-3.40282e38
rdb.user7 = fread(id, 1, 'float');  % rhuser7 +-3.40282e38
rdb.user8 = fread(id, 1, 'float');  % rhuser8 +-3.40282e38
rdb.user9 = fread(id, 1, 'float');  % rhuser9 +-3.40282e38
rdb.user10 = fread(id, 1, 'float');  % rhuser10 +-3.40282e38
rdb.user11 = fread(id, 1, 'float');  % rhuser11 +-3.40282e38
rdb.user12 = fread(id, 1, 'float');  % rhuser12 +-3.40282e38
rdb.user13 = fread(id, 1, 'float');  % rhuser13 +-3.40282e38
rdb.user14 = fread(id, 1, 'float');  % rhuser14 +-3.40282e38
rdb.user15 = fread(id, 1, 'float');  % rhuser15 +-3.40282e38
rdb.user16 = fread(id, 1, 'float');  % rhuser16 +-3.40282e38
rdb.user17 = fread(id, 1, 'float');  % rhuser17 +-3.40282e38
rdb.user18 = fread(id, 1, 'float');  % rhuser18 +-3.40282e38
rdb.user19 = fread(id, 1, 'float');  % rhuser19 +-3.40282e38
rdb.v_type = fread(id, 1, 'long');  % rhvtype bit mask 31 bits
rdb.v_coefxa = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefxb = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefxc = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefxd = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefya = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefyb = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefyc = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefyd = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefza = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.v_coefzb = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.v_coefzc = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.v_coefzd = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.vm_coef1 = fread(id, 1, 'float');  % RX weight for mag image 1 0 - 1
rdb.vm_coef2 = fread(id, 1, 'float');  % RX weight for mag image 2 0 - 1
rdb.vm_coef3 = fread(id, 1, 'float');  % RX weight for mag image 3 0 - 1
rdb.vm_coef4 = fread(id, 1, 'float');  % RX weight for mag image 4 0 - 1
rdb.v_venc = fread(id, 1, 'float');  % RX vel encodeing cm / sec 0.001 - 5000
rdb.spectral_width = fread(id, 1, 'float');  % specwidth filter width kHz 500 - 3355432
rdb.csi_dims = fread(id, 1, 'short');  % spectro
rdb.xcsi = fread(id, 1, 'short');  % rhspecrescsix 2 - 64
rdb.ycsi = fread(id, 1, 'short');  % rhspecrescsiy 2 - 64
rdb.zcsi = fread(id, 1, 'short');  % spectro
rdb.roilenx = fread(id, 1, 'float');  % RX x csi volume dimension
rdb.roileny = fread(id, 1, 'float');  % RX y csi volume dimension
rdb.roilenz = fread(id, 1, 'float');  % RX z csi volume dimension
rdb.roilocx = fread(id, 1, 'float');  % RX x csi volume center
rdb.roilocy = fread(id, 1, 'float');  % RX y csi volume center
rdb.roilocz = fread(id, 1, 'float');  % RX z csi volume center
rdb.numdwell = fread(id, 1, 'float');  % specdwells 0 - 3.40282e38
rdb.ps_command = fread(id, 1, 'long');  % PS internal use only
rdb.ps_mps_r1 = fread(id, 1, 'long');  % PS MPS R1 setting 1 - 7
rdb.ps_mps_r2 = fread(id, 1, 'long');  % PS MPS R2 setting 1 - 30
rdb.ps_mps_tg = fread(id, 1, 'long');  % PS MPS Transmit gain setting 0 - 200
rdb.ps_mps_freq = fread(id, 1, 'long');  % PS MPS Center frequency hz +-3.40282e38
rdb.ps_aps_r1 = fread(id, 1, 'long');  % PS APS R1 setting 1 - 7
rdb.ps_aps_r2 = fread(id, 1, 'long');  % PS APS R2 setting 1 - 30
rdb.ps_aps_tg = fread(id, 1, 'long');  % PS APS Transmit gain setting 0 - 200
rdb.ps_aps_freq = fread(id, 1, 'long');  % PS APS Center frequency hz +-3.40282e38
rdb.ps_scalei = fread(id, 1, 'float');  % PS rational scaling +-3.40282e38
rdb.ps_scaleq = fread(id, 1, 'float');  % PS unused
rdb.ps_snr_warning = fread(id, 1, 'long');  % PS noise test 0=16 1=32 bits 0, 1
rdb.ps_aps_or_mps = fread(id, 1, 'long');  % PS prescan order logic 0 - 5
rdb.ps_mps_bitmap = fread(id, 1, 'long');  % PS bit mask 4 bits
rdb.ps_powerspec = modchar(fread(id, 256, '*char'));  % PS
rdb.ps_filler1 = fread(id, 1, 'long');  % PS filler
rdb.ps_filler2 = fread(id, 1, 'long');  % PS filler
rdb.rec_noise_mean = fread(id, 16, 'float');  % PS mean noise each receiver +-3.40282e38
rdb.rec_noise_std = fread(id, 16, 'float');  % PS noise calc for muti rec +-3.40282e38
rdb.halfecho = fread(id, 1, 'short');  % spectro full, half echo 0, 1
rdb.im_size_y = fread(id, 1, 'short');  % rh???? 0 - 512
rdb.data_collect_type1 = fread(id, 1, 'long');  % rh???? bit mask 31 bits
rdb.freq_scale = fread(id, 1, 'float');  % rh???? freq k-space step +-3.40282e38
rdb.phase_scale = fread(id, 1, 'float');  % rh???? freq k-space step +-3.40282e38
rdb.ovl = fread(id, 1, 'short');  % rhovl - overlaps for MOTSA
rdb.pclin = fread(id, 1, 'short');  % Linear Corr. 0:off, 1:linear, 2:polynomial
rdb.pclinnpts = fread(id, 1, 'short');  % fit number of points
rdb.pclinorder = fread(id, 1, 'short');  % fit order
rdb.pclinavg = fread(id, 1, 'short');  % linear phase corr avg 0:off, 1:on
rdb.pccon = fread(id, 1, 'short');  % Const Corr. 0:off, 1:Ky spec., 2:polyfit(2/ilv), 3:polyfit(1/ilv)
rdb.pcconnpts = fread(id, 1, 'short');  % fit number of points
rdb.pcconorder = fread(id, 1, 'short');  % fit order
rdb.pcextcorr = fread(id, 1, 'short');  % external correction file 0:don't use, 1: use
rdb.pcgraph = fread(id, 1, 'short');  % Phase Correction coef. image 0:off, 1:linear & constant
rdb.pcileave = fread(id, 1, 'short');  % Interleaves to use for correction: 0=all, 1=only first
rdb.hdbestky = fread(id, 1, 'short');  % bestky view for fractional Ky scan
rdb.pcctrl = fread(id, 1, 'short');  % phase correction research control
rdb.pcthrespts = fread(id, 1, 'short');  % 2..512 adjacent points
rdb.pcdiscbeg = fread(id, 1, 'short');  % 0..512 beginning point to discard
rdb.pcdiscmid = fread(id, 1, 'short');  % 0..512 middle point to discard
rdb.pcdiscend = fread(id, 1, 'short');  % 0..512 ending point to discard
rdb.pcthrespct = fread(id, 1, 'short');  % Threshold percentage
rdb.pcspacial = fread(id, 1, 'short');  % Spacial best ref scan index 0..512
rdb.pctemporal = fread(id, 1, 'short');  % Temporal best ref scan index 0..512
rdb.pcspare = fread(id, 1, 'short');  % spare for phase correction
rdb.ileaves = fread(id, 1, 'short');  % Number of interleaves
rdb.kydir = fread(id, 1, 'short');  % Ky traversal dircetion 0: top-down, 1:center out
rdb.alt = fread(id, 1, 'short');  % Alt read sign 0=no, 1=odd/even, 2=pairs
rdb.reps = fread(id, 1, 'short');  % Number of scan repetitions
rdb.ref = fread(id, 1, 'short');  % Ref Scan 0: off 1: on
rdb.pcconnorm = fread(id, 1, 'float');  % Constant S term normalization factor
rdb.pcconfitwt = fread(id, 1, 'float');  % Constant polyfit weighting factor
rdb.pclinnorm = fread(id, 1, 'float');  % Linear S term normalization factor
rdb.pclinfitwt = fread(id, 1, 'float');  % Linear polyfit weighting factor
rdb.pcbestky = fread(id, 1, 'float');  % Best Ky location
rdb.vrgf = fread(id, 1, 'long');  % control word for VRG filter
rdb.vrgfxres = fread(id, 1, 'long');  % control word for VRGF final x resolution
rdb.bp_corr = fread(id, 1, 'long');  % control word for bandpass asymmetry
rdb.recv_freq_s = fread(id, 1, 'float');  % starting frequency (+62.5)
rdb.recv_freq_e = fread(id, 1, 'float');  % ending frequency (-62.5)
rdb.hniter = fread(id, 1, 'long');  % Selects the number of (continued...)
rdb.fast_rec = fread(id, 1, 'long');  % Added for homodyne II, tells if (continued...)
rdb.refframes = fread(id, 1, 'long');  % total # of frames for ref scan
rdb.refframep = fread(id, 1, 'long');  % # of frames per pass for a ref scan
rdb.scnframe = fread(id, 1, 'long');  % total # of frames for a entire scan
rdb.pasframe = fread(id, 1, 'long');  % # of frames per pass
rdb.user_usage_tag = fread(id, 1, 'ulong');  % for spectro
rdb.user_fill_mapMSW = fread(id, 1, 'ulong');  % for spectro
rdb.user_fill_mapLSW = fread(id, 1, 'ulong');  % for Spectro
rdb.user20 = fread(id, 1, 'float');  % all following usercv are for spectro
rdb.user21 = fread(id, 1, 'float');  %
rdb.user22 = fread(id, 1, 'float');  %
rdb.user23 = fread(id, 1, 'float');  %
rdb.user24 = fread(id, 1, 'float');  %
rdb.user25 = fread(id, 1, 'float');  %
rdb.user26 = fread(id, 1, 'float');  %
rdb.user27 = fread(id, 1, 'float');  %
rdb.user28 = fread(id, 1, 'float');  %
rdb.user29 = fread(id, 1, 'float');  %
rdb.user30 = fread(id, 1, 'float');  %
rdb.user31 = fread(id, 1, 'float');  %
rdb.user32 = fread(id, 1, 'float');  %
rdb.user33 = fread(id, 1, 'float');  %
rdb.user34 = fread(id, 1, 'float');  %
rdb.user35 = fread(id, 1, 'float');  %
rdb.user36 = fread(id, 1, 'float');  %
rdb.user37 = fread(id, 1, 'float');  %
rdb.user38 = fread(id, 1, 'float');  %
rdb.user39 = fread(id, 1, 'float');  %
rdb.user40 = fread(id, 1, 'float');  %
rdb.user41 = fread(id, 1, 'float');  %
rdb.user42 = fread(id, 1, 'float');  %
rdb.user43 = fread(id, 1, 'float');  %
rdb.user44 = fread(id, 1, 'float');  %
rdb.user45 = fread(id, 1, 'float');  %
rdb.user46 = fread(id, 1, 'float');  %
rdb.user47 = fread(id, 1, 'float');  %
rdb.user48 = fread(id, 1, 'float');  %
rdb.pcfitorig = fread(id, 1, 'short');  % Adjust view indexes if set so bestky view = 0
rdb.pcshotfirst = fread(id, 1, 'short');  % First view within an echo group used for fit
rdb.pcshotlast = fread(id, 1, 'short');  % Last view within an echo group used for fit
rdb.pcmultegrp = fread(id, 1, 'short');  % If = 1, force pts from other egrps to be used
rdb.pclinfix = fread(id, 1, 'short');  % If = 2, force slope to be set to pclinslope
rdb.pcconfix = fread(id, 1, 'short');  % If = 2, force slope to be set to pcconslope
rdb.pclinslope = fread(id, 1, 'float');  % Value to set lin slope to if forced
rdb.pcconslope = fread(id, 1, 'float');  % Value to set con slope to if forced
rdb.pccoil = fread(id, 1, 'short');  % If 1,2,3,4, use that coil's results for all
rdb.vvsmode = fread(id, 1, 'short');  % Variable view sharing mode
rdb.vvsaimgs = fread(id, 1, 'short');  % number of original images
rdb.vvstr = fread(id, 1, 'short');  % TR in microseconds
rdb.vvsgender = fread(id, 1, 'short');  % gender: male or female
rdb.zip_factor = fread(id, 1, 'short');  % Slice ZIP factor: 0=OFF, 2, or 4
rdb.maxcoef1a = fread(id, 1, 'float');  % Coefficient A for flow image 1
rdb.maxcoef1b = fread(id, 1, 'float');  % Coefficient B for flow image 1
rdb.maxcoef1c = fread(id, 1, 'float');  % Coefficient C for flow image 1
rdb.maxcoef1d = fread(id, 1, 'float');  % Coefficient D for flow image 1
rdb.maxcoef2a = fread(id, 1, 'float');  % Coefficient A for flow image 2
rdb.maxcoef2b = fread(id, 1, 'float');  % Coefficient B for flow image 2
rdb.maxcoef2c = fread(id, 1, 'float');  % Coefficient C for flow image 2
rdb.maxcoef2d = fread(id, 1, 'float');  % Coefficient D for flow image 2
rdb.maxcoef3a = fread(id, 1, 'float');  % Coefficient A for flow image 3
rdb.maxcoef3b = fread(id, 1, 'float');  % Coefficient B for flow image 3
rdb.maxcoef3c = fread(id, 1, 'float');  % Coefficient C for flow image 3
rdb.maxcoef3d = fread(id, 1, 'float');  % Coefficient D for flow image 3
rdb.ut_ctrl = fread(id, 1, 'long');  % System utility control variable
rdb.dp_type = fread(id, 1, 'short');  % EPI II diffusion control cv
if (ver <= 7.0)  % end of rdbm version 7.0
    rdb.excess = fread(id, 423, 'short');  % free space for later expansion
    return
end
rdb.arw = fread(id, 1, 'short');  % Arrhythmia rejection window(percentage:1-100)
rdb.vps = fread(id, 1, 'short');  % View Per Segment for FastCine
rdb.mcReconEnable = fread(id, 1, 'short');  % N-Coil recon map
rdb.fov = fread(id, 1, 'float');  % Auto-NCoil
rdb.te = fread(id, 1, 'long');  % TE for first echo
rdb.te2 = fread(id, 1, 'long');  % TE for second and later echoes
rdb.dfmrbw = fread(id, 1, 'float');  % BW for navigator frames
rdb.dfmctrl = fread(id, 1, 'long');  % Control flag for dfm (0=off, other=on)
rdb.raw_nex = fread(id, 1, 'long');  % Uncombined NEX at start of recon
rdb.navs_per_pass = fread(id, 1, 'long');  % Max. navigator frames in a pass
rdb.dfmxres = fread(id, 1, 'long');  % xres of navigator frames
rdb.dfmptsize = fread(id, 1, 'long');  % point size of navigator frames
rdb.navs_per_view = fread(id, 1, 'long');  % Num. navigators per frame (tag table)
rdb.dfmdebug = fread(id, 1, 'long');  % control flag for dfm debug
rdb.dfmthreshold = fread(id, 1, 'float');  % threshold for navigator correction
rdb.grid_control = fread(id, 1, 'short');  % bit settings controlling gridding
rdb.b0map = fread(id, 1, 'short');  % B0 map enable and map size
rdb.grid_tediff = fread(id, 1, 'short');  % TE difference between b0 map arms
rdb.grid_motion_comp = fread(id, 1, 'short');  % flag to apply motion compensation
rdb.grid_radius_a = fread(id, 1, 'float');  % variable density transition
rdb.grid_radius_b = fread(id, 1, 'float');  % variable density transition
rdb.grid_max_gradient = fread(id, 1, 'float');  % Max gradient amplitude
rdb.grid_max_slew = fread(id, 1, 'float');  % Max slew rate
rdb.grid_scan_fov = fread(id, 1, 'float');  % Rx scan field of view
rdb.grid_a2d_time = fread(id, 1, 'float');  % A to D sample time microsecs
rdb.grid_density_factor = fread(id, 1, 'float');  % change factor for variable density
rdb.grid_display_fov = fread(id, 1, 'float');  % Rx display field of view
if (ver > 8.0)
    rdb.fatwater = fread(id, 1, 'short');  % for Fat and Water Dual Recon
    rdb.fiestamlf = fread(id, 1, 'short');  % MFO FIESTA recon control bit 16bits
end
rdb.app = fread(id, 1, 'short');  % Auto Post-Processing opcode
rdb.rhncoilsel = fread(id, 1, 'short');  % Auto-Ncoil
rdb.rhncoillimit = fread(id, 1, 'short');  % Auto-Ncoil
rdb.app_option = fread(id, 1, 'short');  % Auto Post_processing options
rdb.grad_mode = fread(id, 1, 'short');  % Gradient mode in Gemini project
if (ver > 8.0)
    rdb.pfile_passes = fread(id, 1, 'short');  % Num passes stored in a multi-pass Pfile (0 means 1 pass)
else
    rdb.edb_hdr_excess2 = fread(id, 1, 'short');  % To make sure that next member is 4-byte aligned
end
rdb.asset = fread(id, 1, 'int');  %
rdb.asset_calthresh = fread(id, 1, 'int');  %
rdb.asset_R = fread(id, 1, 'float');  %
rdb.coilno = fread(id, 1, 'int');  %
rdb.asset_phases = fread(id, 1, 'int');  %
rdb.scancent = fread(id, 1, 'float');  % Table position
rdb.position = fread(id, 1, 'int');  % Patient position
rdb.entry = fread(id, 1, 'int');  % Patient entry
rdb.lmhor = fread(id, 1, 'float');  % Landmark
if (ver <= 8.0)  % end of rdbm version 8.0
    rdb.excess = fread(id, 352, 'short');  % free space for later expansion
    return
end
rdb.last_slice_num = fread(id, 1, 'int');  %
rdb.asset_slice_R = fread(id, 1, 'float');  % Slice reduction factor
rdb.asset_slabwrap = fread(id, 1, 'float');  %
rdb.dwnav_coeff = fread(id, 1, 'float');  % Coeff for amount of phase correction
rdb.dwnav_cor = fread(id, 1, 'short');  % Navigator echo correction
rdb.dwnav_view = fread(id, 1, 'short');  % Num of views of nav echoes
rdb.dwnav_corecho = fread(id, 1, 'short');  % Num of nav echoes for actual correction
rdb.dwnav_sview = fread(id, 1, 'short');  % Start view for phase correction process
rdb.dwnav_eview = fread(id, 1, 'short');  % End view for phase correction process
rdb.dwnav_sshot = fread(id, 1, 'short');  % Start shot for delta phase estimation in nav echoes
rdb.dwnav_eshot = fread(id, 1, 'short');  % End shot for delta phase estimation in nav echoes
rdb.win3d_type = fread(id, 1, 'short');  % 0 = Modified Hanning, 1 = modified Tukey
rdb.win3d_apod = fread(id, 1, 'float');  % degree of apodization; 0.0 = boxcar, 1.0=hanning
rdb.win3d_q = fread(id, 1, 'float');  % apodization at ends, 0.0 = max, 1.0 = boxcar
rdb.ime_scic_enable = fread(id, 1, 'short');  % Surface Coil Intensity Correction: 1 if enabled
rdb.clariview_type = fread(id, 1, 'short');  % Type of Clariview/Name of Filter
rdb.ime_scic_edge = fread(id, 1, 'float');  % Edge paramaters for Enhanced Recon
rdb.ime_scic_smooth = fread(id, 1, 'float');  % Smooth paramaters for Enhanced Recon
rdb.ime_scic_focus = fread(id, 1, 'float');  % Focus paramaters for Enhanced Recon
rdb.clariview_edge = fread(id, 1, 'float');  % Edge paramaters for clariview
rdb.clariview_smooth = fread(id, 1, 'float');  % Smooth paramaters for clariview
rdb.clariview_focus = fread(id, 1, 'float');  % Focus paramaters for clariview
rdb.scic_reduction = fread(id, 1, 'float');  % Reduction paramater for SCIC
rdb.scic_gauss = fread(id, 1, 'float');  % Gauss paramater for SCIC
rdb.scic_threshold = fread(id, 1, 'float');  % Threshold paramater for SCIC
rdb.ectricks_no_regions = fread(id, 1, 'long');  % Total no of regions acquired by PSD
rdb.ectricks_input_regions = fread(id, 1, 'long');  % Total no of input regions for reordering
rdb.psc_reuse = fread(id, 1, 'short');  % Header field for smart prescan
rdb.left_blank = fread(id, 1, 'short');  %
rdb.right_blank = fread(id, 1, 'short');  %
rdb.acquire_type = fread(id, 1, 'short');  % Acquire type information from CV
rdb.retro_control = fread(id, 1, 'short');  % Retrosective FSE phase correction control flag. (continued...)
rdb.etl = fread(id, 1, 'short');  % Added for Retrospective FSE phase correction. This (continued...)
if (ver <= 9.0)  % end of rdbm version 9.0
    rdb.excess = fread(id, 300, 'short');  % free space for later expansion
    return
end
rdb.pcref_start = fread(id, 1, 'short');  % 1st view to use for dynamic EPI phase correction.
rdb.pcref_stop = fread(id, 1, 'short');  % Last view to use for dynamic EPI phase correction.
rdb.ref_skip = fread(id, 1, 'short');  % Number of passes to skip for dynamic EPI phase correction.
rdb.extra_frames_top = fread(id, 1, 'short');  % Number of extra frames at top of K-space
rdb.extra_frames_bot = fread(id, 1, 'short');  % Number of extra frames at bottom of K-space
rdb.multiphase_type = fread(id, 1, 'short');  % 0 = INTERLEAVED , 1 = SEQUENTIAL
rdb.nphases = fread(id, 1, 'short');  % Number of phases in a multiphase scan
rdb.pure = fread(id, 1, 'short');  % PURE flag from psd
rdb.pure_scale = fread(id, 1, 'float');  % Recon scale factor ratio for cal scan
rdb.off_data = fread(id, 1, 'int');  % Byte offset to start of raw data (i.e size of POOL_HEADER)
rdb.off_per_pass = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_per_pass of POOL_HEADER
rdb.off_unlock_raw = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_unlock_raw of POOL_HEADER
rdb.off_data_acq_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_data_acq_tab of POOL_HEADER
rdb.off_nex_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_nex_tab of POOL_HEADER
rdb.off_nex_abort_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_nex_abort_tab of POOL_HEADER
rdb.off_tool = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_tool of POOL_HEADER
rdb.off_exam = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_exam of POOL_HEADER
rdb.off_series = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_series of POOL_HEADER
rdb.off_image = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_image of POOL_HEADER
rdb.off_spare_a = fread(id, 1, 'int');  % spare
rdb.off_spare_b = fread(id, 1, 'int');  % spare
rdb.new_wnd_level_flag = fread(id, 1, 'int');  % New WW/WL algo enable/disable flag
rdb.wnd_image_hist_area = fread(id, 1, 'int');  % Image Area %
rdb.wnd_high_hist = fread(id, 1, 'float');  % Histogram Area Top
rdb.wnd_lower_hist = fread(id, 1, 'float');  % Histogram Area Bottom
rdb.pure_filter = fread(id, 1, 'short');  % PURE noise reduction on=1/off=0
rdb.cfg_pure_filter = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_fit_order = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_kernelsize_z = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_kernelsize_xy = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_weight_radius = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_intensity_scale = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_noise_threshold = fread(id, 1, 'short');  % PURE cfg file value
rdb.wienera = fread(id, 1, 'float');  % NB maintain alignment of floats
rdb.wienerb = fread(id, 1, 'float');  %
rdb.wienert2 = fread(id, 1, 'float');  %
rdb.wieneresp = fread(id, 1, 'float');  %
rdb.wiener = fread(id, 1, 'short');  %
rdb.flipfilter = fread(id, 1, 'short');  %
rdb.dbgrecon = fread(id, 1, 'short');  %
rdb.ech2skip = fread(id, 1, 'short');  %
rdb.tricks_type = fread(id, 1, 'int');  % 0 = Subtracted, 1 = Unsubtracted
rdb.lcfiesta_phase = fread(id, 1, 'float');  % LC Fiesta
rdb.lcfiesta = fread(id, 1, 'short');  % LC Fiesta
rdb.herawflt = fread(id, 1, 'short');  % Half echo raw data filter
rdb.herawflt_befnwin = fread(id, 1, 'short');  % Half echo raw data filter
rdb.herawflt_befntran = fread(id, 1, 'short');  % Half echo raw data filter
rdb.herawflt_befamp = fread(id, 1, 'float');  % Half echo raw data filter
rdb.herawflt_hpfamp = fread(id, 1, 'float');  % Half echo raw data filter
rdb.heover = fread(id, 1, 'short');  % Half echo over sampling
rdb.pure_correction_threshold = fread(id, 1, 'short');  % PURE Correction threshold

if ver < 14
    rdb.hdr_ps_autoshim_status = fread(id, 1, 'int'); 
    rdb.excess = fread(id, 222, 'short');  % free space for later expansion
else
    rdb.swiftenable  = fread(id, 1, 'int');
    rdb.numslabs = fread(id, 1, 'short');
    rdb.swiftcoilnos = fread(id, 1, 'ushort');
    rdn.ps_autoshim_status = fread(id,1,'int');
    
    if ver < 14.3
        rdb.excess = fread(id, 218, 'short');  % free space for later expansion
    else
        rdb.dynaplan_numphases = fread(id, 1, 'int');
        rdb.excess = fread(id, 216, 'short');
    end
end

end

function rdb = read_rdb_20(id, ver)
rdb.rdbm_rev = fread(id, 1, 'float');  %
rdb.run_int = fread(id, 1, 'long');  % Rdy pkt Run Number
rdb.scan_seq = fread(id, 1, 'short');  % Rdy pkt Sequence Number
rdb.run_char = modchar(fread(id, 6, '*char'));  % Rdy pkt Run no in char
rdb.scan_date = modchar(fread(id, 10, '*char'));  %
rdb.scan_time = modchar(fread(id, 8, '*char'));  %
rdb.logo = modchar(fread(id, 10, '*char'));  % rdbm used to verify file
rdb.file_contents = fread(id, 1, 'short');  % Data type 0=emp 1=nrec 2=rw 0, 1, 2
rdb.lock_mode = fread(id, 1, 'short');  % unused
rdb.dacq_ctrl = fread(id, 1, 'short');  % rhdacqctrl bit mask 15 bits
rdb.recon_ctrl = fread(id, 1, 'short');  % rhrcctrl bit mask 15 bits
rdb.exec_ctrl = fread(id, 1, 'short');  % rhexecctrl bit mask 15 bits
rdb.scan_type = fread(id, 1, 'short');  % bit mask 15 bits
rdb.data_collect_type = fread(id, 1, 'short');  % rhtype bit mask 15 bits
rdb.data_format = fread(id, 1, 'short');  % rhformat bit mask 15 bits
rdb.recon = fread(id, 1, 'short');  % rhrecon proc-a-son recon 0 - 100
rdb.datacq = fread(id, 1, 'short');  % rhdatacq proc-a-son dacq
rdb.npasses = fread(id, 1, 'short');  % rhnpasses passes for a scan 0 - 256
rdb.npomp = fread(id, 1, 'short');  % rhnpomp pomp group slices 1,2
rdb.nslices = fread(id, 1, 'short');  % rhnslices slices in a pass 0 - 256
rdb.nechoes = fread(id, 1, 'short');  % rhnecho echoes of a slice 1 - 32
rdb.navs = fread(id, 1, 'short');  % rhnavs num of excitiations 1 - 32727
rdb.nframes = fread(id, 1, 'short');  % rhnframes yres 0 - 1024
rdb.baseline_views = fread(id, 1, 'short');  % rhbline baselines 0 - 1028
rdb.hnover = fread(id, 1, 'short');  % rhhnover overscans 0 - 1024
rdb.frame_size = fread(id, 1, 'ushort');  % rhfrsize xres 0 - 32768
rdb.point_size = fread(id, 1, 'short');  % rhptsize 2 - 4
rdb.vquant = fread(id, 1, 'short');  % rhvquant 3d volumes 1
rdb.cheart = fread(id, 1, 'short');  % RX Cine heart phases 1 - 32
rdb.ctr = fread(id, 1, 'float');  % RX Cine TR in sec 0 - 3.40282e38
rdb.ctrr = fread(id, 1, 'float');  % RX Cine RR in sec 0 - 30.0
rdb.initpass = fread(id, 1, 'short');  % rhinitpass allocate passes 0 - 32767
rdb.incrpass = fread(id, 1, 'short');  % rhincrpass tps autopauses 0 - 32767
rdb.method_ctrl = fread(id, 1, 'short');  % rhmethod 0=recon, 1=psd 0, 1
rdb.da_xres = fread(id, 1, 'ushort');  % rhdaxres 0 - 32768
rdb.da_yres = fread(id, 1, 'short');  % rhdayres 0 - 2049
rdb.rc_xres = fread(id, 1, 'short');  % rhrcxres 0 - 1024
rdb.rc_yres = fread(id, 1, 'short');  % rhrcyres 0 - 1024
rdb.im_size = fread(id, 1, 'short');  % rhimsize 0 - 512
rdb.rc_zres = fread(id, 1, 'long');  % power of 2 > rhnslices 0 - 128

rdb.raw_pass_size_deprecated = fread(id, 1, 'ulong');  % rhrawsize 0 - 2147483647
rdb.sspsave_deprecated = fread(id, 1, 'ulong');  % rhsspsave 0 - 2147483647
rdb.udasave_deprecated = fread(id, 1, 'ulong');  % rhudasave 0 - 2147483647

rdb.fermi_radius = fread(id, 1, 'float');  % rhfermr fermi radius 0 - 3.40282e38
rdb.fermi_width = fread(id, 1, 'float');  % rhfermw fermi width 0 - 3.40282e38
rdb.fermi_ecc = fread(id, 1, 'float');  % rhferme fermi excentiricty 0 - 3.40282e38
rdb.clip_min = fread(id, 1, 'float');  % rhclipmin 4x IP limit +-16383
rdb.clip_max = fread(id, 1, 'float');  % rhclipmax 4x IP limit +-16383
rdb.default_offset = fread(id, 1, 'float');  % rhdoffset default offset = 0 +-3.40282e38
rdb.xoff = fread(id, 1, 'float');  % rhxoff scroll img in x +-256
rdb.yoff = fread(id, 1, 'float');  % rhyoff scroll img in y +-256
rdb.nwin = fread(id, 1, 'float');  % rhnwin hecho window width 0 - 256
rdb.ntran = fread(id, 1, 'float');  % rhntran hecho trans width 0 - 256
rdb.scalei = fread(id, 1, 'float');  % PS rhscalei +-3.40282e38
rdb.scaleq = fread(id, 1, 'float');  % PS rhscaleq def = 0 +-3.40282e38
rdb.rotation = fread(id, 1, 'short');  % RX 0 90 180 270 deg 0 - 3
rdb.transpose = fread(id, 1, 'short');  % RX 0, 1 n / y transpose 0 - 1
rdb.kissoff_views = fread(id, 1, 'short');  % rhblank zero image views 0 - 512
rdb.slblank = fread(id, 1, 'short');  % rhslblank slice blank 3d 0 - 128
rdb.gradcoil = fread(id, 1, 'short');  % RX 0=off 1=Schnk 2=Rmr 0 - 2
rdb.ddaover = fread(id, 1, 'short');  % rhddaover unused
rdb.sarr = fread(id, 1, 'short');  % SARR bit mask 15 bits
rdb.fd_tr = fread(id, 1, 'short');  % SARR feeder timing info
rdb.fd_te = fread(id, 1, 'short');  % SARR feeder timing info
rdb.fd_ctrl = fread(id, 1, 'short');  % SARR control of feeder
rdb.algor_num = fread(id, 1, 'short');  % SARR df decimation ratio
rdb.fd_df_dec = fread(id, 1, 'short');  % SARR which feeder algor
rdb = read_multi_rcv(id, rdb); % kluge for RDB_MULTI_RCV_TYPE
rdb.user0 = fread(id, 1, 'float');  % rhuser0 +-3.40282e38
rdb.user1 = fread(id, 1, 'float');  % rhuser1 +-3.40282e38
rdb.user2 = fread(id, 1, 'float');  % rhuser2 +-3.40282e38
rdb.user3 = fread(id, 1, 'float');  % rhuser3 +-3.40282e38
rdb.user4 = fread(id, 1, 'float');  % rhuser4 +-3.40282e38
rdb.user5 = fread(id, 1, 'float');  % rhuser5 +-3.40282e38
rdb.user6 = fread(id, 1, 'float');  % rhuser6 +-3.40282e38
rdb.user7 = fread(id, 1, 'float');  % rhuser7 +-3.40282e38
rdb.user8 = fread(id, 1, 'float');  % rhuser8 +-3.40282e38
rdb.user9 = fread(id, 1, 'float');  % rhuser9 +-3.40282e38
rdb.user10 = fread(id, 1, 'float');  % rhuser10 +-3.40282e38
rdb.user11 = fread(id, 1, 'float');  % rhuser11 +-3.40282e38
rdb.user12 = fread(id, 1, 'float');  % rhuser12 +-3.40282e38
rdb.user13 = fread(id, 1, 'float');  % rhuser13 +-3.40282e38
rdb.user14 = fread(id, 1, 'float');  % rhuser14 +-3.40282e38
rdb.user15 = fread(id, 1, 'float');  % rhuser15 +-3.40282e38
rdb.user16 = fread(id, 1, 'float');  % rhuser16 +-3.40282e38
rdb.user17 = fread(id, 1, 'float');  % rhuser17 +-3.40282e38
rdb.user18 = fread(id, 1, 'float');  % rhuser18 +-3.40282e38
rdb.user19 = fread(id, 1, 'float');  % rhuser19 +-3.40282e38
rdb.v_type = fread(id, 1, 'long');  % rhvtype bit mask 31 bits
rdb.v_coefxa = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefxb = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefxc = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefxd = fread(id, 1, 'float');  % RX x flow direction control 0 - 4
rdb.v_coefya = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefyb = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefyc = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefyd = fread(id, 1, 'float');  % RX y flow direction control 0 - 4
rdb.v_coefza = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.v_coefzb = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.v_coefzc = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.v_coefzd = fread(id, 1, 'float');  % RX z flow direction control 0 - 4
rdb.vm_coef1 = fread(id, 1, 'float');  % RX weight for mag image 1 0 - 1
rdb.vm_coef2 = fread(id, 1, 'float');  % RX weight for mag image 2 0 - 1
rdb.vm_coef3 = fread(id, 1, 'float');  % RX weight for mag image 3 0 - 1
rdb.vm_coef4 = fread(id, 1, 'float');  % RX weight for mag image 4 0 - 1
rdb.v_venc = fread(id, 1, 'float');  % RX vel encodeing cm / sec 0.001 - 5000
rdb.spectral_width = fread(id, 1, 'float');  % specwidth filter width kHz 500 - 3355432
rdb.csi_dims = fread(id, 1, 'short');  % spectro
rdb.xcsi = fread(id, 1, 'short');  % rhspecrescsix 2 - 64
rdb.ycsi = fread(id, 1, 'short');  % rhspecrescsiy 2 - 64
rdb.zcsi = fread(id, 1, 'short');  % spectro
rdb.roilenx = fread(id, 1, 'float');  % RX x csi volume dimension
rdb.roileny = fread(id, 1, 'float');  % RX y csi volume dimension
rdb.roilenz = fread(id, 1, 'float');  % RX z csi volume dimension
rdb.roilocx = fread(id, 1, 'float');  % RX x csi volume center
rdb.roilocy = fread(id, 1, 'float');  % RX y csi volume center
rdb.roilocz = fread(id, 1, 'float');  % RX z csi volume center
rdb.numdwell = fread(id, 1, 'float');  % specdwells 0 - 3.40282e38
rdb.ps_command = fread(id, 1, 'long');  % PS internal use only
rdb.ps_mps_r1 = fread(id, 1, 'long');  % PS MPS R1 setting 1 - 7
rdb.ps_mps_r2 = fread(id, 1, 'long');  % PS MPS R2 setting 1 - 30
rdb.ps_mps_tg = fread(id, 1, 'long');  % PS MPS Transmit gain setting 0 - 200
rdb.ps_mps_freq = fread(id, 1, 'long');  % PS MPS Center frequency hz +-3.40282e38
rdb.ps_aps_r1 = fread(id, 1, 'long');  % PS APS R1 setting 1 - 7
rdb.ps_aps_r2 = fread(id, 1, 'long');  % PS APS R2 setting 1 - 30
rdb.ps_aps_tg = fread(id, 1, 'long');  % PS APS Transmit gain setting 0 - 200
rdb.ps_aps_freq = fread(id, 1, 'long');  % PS APS Center frequency hz +-3.40282e38
rdb.ps_scalei = fread(id, 1, 'float');  % PS rational scaling +-3.40282e38
rdb.ps_scaleq = fread(id, 1, 'float');  % PS unused
rdb.ps_snr_warning = fread(id, 1, 'long');  % PS noise test 0=16 1=32 bits 0, 1
rdb.ps_aps_or_mps = fread(id, 1, 'long');  % PS prescan order logic 0 - 5
rdb.ps_mps_bitmap = fread(id, 1, 'long');  % PS bit mask 4 bits
rdb.ps_powerspec = modchar(fread(id, 256, '*char'));  % PS
rdb.ps_filler1 = fread(id, 1, 'long');  % PS filler
rdb.ps_filler2 = fread(id, 1, 'long');  % PS filler
rdb.rec_noise_mean = fread(id, 16, 'float');  % PS mean noise each receiver +-3.40282e38
rdb.rec_noise_std = fread(id, 16, 'float');  % PS noise calc for muti rec +-3.40282e38
rdb.halfecho = fread(id, 1, 'short');  % spectro full, half echo 0, 1
rdb.im_size_y = fread(id, 1, 'short');  % rh???? 0 - 512
rdb.data_collect_type1 = fread(id, 1, 'long');  % rh???? bit mask 31 bits
rdb.freq_scale = fread(id, 1, 'float');  % rh???? freq k-space step +-3.40282e38
rdb.phase_scale = fread(id, 1, 'float');  % rh???? freq k-space step +-3.40282e38
rdb.ovl = fread(id, 1, 'short');  % rhovl - overlaps for MOTSA
rdb.pclin = fread(id, 1, 'short');  % Linear Corr. 0:off, 1:linear, 2:polynomial
rdb.pclinnpts = fread(id, 1, 'short');  % fit number of points
rdb.pclinorder = fread(id, 1, 'short');  % fit order
rdb.pclinavg = fread(id, 1, 'short');  % linear phase corr avg 0:off, 1:on
rdb.pccon = fread(id, 1, 'short');  % Const Corr. 0:off, 1:Ky spec., 2:polyfit(2/ilv), 3:polyfit(1/ilv)
rdb.pcconnpts = fread(id, 1, 'short');  % fit number of points
rdb.pcconorder = fread(id, 1, 'short');  % fit order
rdb.pcextcorr = fread(id, 1, 'short');  % external correction file 0:don't use, 1: use
rdb.pcgraph = fread(id, 1, 'short');  % Phase Correction coef. image 0:off, 1:linear & constant
rdb.pcileave = fread(id, 1, 'short');  % Interleaves to use for correction: 0=all, 1=only first
rdb.hdbestky = fread(id, 1, 'short');  % bestky view for fractional Ky scan
rdb.pcctrl = fread(id, 1, 'short');  % phase correction research control
rdb.pcthrespts = fread(id, 1, 'short');  % 2..512 adjacent points
rdb.pcdiscbeg = fread(id, 1, 'short');  % 0..512 beginning point to discard
rdb.pcdiscmid = fread(id, 1, 'short');  % 0..512 middle point to discard
rdb.pcdiscend = fread(id, 1, 'short');  % 0..512 ending point to discard
rdb.pcthrespct = fread(id, 1, 'short');  % Threshold percentage
rdb.pcspacial = fread(id, 1, 'short');  % Spacial best ref scan index 0..512
rdb.pctemporal = fread(id, 1, 'short');  % Temporal best ref scan index 0..512
rdb.pcspare = fread(id, 1, 'short');  % spare for phase correction
rdb.ileaves = fread(id, 1, 'short');  % Number of interleaves
rdb.kydir = fread(id, 1, 'short');  % Ky traversal dircetion 0: top-down, 1:center out
rdb.alt = fread(id, 1, 'short');  % Alt read sign 0=no, 1=odd/even, 2=pairs
rdb.reps = fread(id, 1, 'short');  % Number of scan repetitions
rdb.ref = fread(id, 1, 'short');  % Ref Scan 0: off 1: on
rdb.pcconnorm = fread(id, 1, 'float');  % Constant S term normalization factor
rdb.pcconfitwt = fread(id, 1, 'float');  % Constant polyfit weighting factor
rdb.pclinnorm = fread(id, 1, 'float');  % Linear S term normalization factor
rdb.pclinfitwt = fread(id, 1, 'float');  % Linear polyfit weighting factor
rdb.pcbestky = fread(id, 1, 'float');  % Best Ky location
rdb.vrgf = fread(id, 1, 'long');  % control word for VRG filter
rdb.vrgfxres = fread(id, 1, 'long');  % control word for VRGF final x resolution
rdb.bp_corr = fread(id, 1, 'long');  % control word for bandpass asymmetry
rdb.recv_freq_s = fread(id, 1, 'float');  % starting frequency (+62.5)
rdb.recv_freq_e = fread(id, 1, 'float');  % ending frequency (-62.5)
rdb.hniter = fread(id, 1, 'long');  % Selects the number of (continued...)
rdb.fast_rec = fread(id, 1, 'long');  % Added for homodyne II, tells if (continued...)
rdb.refframes = fread(id, 1, 'long');  % total # of frames for ref scan
rdb.refframep = fread(id, 1, 'long');  % # of frames per pass for a ref scan
rdb.scnframe = fread(id, 1, 'long');  % total # of frames for a entire scan
rdb.pasframe = fread(id, 1, 'long');  % # of frames per pass
rdb.user_usage_tag = fread(id, 1, 'ulong');  % for spectro
rdb.user_fill_mapMSW = fread(id, 1, 'ulong');  % for spectro
rdb.user_fill_mapLSW = fread(id, 1, 'ulong');  % for Spectro
rdb.user20 = fread(id, 1, 'float');  % all following usercv are for spectro
rdb.user21 = fread(id, 1, 'float');  %
rdb.user22 = fread(id, 1, 'float');  %
rdb.user23 = fread(id, 1, 'float');  %
rdb.user24 = fread(id, 1, 'float');  %
rdb.user25 = fread(id, 1, 'float');  %
rdb.user26 = fread(id, 1, 'float');  %
rdb.user27 = fread(id, 1, 'float');  %
rdb.user28 = fread(id, 1, 'float');  %
rdb.user29 = fread(id, 1, 'float');  %
rdb.user30 = fread(id, 1, 'float');  %
rdb.user31 = fread(id, 1, 'float');  %
rdb.user32 = fread(id, 1, 'float');  %
rdb.user33 = fread(id, 1, 'float');  %
rdb.user34 = fread(id, 1, 'float');  %
rdb.user35 = fread(id, 1, 'float');  %
rdb.user36 = fread(id, 1, 'float');  %
rdb.user37 = fread(id, 1, 'float');  %
rdb.user38 = fread(id, 1, 'float');  %
rdb.user39 = fread(id, 1, 'float');  %
rdb.user40 = fread(id, 1, 'float');  %
rdb.user41 = fread(id, 1, 'float');  %
rdb.user42 = fread(id, 1, 'float');  %
rdb.user43 = fread(id, 1, 'float');  %
rdb.user44 = fread(id, 1, 'float');  %
rdb.user45 = fread(id, 1, 'float');  %
rdb.user46 = fread(id, 1, 'float');  %
rdb.user47 = fread(id, 1, 'float');  %
rdb.user48 = fread(id, 1, 'float');  %
rdb.pcfitorig = fread(id, 1, 'short');  % Adjust view indexes if set so bestky view = 0
rdb.pcshotfirst = fread(id, 1, 'short');  % First view within an echo group used for fit
rdb.pcshotlast = fread(id, 1, 'short');  % Last view within an echo group used for fit
rdb.pcmultegrp = fread(id, 1, 'short');  % If = 1, force pts from other egrps to be used
rdb.pclinfix = fread(id, 1, 'short');  % If = 2, force slope to be set to pclinslope
rdb.pcconfix = fread(id, 1, 'short');  % If = 2, force slope to be set to pcconslope
rdb.pclinslope = fread(id, 1, 'float');  % Value to set lin slope to if forced
rdb.pcconslope = fread(id, 1, 'float');  % Value to set con slope to if forced
rdb.pccoil = fread(id, 1, 'short');  % If 1,2,3,4, use that coil's results for all
rdb.vvsmode = fread(id, 1, 'short');  % Variable view sharing mode
rdb.vvsaimgs = fread(id, 1, 'short');  % number of original images
rdb.vvstr = fread(id, 1, 'short');  % TR in microseconds
rdb.vvsgender = fread(id, 1, 'short');  % gender: male or female
rdb.zip_factor = fread(id, 1, 'short');  % Slice ZIP factor: 0=OFF, 2, or 4
rdb.maxcoef1a = fread(id, 1, 'float');  % Coefficient A for flow image 1
rdb.maxcoef1b = fread(id, 1, 'float');  % Coefficient B for flow image 1
rdb.maxcoef1c = fread(id, 1, 'float');  % Coefficient C for flow image 1
rdb.maxcoef1d = fread(id, 1, 'float');  % Coefficient D for flow image 1
rdb.maxcoef2a = fread(id, 1, 'float');  % Coefficient A for flow image 2
rdb.maxcoef2b = fread(id, 1, 'float');  % Coefficient B for flow image 2
rdb.maxcoef2c = fread(id, 1, 'float');  % Coefficient C for flow image 2
rdb.maxcoef2d = fread(id, 1, 'float');  % Coefficient D for flow image 2
rdb.maxcoef3a = fread(id, 1, 'float');  % Coefficient A for flow image 3
rdb.maxcoef3b = fread(id, 1, 'float');  % Coefficient B for flow image 3
rdb.maxcoef3c = fread(id, 1, 'float');  % Coefficient C for flow image 3
rdb.maxcoef3d = fread(id, 1, 'float');  % Coefficient D for flow image 3
rdb.ut_ctrl = fread(id, 1, 'long');  % System utility control variable
rdb.dp_type = fread(id, 1, 'short');  % EPI II diffusion control cv
if (ver <= 7.0)  % end of rdbm version 7.0
    rdb.excess = fread(id, 423, 'short');  % free space for later expansion
    return
end
rdb.arw = fread(id, 1, 'short');  % Arrhythmia rejection window(percentage:1-100)
rdb.vps = fread(id, 1, 'short');  % View Per Segment for FastCine
rdb.mcReconEnable = fread(id, 1, 'short');  % N-Coil recon map
rdb.fov = fread(id, 1, 'float');  % Auto-NCoil
rdb.te = fread(id, 1, 'long');  % TE for first echo
rdb.te2 = fread(id, 1, 'long');  % TE for second and later echoes
rdb.dfmrbw = fread(id, 1, 'float');  % BW for navigator frames
rdb.dfmctrl = fread(id, 1, 'long');  % Control flag for dfm (0=off, other=on)
rdb.raw_nex = fread(id, 1, 'long');  % Uncombined NEX at start of recon
rdb.navs_per_pass = fread(id, 1, 'long');  % Max. navigator frames in a pass
rdb.dfmxres = fread(id, 1, 'long');  % xres of navigator frames
rdb.dfmptsize = fread(id, 1, 'long');  % point size of navigator frames
rdb.navs_per_view = fread(id, 1, 'long');  % Num. navigators per frame (tag table)
rdb.dfmdebug = fread(id, 1, 'long');  % control flag for dfm debug
rdb.dfmthreshold = fread(id, 1, 'float');  % threshold for navigator correction
rdb.grid_control = fread(id, 1, 'short');  % bit settings controlling gridding
rdb.b0map = fread(id, 1, 'short');  % B0 map enable and map size
rdb.grid_tediff = fread(id, 1, 'short');  % TE difference between b0 map arms
rdb.grid_motion_comp = fread(id, 1, 'short');  % flag to apply motion compensation
rdb.grid_radius_a = fread(id, 1, 'float');  % variable density transition
rdb.grid_radius_b = fread(id, 1, 'float');  % variable density transition
rdb.grid_max_gradient = fread(id, 1, 'float');  % Max gradient amplitude
rdb.grid_max_slew = fread(id, 1, 'float');  % Max slew rate
rdb.grid_scan_fov = fread(id, 1, 'float');  % Rx scan field of view
rdb.grid_a2d_time = fread(id, 1, 'float');  % A to D sample time microsecs
rdb.grid_density_factor = fread(id, 1, 'float');  % change factor for variable density
rdb.grid_display_fov = fread(id, 1, 'float');  % Rx display field of view
if (ver > 8.0)
    rdb.fatwater = fread(id, 1, 'short');  % for Fat and Water Dual Recon
    rdb.fiestamlf = fread(id, 1, 'short');  % MFO FIESTA recon control bit 16bits
end
rdb.app = fread(id, 1, 'short');  % Auto Post-Processing opcode
rdb.rhncoilsel = fread(id, 1, 'short');  % Auto-Ncoil
rdb.rhncoillimit = fread(id, 1, 'short');  % Auto-Ncoil
rdb.app_option = fread(id, 1, 'short');  % Auto Post_processing options
rdb.grad_mode = fread(id, 1, 'short');  % Gradient mode in Gemini project
if (ver > 8.0)
    rdb.pfile_passes = fread(id, 1, 'short');  % Num passes stored in a multi-pass Pfile (0 means 1 pass)
else
    rdb.edb_hdr_excess2 = fread(id, 1, 'short');  % To make sure that next member is 4-byte aligned
end
rdb.asset = fread(id, 1, 'int');  %
rdb.asset_calthresh = fread(id, 1, 'int');  %
rdb.asset_R = fread(id, 1, 'float');  %
rdb.coilno = fread(id, 1, 'int');  %
rdb.asset_phases = fread(id, 1, 'int');  %
rdb.scancent = fread(id, 1, 'float');  % Table position
rdb.position = fread(id, 1, 'int');  % Patient position
rdb.entry = fread(id, 1, 'int');  % Patient entry
rdb.lmhor = fread(id, 1, 'float');  % Landmark
if (ver <= 8.0)  % end of rdbm version 8.0
    rdb.excess = fread(id, 352, 'short');  % free space for later expansion
    return
end
rdb.last_slice_num = fread(id, 1, 'int');  %
rdb.asset_slice_R = fread(id, 1, 'float');  % Slice reduction factor
rdb.asset_slabwrap = fread(id, 1, 'float');  %
rdb.dwnav_coeff = fread(id, 1, 'float');  % Coeff for amount of phase correction
rdb.dwnav_cor = fread(id, 1, 'short');  % Navigator echo correction
rdb.dwnav_view = fread(id, 1, 'short');  % Num of views of nav echoes
rdb.dwnav_corecho = fread(id, 1, 'short');  % Num of nav echoes for actual correction
rdb.dwnav_sview = fread(id, 1, 'short');  % Start view for phase correction process
rdb.dwnav_eview = fread(id, 1, 'short');  % End view for phase correction process
rdb.dwnav_sshot = fread(id, 1, 'short');  % Start shot for delta phase estimation in nav echoes
rdb.dwnav_eshot = fread(id, 1, 'short');  % End shot for delta phase estimation in nav echoes
rdb.win3d_type = fread(id, 1, 'short');  % 0 = Modified Hanning, 1 = modified Tukey
rdb.win3d_apod = fread(id, 1, 'float');  % degree of apodization; 0.0 = boxcar, 1.0=hanning
rdb.win3d_q = fread(id, 1, 'float');  % apodization at ends, 0.0 = max, 1.0 = boxcar
rdb.ime_scic_enable = fread(id, 1, 'short');  % Surface Coil Intensity Correction: 1 if enabled
rdb.clariview_type = fread(id, 1, 'short');  % Type of Clariview/Name of Filter
rdb.ime_scic_edge = fread(id, 1, 'float');  % Edge paramaters for Enhanced Recon
rdb.ime_scic_smooth = fread(id, 1, 'float');  % Smooth paramaters for Enhanced Recon
rdb.ime_scic_focus = fread(id, 1, 'float');  % Focus paramaters for Enhanced Recon
rdb.clariview_edge = fread(id, 1, 'float');  % Edge paramaters for clariview
rdb.clariview_smooth = fread(id, 1, 'float');  % Smooth paramaters for clariview
rdb.clariview_focus = fread(id, 1, 'float');  % Focus paramaters for clariview
rdb.scic_reduction = fread(id, 1, 'float');  % Reduction paramater for SCIC
rdb.scic_gauss = fread(id, 1, 'float');  % Gauss paramater for SCIC
rdb.scic_threshold = fread(id, 1, 'float');  % Threshold paramater for SCIC
rdb.ectricks_no_regions = fread(id, 1, 'long');  % Total no of regions acquired by PSD
rdb.ectricks_input_regions = fread(id, 1, 'long');  % Total no of input regions for reordering
rdb.psc_reuse = fread(id, 1, 'short');  % Header field for smart prescan
rdb.left_blank = fread(id, 1, 'short');  %
rdb.right_blank = fread(id, 1, 'short');  %
rdb.acquire_type = fread(id, 1, 'short');  % Acquire type information from CV
rdb.retro_control = fread(id, 1, 'short');  % Retrosective FSE phase correction control flag. (continued...)
rdb.etl = fread(id, 1, 'short');  % Added for Retrospective FSE phase correction. This (continued...)
if (ver <= 9.0)  % end of rdbm version 9.0
    rdb.excess = fread(id, 300, 'short');  % free space for later expansion
    return
end
rdb.pcref_start = fread(id, 1, 'short');  % 1st view to use for dynamic EPI phase correction.
rdb.pcref_stop = fread(id, 1, 'short');  % Last view to use for dynamic EPI phase correction.
rdb.ref_skip = fread(id, 1, 'short');  % Number of passes to skip for dynamic EPI phase correction.
rdb.extra_frames_top = fread(id, 1, 'short');  % Number of extra frames at top of K-space
rdb.extra_frames_bot = fread(id, 1, 'short');  % Number of extra frames at bottom of K-space
rdb.multiphase_type = fread(id, 1, 'short');  % 0 = INTERLEAVED , 1 = SEQUENTIAL
rdb.nphases = fread(id, 1, 'short');  % Number of phases in a multiphase scan
rdb.pure = fread(id, 1, 'short');  % PURE flag from psd
rdb.pure_scale = fread(id, 1, 'float');  % Recon scale factor ratio for cal scan
rdb.off_data = fread(id, 1, 'int');  % Byte offset to start of raw data (i.e size of POOL_HEADER)
rdb.off_per_pass = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_per_pass of POOL_HEADER
rdb.off_unlock_raw = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_unlock_raw of POOL_HEADER
rdb.off_data_acq_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_data_acq_tab of POOL_HEADER
rdb.off_nex_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_nex_tab of POOL_HEADER
rdb.off_nex_abort_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_nex_abort_tab of POOL_HEADER
rdb.off_tool = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_tool of POOL_HEADER
rdb.off_exam = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_exam of POOL_HEADER
rdb.off_series = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_series of POOL_HEADER
rdb.off_image = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_image of POOL_HEADER
rdb.off_ps = fread(id, 1, 'int');  % spare
rdb.off_spare_b = fread(id, 1, 'int');  % spare
rdb.new_wnd_level_flag = fread(id, 1, 'int');  % New WW/WL algo enable/disable flag
rdb.wnd_image_hist_area = fread(id, 1, 'int');  % Image Area %
rdb.wnd_high_hist = fread(id, 1, 'float');  % Histogram Area Top
rdb.wnd_lower_hist = fread(id, 1, 'float');  % Histogram Area Bottom
rdb.pure_filter = fread(id, 1, 'short');  % PURE noise reduction on=1/off=0
rdb.cfg_pure_filter = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_fit_order = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_kernelsize_z = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_kernelsize_xy = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_weight_radius = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_intensity_scale = fread(id, 1, 'short');  % PURE cfg file value
rdb.cfg_pure_noise_threshold = fread(id, 1, 'short');  % PURE cfg file value
rdb.wienera = fread(id, 1, 'float');  % NB maintain alignment of floats
rdb.wienerb = fread(id, 1, 'float');  %
rdb.wienert2 = fread(id, 1, 'float');  %
rdb.wieneresp = fread(id, 1, 'float');  %
rdb.wiener = fread(id, 1, 'short');  %
rdb.flipfilter = fread(id, 1, 'short');  %
rdb.dbgrecon = fread(id, 1, 'short');  %
rdb.ech2skip = fread(id, 1, 'short');  %
rdb.tricks_type = fread(id, 1, 'int');  % 0 = Subtracted, 1 = Unsubtracted
rdb.lcfiesta_phase = fread(id, 1, 'float');  % LC Fiesta
rdb.lcfiesta = fread(id, 1, 'short');  % LC Fiesta
rdb.herawflt = fread(id, 1, 'short');  % Half echo raw data filter
rdb.herawflt_befnwin = fread(id, 1, 'short');  % Half echo raw data filter
rdb.herawflt_befntran = fread(id, 1, 'short');  % Half echo raw data filter
rdb.herawflt_befamp = fread(id, 1, 'float');  % Half echo raw data filter
rdb.herawflt_hpfamp = fread(id, 1, 'float');  % Half echo raw data filter
rdb.heover = fread(id, 1, 'short');  % Half echo over sampling
rdb.pure_correction_threshold = fread(id, 1, 'short');  % PURE Correction threshold
    
rdb.swiftenable  = fread(id, 1, 'int');
rdb.numslabs = fread(id, 1, 'short');
rdb.numCoilConfigs = fread(id, 1, 'ushort');
rdb.ps_autoshim_status = fread(id, 1, 'int');
rdb.dynaplan_numphases = fread(id, 1, 'int');
rdb.medal_cfg = fread(id, 1, 'short');
rdb.medal_nstack = fread(id, 1, 'short');
rdb.edal_echo_order = fread(id, 1, 'short');
rdb.medal_kernel_up = fread(id, 1, 'short');
rdb.medal_kernel_down = fread(id, 1, 'short');
rdb.medal_kernel_smooth = fread(id, 1, 'short');
rdb.medal_start = fread(id, 1, 'short');
rdb.medal_end = fread(id, 1, 'short');
rdb.rcideal = fread(id, 1, 'int');
rdb.rcdixproc = fread(id, 1, 'int');
rdb.df = fread(id, 1, 'float');
rdb.bw = fread(id, 1, 'float');
rdb.te1 = fread(id, 1, 'float');
rdb.esp = fread(id, 1, 'float');
rdb.feextra = fread(id, 1, 'int');
rdb.raw_pass_size = fread(id, 1, 'ubit64');
rdb.sspsave = fread(id, 1, 'ubit64');
rdb.udasave = fread(id, 1, 'ubit64');
rdb.vibrant = fread(id, 1, 'short');
rdb.asset_torso = fread(id, 1, 'short');
rdb.asset_alt_cal = fread(id, 1, 'int');
rdb.hdr_kacq_uid = fread(id, 1, 'int');
rdb.cttEntry = modchar(fread(id, 944, '*char'));  
rdb.psc_ta = fread(id, 1, 'long');
rdb.disk_acq_ctrl = fread(id, 1, 'long');
rdb.asset_localTx = fread(id, 1, 'int');
rdb.scale3d = fread(id, 1, 'float');
rdb.broad_band_select = fread(id, 1, 'int');
rdb.scanner_mode = fread(id, 1, 'short');
rdb.padding_1 = fread(id, 1, 'short');
rdb.excess = fread(id, 716, 'short');

   
end


function exam = read_exam(id, ver);
if (ver < 11.0)
    exam.ex_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Exam
    exam.ex_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    exam.ex_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Exam
    fseek(id, 1, 0);  % for struct byte alignment
    exam.ex_no = fread(id, 1, 'ushort');  % Exam Number
    exam.hospname = modchar(fread(id, 33, '*char'));  % Hospital Name
    fseek(id, 1, 0);  % for struct byte alignment
    exam.detect = fread(id, 1, 'short');  % Detector Type
    fseek(id, 2, 0);  % for struct byte alignment
    exam.numcells = fread(id, 1, 'int');  % Number of cells in det
    exam.zerocell = fread(id, 1, 'float');  % Cell number at theta
    exam.cellspace = fread(id, 1, 'float');  % Cell spacing
    exam.srctodet = fread(id, 1, 'float');  % Distance from source to detector
    exam.srctoiso = fread(id, 1, 'float');  % Distance from source to iso
    exam.tubetyp = fread(id, 1, 'short');  % Tube type
    exam.dastyp = fread(id, 1, 'short');  % DAS type
    exam.num_dcnk = fread(id, 1, 'short');  % Number of Decon Kernals
    exam.dcn_len = fread(id, 1, 'short');  % Number of elements in a Decon Kernal
    exam.dcn_density = fread(id, 1, 'short');  % Decon Kernal density
    exam.dcn_stepsize = fread(id, 1, 'short');  % Decon Kernal stepsize
    exam.dcn_shiftcnt = fread(id, 1, 'short');  % Decon Kernal Shift Count
    fseek(id, 2, 0);  % for struct byte alignment
    exam.magstrength = fread(id, 1, 'int');  % Magnet strength (in gauss)
    exam.patid = modchar(fread(id, 13, '*char'));  % Patient ID for this Exam
    exam.patname = modchar(fread(id, 25, '*char'));  % Patient Name
    exam.patage = fread(id, 1, 'short');  % Patient Age (years, months or days)
    exam.patian = fread(id, 1, 'short');  % Patient Age Notation
    exam.patsex = fread(id, 1, 'short');  % Patient Sex
    exam.patweight = fread(id, 1, 'int');  % Patient Weight
    exam.trauma = fread(id, 1, 'short');  % Trauma Flag
    exam.hist = modchar(fread(id, 61, '*char'));  % Patient History
    exam.reqnum = modchar(fread(id, 13, '*char'));  % Requisition Number
    exam.ex_datetime = fread(id, 1, 'int');  % Exam date/time stamp
    exam.refphy = modchar(fread(id, 33, '*char'));  % Referring Physician
    exam.diagrad = modchar(fread(id, 33, '*char'));  % Diagnostician/Radiologist
    exam.op = modchar(fread(id, 4, '*char'));  % Operator
    exam.ex_desc = modchar(fread(id, 23, '*char'));  % Exam Description
    exam.ex_typ = modchar(fread(id, 3, '*char'));  % Exam Type
    exam.ex_format = fread(id, 1, 'short');  % Exam Format
    exam.padding1 = modchar(fread(id, 4, '*char'));  % Padding for 32bit alignment
    fseek(id, 2, 0);  % for struct byte alignment
    exam.firstaxtime = fread(id, 1, 'double');  % Start time(secs) of first axial in exam
    exam.ex_sysid = modchar(fread(id, 9, '*char'));  % Creator Suite and Host
    fseek(id, 3, 0);  % for struct byte alignment
    exam.ex_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    exam.protocolflag = fread(id, 1, 'short');  % Non-Zero indicates Protocol Exam
    exam.ex_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    fseek(id, 1, 0);  % for struct byte alignment
    exam.ex_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    exam.ex_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    exam.ex_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    exam.ex_checksum = fread(id, 1, 'ulong');  % Exam Record Checksum
    exam.ex_complete = fread(id, 1, 'long');  % Exam Complete Flag
    exam.ex_seriesct = fread(id, 1, 'long');  % Last Series Number Used
    exam.ex_numarch = fread(id, 1, 'long');  % Number of Series Archived
    exam.ex_numseries = fread(id, 1, 'long');  % Number of Series Existing
    exam.ex_series = read_vartype(id);  % kluge for VARTYPE
    exam.ex_numunser = fread(id, 1, 'long');  % Number of Unstored Series
    exam.ex_unseries = read_vartype(id);  % kluge for VARTYPE
    exam.ex_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Series
    exam.ex_toarchive = read_vartype(id);  % kluge for VARTYPE
    exam.ex_prospcnt = fread(id, 1, 'long');  % Number of Prospective/Scout Series
    exam.ex_prosp = read_vartype(id);  % kluge for VARTYPE
    exam.ex_modelnum = fread(id, 1, 'long');  % Last Model Number used
    exam.ex_modelcnt = fread(id, 1, 'long');  % Number of ThreeD Models
    exam.ex_models = read_vartype(id);  % kluge for VARTYPE
    exam.ex_stat = fread(id, 1, 'short');  % Patient Status
    exam.uniq_sys_id = modchar(fread(id, 16, '*char'));  % Unique System ID
    exam.service_id = modchar(fread(id, 16, '*char'));  % Unique Service ID
    exam.mobile_loc = modchar(fread(id, 4, '*char'));  % Mobile Location Number
    exam.study_uid = moduid(fread(id, 32, 'uint8'));  % Study Entity Unique ID
    exam.study_status = fread(id, 1, 'short');  % indicates if study has complete info(DICOM/genesis)
    exam.refsopcuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    exam.refsopiuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    exam.patnameff = modchar(fread(id, 65, '*char'));  % FF Patient Name
    exam.patidff = modchar(fread(id, 65, '*char'));  % FF Patient ID
    exam.reqnumff = modchar(fread(id, 17, '*char'));  % FF Requisition No
    exam.dateofbirth = modchar(fread(id, 9, '*char'));  % Date of Birth
    exam.mwlstudyuid = moduid(fread(id, 32, 'uint8'));  % Genesis Exam UID
    exam.mwlstudyid = modchar(fread(id, 16, '*char'));  % Genesis Exam No
    exam.ex_padding = modchar(fread(id, 248, '*char'));  % Spare Space
    exam.padding2 = modchar(fread(id, 4, '*char'));  % Padding for 32bit alignment
elseif (ver == 11.0)
    exam.firstaxtime = fread(id, 1, 'double');  % Start time(secs) of first axial in exam
    exam.ex_series = read_vartype(id);  % kluge for VARTYPE
    exam.ex_unseries = read_vartype(id);  % kluge for VARTYPE
    exam.ex_toarchive = read_vartype(id);  % kluge for VARTYPE
    exam.ex_prosp = read_vartype(id);  % kluge for VARTYPE
    exam.ex_models = read_vartype(id);  % kluge for VARTYPE
    exam.zerocell = fread(id, 1, 'float');  % Cell number at theta
    exam.cellspace = fread(id, 1, 'float');  % Cell spacing
    exam.srctodet = fread(id, 1, 'float');  % Distance from source to detector
    exam.srctoiso = fread(id, 1, 'float');  % Distance from source to iso
    exam.ex_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    exam.ex_complete = fread(id, 1, 'long');  % Exam Complete Flag
    exam.ex_seriesct = fread(id, 1, 'long');  % Last Series Number Used
    exam.ex_numarch = fread(id, 1, 'long');  % Number of Series Archived
    exam.ex_numseries = fread(id, 1, 'long');  % Number of Series Existing
    exam.ex_numunser = fread(id, 1, 'long');  % Number of Unstored Series
    exam.ex_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Series
    exam.ex_prospcnt = fread(id, 1, 'long');  % Number of Prospective/Scout Series
    exam.ex_modelnum = fread(id, 1, 'long');  % Last Model Number used
    exam.ex_modelcnt = fread(id, 1, 'long');  % Number of ThreeD Models
    exam.ex_checksum = fread(id, 1, 'ulong');  % Exam Record Checksum
    exam.numcells = fread(id, 1, 'int');  % Number of cells in det
    exam.magstrength = fread(id, 1, 'int');  % Magnet strength (in gauss)
    exam.patweight = fread(id, 1, 'int');  % Patient Weight
    exam.ex_datetime = fread(id, 1, 'int');  % Exam date/time stamp
    exam.ex_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    exam.ex_no = fread(id, 1, 'ushort');  % Exam Number
    exam.ex_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    exam.detect = fread(id, 1, 'short');  % Detector Type
    exam.tubetyp = fread(id, 1, 'short');  % Tube type
    exam.dastyp = fread(id, 1, 'short');  % DAS type
    exam.num_dcnk = fread(id, 1, 'short');  % Number of Decon Kernals
    exam.dcn_len = fread(id, 1, 'short');  % Number of elements in a Decon Kernal
    exam.dcn_density = fread(id, 1, 'short');  % Decon Kernal density
    exam.dcn_stepsize = fread(id, 1, 'short');  % Decon Kernal stepsize
    exam.dcn_shiftcnt = fread(id, 1, 'short');  % Decon Kernal Shift Count
    exam.patage = fread(id, 1, 'short');  % Patient Age (years, months or days)
    exam.patian = fread(id, 1, 'short');  % Patient Age Notation
    exam.patsex = fread(id, 1, 'short');  % Patient Sex
    exam.ex_format = fread(id, 1, 'short');  % Exam Format
    exam.trauma = fread(id, 1, 'short');  % Trauma Flag
    exam.protocolflag = fread(id, 1, 'short');  % Non-Zero indicates Protocol Exam
    exam.study_status = fread(id, 1, 'short');  % indicates if study has complete info(DICOM/genesis)
    exam.padding = fread(id, 3, 'short');  %
    exam.hist = modchar(fread(id, 61, '*char'));  % Patient History
    exam.reqnum = modchar(fread(id, 13, '*char'));  % Requisition Number
    exam.refphy = modchar(fread(id, 33, '*char'));  % Referring Physician
    exam.diagrad = modchar(fread(id, 33, '*char'));  % Diagnostician/Radiologist
    exam.op = modchar(fread(id, 4, '*char'));  % Operator
    exam.ex_desc = modchar(fread(id, 65, '*char'));  % Exam Description
    exam.ex_typ = modchar(fread(id, 3, '*char'));  % Exam Type
    exam.ex_sysid = modchar(fread(id, 9, '*char'));  % Creator Suite and Host
    exam.ex_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    exam.ex_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Exam
    exam.hospname = modchar(fread(id, 33, '*char'));  % Hospital Name
    exam.patid = modchar(fread(id, 13, '*char'));  % Patient ID for this Exam
    exam.patname = modchar(fread(id, 25, '*char'));  % Patient Name
    exam.ex_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Exam
    exam.ex_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    exam.ex_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    exam.uniq_sys_id = modchar(fread(id, 16, '*char'));  % Unique System ID
    exam.service_id = modchar(fread(id, 16, '*char'));  % Unique Service ID
    exam.mobile_loc = modchar(fread(id, 4, '*char'));  % Mobile Location Number
    exam.study_uid = moduid(fread(id, 32, 'uint8'));  % Study Entity Unique ID
    exam.refsopcuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    exam.refsopiuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    exam.patnameff = modchar(fread(id, 65, '*char'));  % FF Patient Name
    exam.patidff = modchar(fread(id, 65, '*char'));  % FF Patient ID
    exam.reqnumff = modchar(fread(id, 17, '*char'));  % FF Requisition No
    exam.dateofbirth = modchar(fread(id, 9, '*char'));  % Date of Birth
    exam.mwlstudyuid = moduid(fread(id, 32, 'uint8'));  % Genesis Exam UID
    exam.mwlstudyid = modchar(fread(id, 16, '*char'));  % Genesis Exam No
    exam.ex_padding = modchar(fread(id, 222, '*char'));  % Spare Space    
    
elseif (ver >= 14.0 & ver < 20)
    exam.firstaxtime = fread(id, 1, 'double');  % Start time(secs) of first axial in exam
    exam.double_padding = modchar(fread(id, 72, 'char')); % double padding
    
    exam.zerocell = fread(id, 1, 'float');  % Cell number at theta
    exam.cellspace = fread(id, 1, 'float');  % Cell spacing
    exam.srctodet = fread(id, 1, 'float');  % Distance from source to detector
    exam.srctoiso = fread(id, 1, 'float');  % Distance from source to iso    
    exam.float_padding = modchar(fread(id, 32, 'char')); % float padding
    
    exam.ex_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    exam.ex_complete = fread(id, 1, 'long');  % Exam Complete Flag
    exam.ex_seriesct = fread(id, 1, 'long');  % Last Series Number Used
    exam.ex_numarch = fread(id, 1, 'long');  % Number of Series Archived
    exam.ex_numseries = fread(id, 1, 'long');  % Number of Series Existing
    exam.ex_numunser = fread(id, 1, 'long');  % Number of Unstored Series
    exam.ex_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Series
    exam.ex_prospcnt = fread(id, 1, 'long');  % Number of Prospective/Scout Series
    exam.ex_modelnum = fread(id, 1, 'long');  % Last Model Number used
    exam.ex_modelcnt = fread(id, 1, 'long');  % Number of ThreeD Models
    exam.ex_checksum = fread(id, 1, 'ulong');  % Exam Record Checksum
    exam.long_padding = modchar(fread(id, 32, 'char')); % long padding
    
    exam.numcells = fread(id, 1, 'int');  % Number of cells in det
    exam.magstrength = fread(id, 1, 'int');  % Magnet strength (in gauss)
    exam.patweight = fread(id, 1, 'int');  % Patient Weight
    exam.ex_datetime = fread(id, 1, 'int');  % Exam date/time stamp
    exam.ex_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    exam.int_padding = modchar(fread(id, 48, 'char')); % long padding
        
    exam.ex_no = fread(id, 1, 'ushort');  % Exam Number
    exam.ex_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    exam.detect = fread(id, 1, 'short');  % Detector Type
    exam.tubetyp = fread(id, 1, 'short');  % Tube type
    exam.dastyp = fread(id, 1, 'short');  % DAS type
    exam.num_dcnk = fread(id, 1, 'short');  % Number of Decon Kernals
    exam.dcn_len = fread(id, 1, 'short');  % Number of elements in a Decon Kernal
    exam.dcn_density = fread(id, 1, 'short');  % Decon Kernal density
    exam.dcn_stepsize = fread(id, 1, 'short');  % Decon Kernal stepsize
    exam.dcn_shiftcnt = fread(id, 1, 'short');  % Decon Kernal Shift Count
    exam.patage = fread(id, 1, 'short');  % Patient Age (years, months or days)
    exam.patian = fread(id, 1, 'short');  % Patient Age Notation
    exam.patsex = fread(id, 1, 'short');  % Patient Sex
    exam.ex_format = fread(id, 1, 'short');  % Exam Format
    exam.trauma = fread(id, 1, 'short');  % Trauma Flag
    exam.protocolflag = fread(id, 1, 'short');  % Non-Zero indicates Protocol Exam
    exam.study_status = fread(id, 1, 'short');  % indicates if study has complete info(DICOM/genesis)
    exam.short_padding = modchar(fread(id, 22, 'char'));  %
    
    exam.hist = modchar(fread(id, 61, '*char'));  % Patient History
    exam.reqnum = modchar(fread(id, 13, '*char'));  % Requisition Number
    exam.refphy = modchar(fread(id, 33, '*char'));  % Referring Physician
    exam.diagrad = modchar(fread(id, 33, '*char'));  % Diagnostician/Radiologist
    exam.op = modchar(fread(id, 4, '*char'));  % Operator
    exam.ex_desc = modchar(fread(id, 65, '*char'));  % Exam Description
    exam.ex_typ = modchar(fread(id, 3, '*char'));  % Exam Type
    exam.ex_sysid = modchar(fread(id, 9, '*char'));  % Creator Suite and Host
    exam.ex_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    exam.ex_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Exam
    exam.hospname = modchar(fread(id, 33, '*char'));  % Hospital Name
    exam.patid = modchar(fread(id, 13, '*char'));  % Patient ID for this Exam
    exam.patname = modchar(fread(id, 25, '*char'));  % Patient Name
    exam.ex_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Exam
    exam.ex_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    exam.ex_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    exam.uniq_sys_id = modchar(fread(id, 16, '*char'));  % Unique System ID
    exam.service_id = modchar(fread(id, 16, '*char'));  % Unique Service ID
    exam.mobile_loc = modchar(fread(id, 4, '*char'));  % Mobile Location Number
    exam.study_uid = moduid(fread(id, 32, 'uint8'));  % Study Entity Unique ID
    exam.refsopcuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    exam.refsopiuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    exam.patnameff = modchar(fread(id, 65, '*char'));  % FF Patient Name
    exam.patidff = modchar(fread(id, 65, '*char'));  % FF Patient ID
    exam.reqnumff = modchar(fread(id, 17, '*char'));  % FF Requisition No
    exam.dateofbirth = modchar(fread(id, 9, '*char'));  % Date of Birth
    exam.mwlstudyuid = moduid(fread(id, 32, 'uint8'));  % Genesis Exam UID
    exam.mwlstudyid = modchar(fread(id, 16, '*char'));  % Genesis Exam No
    exam.ex_padding = modchar(fread(id, 62, '*char'));  % ex_padding 
    
    
elseif (ver >= 20)   % 20.0
    exam.firstaxtime = fread(id, 1, 'double');  % Start time(secs) of first axial in exam
    exam.double_padding = modchar(fread(id, 248, 'char')); % double padding
    
    exam.zerocell = fread(id, 1, 'float');  % Cell number at theta
    exam.cellspace = fread(id, 1, 'float');  % Cell spacing
    exam.srctodet = fread(id, 1, 'float');  % Distance from source to detector
    exam.srctoiso = fread(id, 1, 'float');  % Distance from source to iso    
    exam.float_padding = modchar(fread(id, 128, 'char')); % float padding
    
    exam.ex_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    exam.ex_complete = fread(id, 1, 'long');  % Exam Complete Flag
    exam.ex_seriesct = fread(id, 1, 'long');  % Last Series Number Used
    exam.ex_numarch = fread(id, 1, 'long');  % Number of Series Archived
    exam.ex_numseries = fread(id, 1, 'long');  % Number of Series Existing
    exam.ex_numunser = fread(id, 1, 'long');  % Number of Unstored Series
    exam.ex_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Series
    exam.ex_prospcnt = fread(id, 1, 'long');  % Number of Prospective/Scout Series
    exam.ex_modelnum = fread(id, 1, 'long');  % Last Model Number used
    exam.ex_modelcnt = fread(id, 1, 'long');  % Number of ThreeD Models    
    exam.int_padding1 = modchar(fread(id, 128, 'char')); % int padding1
    
    %  Changed from long to int for 64-bit recon      
    exam.numcells = fread(id, 1, 'int');  % Number of cells in det
    exam.magstrength = fread(id, 1, 'int');  % Magnet strength (in gauss)
    exam.patweight = fread(id, 1, 'int');  % Patient Weight
    exam.ex_datetime = fread(id, 1, 'int');  % Exam date/time stamp
    exam.ex_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change    
    exam.int_padding2 = modchar(fread(id, 108, 'char')); % int padding2
    
    exam.ex_no = fread(id, 1, 'ushort');  % Exam Number
    exam.ex_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    exam.detect = fread(id, 1, 'short');  % Detector Type
    exam.tubetyp = fread(id, 1, 'short');  % Tube type
    exam.dastyp = fread(id, 1, 'short');  % DAS type
    exam.num_dcnk = fread(id, 1, 'short');  % Number of Decon Kernals
    exam.dcn_len = fread(id, 1, 'short');  % Number of elements in a Decon Kernal
    exam.dcn_density = fread(id, 1, 'short');  % Decon Kernal density
    exam.dcn_stepsize = fread(id, 1, 'short');  % Decon Kernal stepsize
    exam.dcn_shiftcnt = fread(id, 1, 'short');  % Decon Kernal Shift Count
    exam.patage = fread(id, 1, 'short');  % Patient Age (years, months or days)
    exam.patian = fread(id, 1, 'short');  % Patient Age Notation
    exam.patsex = fread(id, 1, 'short');  % Patient Sex
    exam.ex_format = fread(id, 1, 'short');  % Exam Format
    exam.trauma = fread(id, 1, 'short');  % Trauma Flag
    exam.protocolflag = fread(id, 1, 'short');  % Non-Zero indicates Protocol Exam
    exam.study_status = fread(id, 1, 'short');  % indicates if study has complete info(DICOM/genesis)    
    exam.short_padding = modchar(fread(id, 70, 'char')); % short padding

    exam.hist = modchar(fread(id, 257, '*char'));  % Patient History
    exam.refphy = modchar(fread(id, 65, '*char'));  % Referring Physician
    exam.diagrad = modchar(fread(id, 65, '*char'));  % Diagnostician/Radiologist
    
    % moved here from series struct
    exam.operator_new = modchar(fread(id, 65, '*char'));  % Operator_new
    exam.ex_desc = modchar(fread(id, 65, '*char'));  % Exam Description
    exam.ex_typ = modchar(fread(id, 3, '*char'));  % Exam Type
    exam.ex_sysid = modchar(fread(id, 9, '*char'));  % Creator Suite and Host
    exam.ex_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    exam.ex_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Exam
    exam.hospname = modchar(fread(id, 33, '*char'));  % Hospital Name       
    exam.ex_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Exam
    exam.ex_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    exam.ex_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    exam.uniq_sys_id = modchar(fread(id, 16, '*char'));  % Unique System ID
    exam.service_id = modchar(fread(id, 16, '*char'));  % Unique Service ID
    exam.mobile_loc = modchar(fread(id, 4, '*char'));  % Mobile Location Number
    exam.study_uid = moduid(fread(id, 32, 'uint8'));  % Study Entity Unique ID
    exam.refsopcuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    exam.refsopiuid = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    
    % Part of Ref Study Seq    
    exam.patnameff = modchar(fread(id, 65, '*char'));  % FF Patient Name
    exam.patidff = modchar(fread(id, 65, '*char'));  % FF Patient ID
    exam.reqnumff = modchar(fread(id, 17, '*char'));  % FF Requisition No
    exam.dateofbirth = modchar(fread(id, 9, '*char'));  % Date of Birth
    exam.mwlstudyuid = moduid(fread(id, 32, 'uint8'));  % Genesis Exam UID
    exam.mwlstudyid = modchar(fread(id, 16, '*char'));  % Genesis Exam No
    exam.ex_padding = modchar(fread(id, 240, '*char'));  % Spare Space
end
end

function series = read_series(id, ver)
if (ver < 9.0)
    series.se_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Series
    series.se_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    series.se_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Series
    fseek(id, 1, 0);  % for struct byte alignment
    series.se_exno = fread(id, 1, 'ushort');  % Exam Number
    series.se_no = fread(id, 1, 'short');  % Series Number
    series.se_datetime = fread(id, 1, 'int');  % Allocation Series Data/Time stamp
    series.se_actual_dt = fread(id, 1, 'int');  % Actual Series Data/Time stamp
    series.se_desc = modchar(fread(id, 30, '*char'));  % Series Description
    series.pr_sysid = modchar(fread(id, 9, '*char'));  % Primary Receiver Suite and Host
    series.pansysid = modchar(fread(id, 9, '*char'));  % Archiver Suite and Host
    series.se_typ = fread(id, 1, 'short');  % Series Type
    series.se_source = fread(id, 1, 'short');  % Series from which prescribed
    series.se_plane = fread(id, 1, 'short');  % Most-like Plane (for L/S)
    series.scan_type = fread(id, 1, 'short');  % Scout or Axial (for CT)
    series.position = fread(id, 1, 'int');  % Patient Position
    series.entry = fread(id, 1, 'int');  % Patient Entry
    series.anref = modchar(fread(id, 3, '*char'));  % Anatomical reference
    fseek(id, 1, 0);  % for struct byte alignment
    series.lmhor = fread(id, 1, 'float');  % Horizontal Landmark
    series.prtcl = modchar(fread(id, 25, '*char'));  % Scan Protocol Name
    fseek(id, 1, 0);  % for struct byte alignment
    series.se_contrast = fread(id, 1, 'short');  % Non-zero if > 0 image used contrast(L/S)
    series.start_ras = modchar(fread(id, 1, '*char'));  % RAS letter for first scan location (L/S)
    fseek(id, 3, 0);  % for struct byte alignment
    series.start_loc = fread(id, 1, 'float');  % First scan location (L/S)
    series.end_ras = modchar(fread(id, 1, '*char'));  % RAS letter for last scan location (L/S)
    fseek(id, 3, 0);  % for struct byte alignment
    series.end_loc = fread(id, 1, 'float');  % Last scan location (L/S)
    series.se_pseq = fread(id, 1, 'short');  % Last Pulse Sequence Used (L/S)
    series.se_sortorder = fread(id, 1, 'short');  % Image Sort Order (L/S)
    series.se_lndmrkcnt = fread(id, 1, 'int');  % Landmark Counter
    series.se_nacq = fread(id, 1, 'short');  % Number of Acquisitions
    series.xbasest = fread(id, 1, 'short');  % Starting number for baselines
    series.xbaseend = fread(id, 1, 'short');  % Ending number for baselines
    series.xenhst = fread(id, 1, 'short');  % Starting number for enhanced scans
    series.xenhend = fread(id, 1, 'short');  % Ending number for enhanced scans
    fseek(id, 2, 0);  % for struct byte alignment
    series.se_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    series.se_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    fseek(id, 3, 0);  % for struct byte alignment
    series.se_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    series.se_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    series.se_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    series.se_pds_a = fread(id, 1, 'float');  % PixelData size - as stored
    series.se_pds_c = fread(id, 1, 'float');  % PixelData size - Compressed
    series.se_pds_u = fread(id, 1, 'float');  % PixelData size - UnCompressed
    series.se_checksum = fread(id, 1, 'ulong');  % Series Record checksum
    series.se_complete = fread(id, 1, 'long');  % Series Complete Flag
    series.se_numarch = fread(id, 1, 'long');  % Number of Images Archived
    series.se_imagect = fread(id, 1, 'long');  % Last Image Number Used
    series.se_numimages = fread(id, 1, 'long');  % Number of Images Existing
    series.se_images = read_vartype(id);  % kluge for VARTYPE
    series.se_numunimg = fread(id, 1, 'long');  % Number of Unstored Images
    series.se_unimages = read_vartype(id);  % kluge for VARTYPE
    series.se_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Images
    series.se_toarchive = read_vartype(id);  % kluge for VARTYPE
    series.echo1_alpha = fread(id, 1, 'float');  % Echo 1 Alpha Value
    series.echo1_beta = fread(id, 1, 'float');  % Echo 1 Beta Value
    series.echo1_window = fread(id, 1, 'ushort');  % Echo 1 Window Value
    series.echo1_level = fread(id, 1, 'short');  % Echo 1 Level Value
    series.echo2_alpha = fread(id, 1, 'float');  % Echo 2 Alpha Value
    series.echo2_beta = fread(id, 1, 'float');  % Echo 2 Beta Value
    series.echo2_window = fread(id, 1, 'ushort');  % Echo 2 Window Value
    series.echo2_level = fread(id, 1, 'short');  % Echo 2 Level Value
    series.echo3_alpha = fread(id, 1, 'float');  % Echo 3 Alpha Value
    series.echo3_beta = fread(id, 1, 'float');  % Echo 3 Beta Value
    series.echo3_window = fread(id, 1, 'ushort');  % Echo 3 Window Value
    series.echo3_level = fread(id, 1, 'short');  % Echo 3 Level Value
    series.echo4_alpha = fread(id, 1, 'float');  % Echo 4 Alpha Value
    series.echo4_beta = fread(id, 1, 'float');  % Echo 4 Beta Value
    series.echo4_window = fread(id, 1, 'ushort');  % Echo 4 Window Value
    series.echo4_level = fread(id, 1, 'short');  % Echo 4 Level Value
    series.echo5_alpha = fread(id, 1, 'float');  % Echo 5 Alpha Value
    series.echo5_beta = fread(id, 1, 'float');  % Echo 5 Beta Value
    series.echo5_window = fread(id, 1, 'ushort');  % Echo 5 Window Value
    series.echo5_level = fread(id, 1, 'short');  % Echo 5 Level Value
    series.echo6_alpha = fread(id, 1, 'float');  % Echo 6 Alpha Value
    series.echo6_beta = fread(id, 1, 'float');  % Echo 6 Beta Value
    series.echo6_window = fread(id, 1, 'ushort');  % Echo 6 Window Value
    series.echo6_level = fread(id, 1, 'short');  % Echo 6 Level Value
    series.echo7_alpha = fread(id, 1, 'float');  % Echo 7 Alpha Value
    series.echo7_beta = fread(id, 1, 'float');  % Echo 7 Beta Value
    series.echo7_window = fread(id, 1, 'ushort');  % Echo 7 Window Value
    series.echo7_level = fread(id, 1, 'short');  % Echo 7 Level Value
    series.echo8_alpha = fread(id, 1, 'float');  % Echo 8 Alpha Value
    series.echo8_beta = fread(id, 1, 'float');  % Echo 8 Beta Value
    series.echo8_window = fread(id, 1, 'ushort');  % Echo 8 Window Value
    series.echo8_level = fread(id, 1, 'short');  % Echo 8 Level Value
    series.series_uid = moduid(fread(id, 32, 'uint8'));  % Series Entity Unique ID
    series.landmark_uid = moduid(fread(id, 32, 'uint8'));  % Landmark Unique ID
    series.equipmnt_uid = moduid(fread(id, 32, 'uint8'));  % Equipment Unique ID
    series.refsopcuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    series.refsopiuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    series.schacitval = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Value
    series.schacitdesc = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Description
    series.schacitmea = modchar(fread(id, 64, '*char'));  % Sched Proc Action Item Seq - Meaning
    series.schprocstdesc = modchar(fread(id, 64, '*char'));  % Sched Proc Step Desc
    series.schprocstid = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 1
    series.reqprocstid = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 1
    series.perprocstid = modchar(fread(id, 16, '*char'));  % PPS ID
    series.perprocstdesc = modchar(fread(id, 64, '*char'));  % PPS Description
    series.table_entry = fread(id, 1, 'short');  % Table position for nMR and iMR
    series.SwingAngle = fread(id, 1, 'short');  % nMR - Swing Angle
    series.LateralOffset = fread(id, 1, 'short');  % nMR - Offset
    series.reqprocstid2 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 2
    series.reqprocstid3 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 3
    series.schprocstid2 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 2
    series.schprocstid3 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 3
    series.refImgUID = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image
    series.GradientCoil = fread(id, 1, 'short');  % Gradient Coil Selection
    series.se_padding = modchar(fread(id, 148, '*char'));  % Spare Space
elseif (ver >= 9.0 && ver < 11.0)
    series.se_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Series
    series.se_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    series.se_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Series
    fseek(id, 1, 0);  % for struct byte alignment
    series.se_exno = fread(id, 1, 'ushort');  % Exam Number
    series.se_no = fread(id, 1, 'short');  % Series Number
    series.se_datetime = fread(id, 1, 'int');  % Allocation Series Data/Time stamp
    series.se_actual_dt = fread(id, 1, 'int');  % Actual Series Data/Time stamp
    series.se_desc = modchar(fread(id, 30, '*char'));  % Series Description
    series.pr_sysid = modchar(fread(id, 9, '*char'));  % Primary Receiver Suite and Host
    series.pansysid = modchar(fread(id, 9, '*char'));  % Archiver Suite and Host
    series.se_typ = fread(id, 1, 'short');  % Series Type
    series.se_source = fread(id, 1, 'short');  % Series from which prescribed
    series.se_plane = fread(id, 1, 'short');  % Most-like Plane (for L/S)
    series.scan_type = fread(id, 1, 'short');  % Scout or Axial (for CT)
    series.position = fread(id, 1, 'int');  % Patient Position
    series.entry = fread(id, 1, 'int');  % Patient Entry
    series.anref = modchar(fread(id, 3, '*char'));  % Anatomical reference
    fseek(id, 1, 0);  % for struct byte alignment
    series.lmhor = fread(id, 1, 'float');  % Horizontal Landmark
    series.prtcl = modchar(fread(id, 25, '*char'));  % Scan Protocol Name
    fseek(id, 1, 0);  % for struct byte alignment
    series.se_contrast = fread(id, 1, 'short');  % Non-zero if > 0 image used contrast(L/S)
    series.start_ras = modchar(fread(id, 1, '*char'));  % RAS letter for first scan location (L/S)
    fseek(id, 3, 0);  % for struct byte alignment
    series.start_loc = fread(id, 1, 'float');  % First scan location (L/S)
    series.end_ras = modchar(fread(id, 1, '*char'));  % RAS letter for last scan location (L/S)
    fseek(id, 3, 0);  % for struct byte alignment
    series.end_loc = fread(id, 1, 'float');  % Last scan location (L/S)
    series.se_pseq = fread(id, 1, 'short');  % Last Pulse Sequence Used (L/S)
    series.se_sortorder = fread(id, 1, 'short');  % Image Sort Order (L/S)
    series.se_lndmrkcnt = fread(id, 1, 'int');  % Landmark Counter
    series.se_nacq = fread(id, 1, 'short');  % Number of Acquisitions
    series.xbasest = fread(id, 1, 'short');  % Starting number for baselines
    series.xbaseend = fread(id, 1, 'short');  % Ending number for baselines
    series.xenhst = fread(id, 1, 'short');  % Starting number for enhanced scans
    series.xenhend = fread(id, 1, 'short');  % Ending number for enhanced scans
    fseek(id, 2, 0);  % for struct byte alignment
    series.se_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    series.se_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    fseek(id, 3, 0);  % for struct byte alignment
    series.se_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    series.se_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    series.se_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    series.se_pds_a = fread(id, 1, 'float');  % PixelData size - as stored
    series.se_pds_c = fread(id, 1, 'float');  % PixelData size - Compressed
    series.se_pds_u = fread(id, 1, 'float');  % PixelData size - UnCompressed
    series.se_checksum = fread(id, 1, 'ulong');  % Series Record checksum
    series.se_complete = fread(id, 1, 'long');  % Series Complete Flag
    series.se_numarch = fread(id, 1, 'long');  % Number of Images Archived
    series.se_imagect = fread(id, 1, 'long');  % Last Image Number Used
    series.se_numimages = fread(id, 1, 'long');  % Number of Images Existing
    series.se_images = read_vartype(id);  % kluge for VARTYPE
    series.se_numunimg = fread(id, 1, 'long');  % Number of Unstored Images
    series.se_unimages = read_vartype(id);  % kluge for VARTYPE
    series.se_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Images
    series.se_toarchive = read_vartype(id);  % kluge for VARTYPE
    series.echo1_alpha = fread(id, 1, 'float');  % Echo 1 Alpha Value
    series.echo1_beta = fread(id, 1, 'float');  % Echo 1 Beta Value
    series.echo1_window = fread(id, 1, 'ushort');  % Echo 1 Window Value
    series.echo1_level = fread(id, 1, 'short');  % Echo 1 Level Value
    series.echo2_alpha = fread(id, 1, 'float');  % Echo 2 Alpha Value
    series.echo2_beta = fread(id, 1, 'float');  % Echo 2 Beta Value
    series.echo2_window = fread(id, 1, 'ushort');  % Echo 2 Window Value
    series.echo2_level = fread(id, 1, 'short');  % Echo 2 Level Value
    series.echo3_alpha = fread(id, 1, 'float');  % Echo 3 Alpha Value
    series.echo3_beta = fread(id, 1, 'float');  % Echo 3 Beta Value
    series.echo3_window = fread(id, 1, 'ushort');  % Echo 3 Window Value
    series.echo3_level = fread(id, 1, 'short');  % Echo 3 Level Value
    series.echo4_alpha = fread(id, 1, 'float');  % Echo 4 Alpha Value
    series.echo4_beta = fread(id, 1, 'float');  % Echo 4 Beta Value
    series.echo4_window = fread(id, 1, 'ushort');  % Echo 4 Window Value
    series.echo4_level = fread(id, 1, 'short');  % Echo 4 Level Value
    series.echo5_alpha = fread(id, 1, 'float');  % Echo 5 Alpha Value
    series.echo5_beta = fread(id, 1, 'float');  % Echo 5 Beta Value
    series.echo5_window = fread(id, 1, 'ushort');  % Echo 5 Window Value
    series.echo5_level = fread(id, 1, 'short');  % Echo 5 Level Value
    series.echo6_alpha = fread(id, 1, 'float');  % Echo 6 Alpha Value
    series.echo6_beta = fread(id, 1, 'float');  % Echo 6 Beta Value
    series.echo6_window = fread(id, 1, 'ushort');  % Echo 6 Window Value
    series.echo6_level = fread(id, 1, 'short');  % Echo 6 Level Value
    series.echo7_alpha = fread(id, 1, 'float');  % Echo 7 Alpha Value
    series.echo7_beta = fread(id, 1, 'float');  % Echo 7 Beta Value
    series.echo7_window = fread(id, 1, 'ushort');  % Echo 7 Window Value
    series.echo7_level = fread(id, 1, 'short');  % Echo 7 Level Value
    series.echo8_alpha = fread(id, 1, 'float');  % Echo 8 Alpha Value
    series.echo8_beta = fread(id, 1, 'float');  % Echo 8 Beta Value
    series.echo8_window = fread(id, 1, 'ushort');  % Echo 8 Window Value
    series.echo8_level = fread(id, 1, 'short');  % Echo 8 Level Value
    series.series_uid = moduid(fread(id, 32, 'uint8'));  % Series Entity Unique ID
    series.landmark_uid = moduid(fread(id, 32, 'uint8'));  % Landmark Unique ID
    series.equipmnt_uid = moduid(fread(id, 32, 'uint8'));  % Equipment Unique ID
    series.refsopcuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    series.refsopiuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    series.schacitval = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Value
    series.schacitdesc = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Description
    series.schacitmea = modchar(fread(id, 64, '*char'));  % Sched Proc Action Item Seq - Meaning
    series.schprocstdesc = modchar(fread(id, 64, '*char'));  % Sched Proc Step Desc
    series.schprocstid = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 1
    series.reqprocstid = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 1
    series.perprocstid = modchar(fread(id, 16, '*char'));  % PPS ID
    series.perprocstdesc = modchar(fread(id, 64, '*char'));  % PPS Description
    series.reqprocstid2 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 2
    series.reqprocstid3 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 3
    series.schprocstid2 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 2
    series.schprocstid3 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 3
    series.refImgUID = moduid(fread(id, 128, 'uint8'));  % Dicom Reference Image
    series.table_entry = fread(id, 1, 'short');  % Table position for nMR and iMR
    series.SwingAngle = fread(id, 1, 'short');  % nMR - Swing Angle
    series.LateralOffset = fread(id, 1, 'short');  % nMR - Offset
    series.GradientCoil = fread(id, 1, 'short');  % Gradient Coil Selection
    series.se_subtype = fread(id, 1, 'short');  % supplements se_typ, see DICOM (0008,0008) //GSAge04506
    series.BWRT = fread(id, 1, 'short');  % for fMRI till ExptTimePts
    series.PdgmStr = modchar(fread(id, 64, '*char'));  %
    series.PdgmDesc = modchar(fread(id, 256, '*char'));  %
    series.PdgmUID = moduid(fread(id, 64, 'uint8'));  %
    series.ExpType = fread(id, 1, 'int');  %
    series.TrRest = fread(id, 1, 'int');  %
    series.TrActive = fread(id, 1, 'int');  %
    series.DumAcq = fread(id, 1, 'int');  %
    series.ApplName = modchar(fread(id, 16, '*char'));  %
    series.ApplVer = modchar(fread(id, 16, '*char'));  %
    series.ExptTimePts = fread(id, 1, 'int');  %
    series.se_padding = modchar(fread(id, 120, '*char'));  %
elseif  (ver == 11.0)
    series.se_images = read_vartype(id);  % kluge for VARTYPE
    series.se_unimages = read_vartype(id);  % kluge for VARTYPE
    series.se_toarchive = read_vartype(id);  % kluge for VARTYPE
    series.se_pds_a = fread(id, 1, 'float');  % PixelData size - as stored
    series.se_pds_c = fread(id, 1, 'float');  % PixelData size - Compressed
    series.se_pds_u = fread(id, 1, 'float');  % PixelData size - UnCompressed
    series.lmhor = fread(id, 1, 'float');  % Horizontal Landmark
    series.start_loc = fread(id, 1, 'float');  % First scan location (L/S)
    series.end_loc = fread(id, 1, 'float');  % Last scan location (L/S)
    series.echo1_alpha = fread(id, 1, 'float');  % Echo 1 Alpha Value
    series.echo1_beta = fread(id, 1, 'float');  % Echo 1 Beta Value
    series.echo2_alpha = fread(id, 1, 'float');  % Echo 2 Alpha Value
    series.echo2_beta = fread(id, 1, 'float');  % Echo 2 Beta Value
    series.echo3_alpha = fread(id, 1, 'float');  % Echo 3 Alpha Value
    series.echo3_beta = fread(id, 1, 'float');  % Echo 3 Beta Value
    series.echo4_alpha = fread(id, 1, 'float');  % Echo 4 Alpha Value
    series.echo4_beta = fread(id, 1, 'float');  % Echo 4 Beta Value
    series.echo5_alpha = fread(id, 1, 'float');  % Echo 5 Alpha Value
    series.echo5_beta = fread(id, 1, 'float');  % Echo 5 Beta Value
    series.echo6_alpha = fread(id, 1, 'float');  % Echo 6 Alpha Value
    series.echo6_beta = fread(id, 1, 'float');  % Echo 6 Beta Value
    series.echo7_alpha = fread(id, 1, 'float');  % Echo 7 Alpha Value
    series.echo7_beta = fread(id, 1, 'float');  % Echo 7 Beta Value
    series.echo8_alpha = fread(id, 1, 'float');  % Echo 8 Alpha Value
    series.echo8_beta = fread(id, 1, 'float');  % Echo 8 Beta Value
    series.se_checksum = fread(id, 1, 'ulong');  % Series Record checksum
    series.se_complete = fread(id, 1, 'long');  % Series Complete Flag
    series.se_numarch = fread(id, 1, 'long');  % Number of Images Archived
    series.se_imagect = fread(id, 1, 'long');  % Last Image Number Used
    series.se_numimages = fread(id, 1, 'long');  % Number of Images Existing
    series.se_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    series.se_numunimg = fread(id, 1, 'long');  % Number of Unstored Images
    series.se_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Images
    series.se_datetime = fread(id, 1, 'int');  % Allocation Series Data/Time stamp
    series.se_actual_dt = fread(id, 1, 'int');  % Actual Series Data/Time stamp
    series.position = fread(id, 1, 'int');  % Patient Position
    series.entry = fread(id, 1, 'int');  % Patient Entry
    series.se_lndmrkcnt = fread(id, 1, 'int');  % Landmark Counter
    series.se_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    series.ExpType = fread(id, 1, 'int');  %
    series.TrRest = fread(id, 1, 'int');  %
    series.TrActive = fread(id, 1, 'int');  %
    series.DumAcq = fread(id, 1, 'int');  %
    series.ExptTimePts = fread(id, 1, 'int');  %
    series.se_exno = fread(id, 1, 'ushort');  % Exam Number
    series.echo1_window = fread(id, 1, 'ushort');  % Echo 1 Window Value
    series.echo2_window = fread(id, 1, 'ushort');  % Echo 2 Window Value
    series.echo3_window = fread(id, 1, 'ushort');  % Echo 3 Window Value
    series.echo4_window = fread(id, 1, 'ushort');  % Echo 4 Window Value
    series.echo5_window = fread(id, 1, 'ushort');  % Echo 5 Window Value
    series.echo6_window = fread(id, 1, 'ushort');  % Echo 6 Window Value
    series.echo7_window = fread(id, 1, 'ushort');  % Echo 7 Window Value
    series.echo8_window = fread(id, 1, 'ushort');  % Echo 8 Window Value
    series.echo8_level = fread(id, 1, 'short');  % Echo 8 Level Value
    series.echo7_level = fread(id, 1, 'short');  % Echo 7 Level Value
    series.echo6_level = fread(id, 1, 'short');  % Echo 6 Level Value
    series.echo5_level = fread(id, 1, 'short');  % Echo 5 Level Value
    series.echo4_level = fread(id, 1, 'short');  % Echo 4 Level Value
    series.echo3_level = fread(id, 1, 'short');  % Echo 3 Level Value
    series.echo2_level = fread(id, 1, 'short');  % Echo 2 Level Value
    series.echo1_level = fread(id, 1, 'short');  % Echo 1 Level Value
    series.se_no = fread(id, 1, 'short');  % Series Number
    series.se_typ = fread(id, 1, 'short');  % Series Type
    series.se_source = fread(id, 1, 'short');  % Series from which prescribed
    series.se_plane = fread(id, 1, 'short');  % Most-like Plane (for L/S)
    series.scan_type = fread(id, 1, 'short');  % Scout or Axial (for CT)
    series.se_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    series.se_contrast = fread(id, 1, 'short');  % Non-zero if > 0 image used contrast(L/S)
    series.se_pseq = fread(id, 1, 'short');  % Last Pulse Sequence Used (L/S)
    series.se_sortorder = fread(id, 1, 'short');  % Image Sort Order (L/S)
    series.se_nacq = fread(id, 1, 'short');  % Number of Acquisitions
    series.xbasest = fread(id, 1, 'short');  % Starting number for baselines
    series.xbaseend = fread(id, 1, 'short');  % Ending number for baselines
    series.xenhst = fread(id, 1, 'short');  % Starting number for enhanced scans
    series.xenhend = fread(id, 1, 'short');  % Ending number for enhanced scans
    series.table_entry = fread(id, 1, 'short');  % Table position for nMR and iMR
    series.SwingAngle = fread(id, 1, 'short');  % nMR - Swing Angle
    series.LateralOffset = fread(id, 1, 'short');  % nMR - Offset
    series.GradientCoil = fread(id, 1, 'short');  % Gradient Coil Selection
    series.se_subtype = fread(id, 1, 'short');  % supplements se_typ, see DICOM (0008,0008) //GSAge04506
    series.BWRT = fread(id, 1, 'short');  % for fMRI till ExptTimePts
    series.assetcal_serno = fread(id, 1, 'short');  % Calibration Series number
    series.assetcal_scnno = fread(id, 1, 'short');  % Calibration Scan number
    series.content_qualifn = fread(id, 1, 'short');  % PRODUCT/RESEARCH/SERVICE
    series.purecal_serno = fread(id, 1, 'short');  % Calibration Series number
    series.purecal_scnno = fread(id, 1, 'short');  % Calibration Scan number
    series.short_padding = fread(id, 2, 'short');  %
    series.se_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    series.se_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    series.se_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Series
    series.se_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    series.se_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Series
    series.se_desc = modchar(fread(id, 65, '*char'));  % Series Description
    series.pr_sysid = modchar(fread(id, 9, '*char'));  % Primary Receiver Suite and Host
    series.pansysid = modchar(fread(id, 9, '*char'));  % Archiver Suite and Host
    series.anref = modchar(fread(id, 3, '*char'));  % Anatomical reference
    series.prtcl = modchar(fread(id, 25, '*char'));  % Scan Protocol Name
    series.start_ras = modchar(fread(id, 1, '*char'));  % RAS letter for first scan location (L/S)
    series.end_ras = modchar(fread(id, 1, '*char'));  % RAS letter for last scan location (L/S)
    series.series_uid = moduid(fread(id, 32, 'uint8'));  % Series Entity Unique ID
    series.landmark_uid = moduid(fread(id, 32, 'uint8'));  % Landmark Unique ID
    series.equipmnt_uid = moduid(fread(id, 32, 'uint8'));  % Equipment Unique ID
    series.refsopcuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    series.refsopiuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    series.schacitval = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Value
    series.schacitdesc = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Description
    series.schacitmea = modchar(fread(id, 64, '*char'));  % Sched Proc Action Item Seq - Meaning
    series.schprocstdesc = modchar(fread(id, 65, '*char'));  % Sched Proc Step Desc
    series.schprocstid = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 1
    series.reqprocstid = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 1
    series.perprocstid = modchar(fread(id, 16, '*char'));  % PPS ID
    series.perprocstdesc = modchar(fread(id, 65, '*char'));  % PPS Description
    series.reqprocstid2 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 2
    series.reqprocstid3 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 3
    series.schprocstid2 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 2
    series.schprocstid3 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 3
    series.refImgUID = moduid(fread(id, 128, 'uint8'));  % Dicom Reference Image
    series.PdgmStr = modchar(fread(id, 64, '*char'));  %
    series.PdgmDesc = modchar(fread(id, 256, '*char'));  %
    series.PdgmUID = moduid(fread(id, 64, 'uint8'));  %
    series.ApplName = modchar(fread(id, 16, '*char'));  %
    series.ApplVer = modchar(fread(id, 16, '*char'));  %
    series.asset_appl = modchar(fread(id, 12, '*char'));  % Asset application name
    series.scic_a = modchar(fread(id, 32, '*char'));  % Scic_a values from CoilConfig.cfg
    series.scic_s = modchar(fread(id, 32, '*char'));  % Scic_s values from CoilConfig.cfg
    series.scic_c = modchar(fread(id, 32, '*char'));  % Scic_c values from CoilConfig.cfg
    series.pure_cfg_params = modchar(fread(id, 64, '*char'));  % PURE Config Parameters from pure.cfg
    series.se_padding = modchar(fread(id, 423, '*char'));  % Spare Space

elseif (ver >= 14.0 & ver < 20) % 14.0 & 14.3
    % series.double_padding = modchar(fread(id, 8, 'char')); % double padding
    series.double_padding = modchar(fread(id, 56, 'char')); % by jz
    
    series.se_pds_a = fread(id, 1, 'float');  % PixelData size - as stored
    series.se_pds_c = fread(id, 1, 'float');  % PixelData size - Compressed
    series.se_pds_u = fread(id, 1, 'float');  % PixelData size - UnCompressed
    series.lmhor = fread(id, 1, 'float');  % Horizontal Landmark
    series.start_loc = fread(id, 1, 'float');  % First scan location (L/S)
    series.end_loc = fread(id, 1, 'float');  % Last scan location (L/S)
    series.echo1_alpha = fread(id, 1, 'float');  % Echo 1 Alpha Value
    series.echo1_beta = fread(id, 1, 'float');  % Echo 1 Beta Value
    series.echo2_alpha = fread(id, 1, 'float');  % Echo 2 Alpha Value
    series.echo2_beta = fread(id, 1, 'float');  % Echo 2 Beta Value
    series.echo3_alpha = fread(id, 1, 'float');  % Echo 3 Alpha Value
    series.echo3_beta = fread(id, 1, 'float');  % Echo 3 Beta Value
    series.echo4_alpha = fread(id, 1, 'float');  % Echo 4 Alpha Value
    series.echo4_beta = fread(id, 1, 'float');  % Echo 4 Beta Value
    series.echo5_alpha = fread(id, 1, 'float');  % Echo 5 Alpha Value
    series.echo5_beta = fread(id, 1, 'float');  % Echo 5 Beta Value
    series.echo6_alpha = fread(id, 1, 'float');  % Echo 6 Alpha Value
    series.echo6_beta = fread(id, 1, 'float');  % Echo 6 Beta Value
    series.echo7_alpha = fread(id, 1, 'float');  % Echo 7 Alpha Value
    series.echo7_beta = fread(id, 1, 'float');  % Echo 7 Beta Value
    series.echo8_alpha = fread(id, 1, 'float');  % Echo 8 Alpha Value
    series.echo8_beta = fread(id, 1, 'float');  % Echo 8 Beta Value    
    series.float_padding = modchar(fread(id, 32, 'char')); % float padding
    
    series.checksum = fread(id, 1, 'long');  % Checksum 
    series.se_complete = fread(id, 1, 'long');  % Series Complete Flag
    series.se_numarch = fread(id, 1, 'long');  % Number of Images Archived
    series.se_imagect = fread(id, 1, 'long');  % Last Image Number Used
    series.se_numimages = fread(id, 1, 'long');  % Number of Images Existing
    series.se_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    series.se_numunimg = fread(id, 1, 'long');  % Number of Unstored Images
    series.se_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Images    
    series.long_padding = modchar(fread(id, 32, 'char')); % long padding
    
    series.se_datetime = fread(id, 1, 'int');  % Allocation Series Data/Time stamp
    series.se_actual_dt = fread(id, 1, 'int');  % Actual Series Data/Time stamp
    series.position = fread(id, 1, 'int');  % Patient Position
    series.entry = fread(id, 1, 'int');  % Patient Entry
    series.se_lndmrkcnt = fread(id, 1, 'int');  % Landmark Counter
    series.se_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    series.ExpType = fread(id, 1, 'int');  %
    series.TrRest = fread(id, 1, 'int');  %
    series.TrActive = fread(id, 1, 'int');  %
    series.DumAcq = fread(id, 1, 'int');  %
    series.ExptTimePts = fread(id, 1, 'int');  %    
    series.int_padding = modchar(fread(id, 64, 'char')); % int padding
    
    series.se_exno = fread(id, 1, 'ushort');  % Exam Number
    series.echo1_window = fread(id, 1, 'ushort');  % Echo 1 Window Value
    series.echo2_window = fread(id, 1, 'ushort');  % Echo 2 Window Value
    series.echo3_window = fread(id, 1, 'ushort');  % Echo 3 Window Value
    series.echo4_window = fread(id, 1, 'ushort');  % Echo 4 Window Value
    series.echo5_window = fread(id, 1, 'ushort');  % Echo 5 Window Value
    series.echo6_window = fread(id, 1, 'ushort');  % Echo 6 Window Value
    series.echo7_window = fread(id, 1, 'ushort');  % Echo 7 Window Value
    series.echo8_window = fread(id, 1, 'ushort');  % Echo 8 Window Value
    series.echo8_level = fread(id, 1, 'short');  % Echo 8 Level Value
    series.echo7_level = fread(id, 1, 'short');  % Echo 7 Level Value
    series.echo6_level = fread(id, 1, 'short');  % Echo 6 Level Value
    series.echo5_level = fread(id, 1, 'short');  % Echo 5 Level Value
    series.echo4_level = fread(id, 1, 'short');  % Echo 4 Level Value
    series.echo3_level = fread(id, 1, 'short');  % Echo 3 Level Value
    series.echo2_level = fread(id, 1, 'short');  % Echo 2 Level Value
    series.echo1_level = fread(id, 1, 'short');  % Echo 1 Level Value
    series.se_no = fread(id, 1, 'short');  % Series Number
    series.se_typ = fread(id, 1, 'short');  % Series Type
    series.se_source = fread(id, 1, 'short');  % Series from which prescribed
    series.se_plane = fread(id, 1, 'short');  % Most-like Plane (for L/S)
    series.scan_type = fread(id, 1, 'short');  % Scout or Axial (for CT)
    series.se_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    series.se_contrast = fread(id, 1, 'short');  % Non-zero if > 0 image used contrast(L/S)
    series.se_pseq = fread(id, 1, 'short');  % Last Pulse Sequence Used (L/S)
    series.se_sortorder = fread(id, 1, 'short');  % Image Sort Order (L/S)
    series.se_nacq = fread(id, 1, 'short');  % Number of Acquisitions
    series.xbasest = fread(id, 1, 'short');  % Starting number for baselines
    series.xbaseend = fread(id, 1, 'short');  % Ending number for baselines
    series.xenhst = fread(id, 1, 'short');  % Starting number for enhanced scans
    series.xenhend = fread(id, 1, 'short');  % Ending number for enhanced scans
    series.table_entry = fread(id, 1, 'short');  % Table position for nMR and iMR
    series.SwingAngle = fread(id, 1, 'short');  % nMR - Swing Angle
    series.LateralOffset = fread(id, 1, 'short');  % nMR - Offset
    series.GradientCoil = fread(id, 1, 'short');  % Gradient Coil Selection
    series.se_subtype = fread(id, 1, 'short');  % supplements se_typ, see DICOM (0008,0008) //GSAge04506
    series.BWRT = fread(id, 1, 'short');  % for fMRI till ExptTimePts
    series.assetcal_serno = fread(id, 1, 'short');  % Calibration Series number
    series.assetcal_scnno = fread(id, 1, 'short');  % Calibration Scan number
    series.content_qualifn = fread(id, 1, 'short');  % PRODUCT/RESEARCH/SERVICE
    series.purecal_serno = fread(id, 1, 'short');  % Calibration Series number
    series.purecal_scnno = fread(id, 1, 'short');  % Calibration Scan number
    series.short_padding = modchar(fread(id, 52, 'char'));  % Short padding
    
    series.se_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    series.se_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    series.se_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Series
    series.se_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    series.se_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Series
    series.se_desc = modchar(fread(id, 65, '*char'));  % Series Description
    series.pr_sysid = modchar(fread(id, 9, '*char'));  % Primary Receiver Suite and Host
    series.pansysid = modchar(fread(id, 9, '*char'));  % Archiver Suite and Host
    series.anref = modchar(fread(id, 3, '*char'));  % Anatomical reference
    series.prtcl = modchar(fread(id, 25, '*char'));  % Scan Protocol Name
    series.start_ras = modchar(fread(id, 1, '*char'));  % RAS letter for first scan location (L/S)
    series.end_ras = modchar(fread(id, 1, '*char'));  % RAS letter for last scan location (L/S)
    series.series_uid = moduid(fread(id, 32, 'uint8'));  % Series Entity Unique ID
    series.landmark_uid = moduid(fread(id, 32, 'uint8'));  % Landmark Unique ID
    series.equipmnt_uid = moduid(fread(id, 32, 'uint8'));  % Equipment Unique ID
    series.refsopcuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    series.refsopiuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    series.schacitval = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Value
    series.schacitdesc = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Description
    series.schacitmea = modchar(fread(id, 64, '*char'));  % Sched Proc Action Item Seq - Meaning
    series.schprocstdesc = modchar(fread(id, 65, '*char'));  % Sched Proc Step Desc
    series.schprocstid = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 1
    series.reqprocstid = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 1
    series.perprocstid = modchar(fread(id, 16, '*char'));  % PPS ID
    series.perprocstdesc = modchar(fread(id, 65, '*char'));  % PPS Description
    series.reqprocstid2 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 2
    series.reqprocstid3 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 3
    series.schprocstid2 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 2
    series.schprocstid3 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 3
    series.refImgUID0 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image0
    series.refImgUID1 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image1
    series.refImgUID2 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image2
    series.refImgUID3 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image3   
    series.PdgmStr = modchar(fread(id, 64, '*char'));  %
    series.PdgmDesc = modchar(fread(id, 256, '*char'));  %
    series.PdgmUID = moduid(fread(id, 64, 'uint8'));  %
    series.ApplName = modchar(fread(id, 16, '*char'));  %
    series.ApplVer = modchar(fread(id, 16, '*char'));  %
    series.asset_appl = modchar(fread(id, 12, '*char'));  % Asset application name
    series.scic_a = modchar(fread(id, 32, '*char'));  % Scic_a values from CoilConfig.cfg
    series.scic_s = modchar(fread(id, 32, '*char'));  % Scic_s values from CoilConfig.cfg
    series.scic_c = modchar(fread(id, 32, '*char'));  % Scic_c values from CoilConfig.cfg
    series.pure_cfg_params = modchar(fread(id, 64, '*char'));  % PURE Config Parameters from pure.cfg    
    series.se_padding = modchar(fread(id, 215, '*char'));  % Spare Space
    
elseif (ver >= 20)   % 20.0 
    series.double_padding = modchar(fread(id, 256, 'char')); % double padding
    
    series.se_pds_a = fread(id, 1, 'float');  % PixelData size - as stored
    series.se_pds_c = fread(id, 1, 'float');  % PixelData size - Compressed
    series.se_pds_u = fread(id, 1, 'float');  % PixelData size - UnCompressed
    series.lmhor = fread(id, 1, 'float');  % Horizontal Landmark
    series.start_loc = fread(id, 1, 'float');  % First scan location (L/S)
    series.end_loc = fread(id, 1, 'float');  % Last scan location (L/S)
    series.echo1_alpha = fread(id, 1, 'float');  % Echo 1 Alpha Value
    series.echo1_beta = fread(id, 1, 'float');  % Echo 1 Beta Value
    series.echo2_alpha = fread(id, 1, 'float');  % Echo 2 Alpha Value
    series.echo2_beta = fread(id, 1, 'float');  % Echo 2 Beta Value
    series.echo3_alpha = fread(id, 1, 'float');  % Echo 3 Alpha Value
    series.echo3_beta = fread(id, 1, 'float');  % Echo 3 Beta Value
    series.echo4_alpha = fread(id, 1, 'float');  % Echo 4 Alpha Value
    series.echo4_beta = fread(id, 1, 'float');  % Echo 4 Beta Value
    series.echo5_alpha = fread(id, 1, 'float');  % Echo 5 Alpha Value
    series.echo5_beta = fread(id, 1, 'float');  % Echo 5 Beta Value
    series.echo6_alpha = fread(id, 1, 'float');  % Echo 6 Alpha Value
    series.echo6_beta = fread(id, 1, 'float');  % Echo 6 Beta Value
    series.echo7_alpha = fread(id, 1, 'float');  % Echo 7 Alpha Value
    series.echo7_beta = fread(id, 1, 'float');  % Echo 7 Beta Value
    series.echo8_alpha = fread(id, 1, 'float');  % Echo 8 Alpha Value
    series.echo8_beta = fread(id, 1, 'float');  % Echo 8 Beta Value    
    series.float_padding = modchar(fread(id, 128, 'char')); % float padding
    
    series.se_complete = fread(id, 1, 'long');  % Series Complete Flag
    series.se_numarch = fread(id, 1, 'long');  % Number of Images Archived
    series.se_imagect = fread(id, 1, 'long');  % Last Image Number Used
    series.se_numimages = fread(id, 1, 'long');  % Number of Images Existing
    series.se_delta_cnt = fread(id, 1, 'long');  % Indicates number of updates to header
    series.se_numunimg = fread(id, 1, 'long');  % Number of Unstored Images
    series.se_toarchcnt = fread(id, 1, 'long');  % Number of Unarchived Images    
    series.int_padding1 = modchar(fread(id, 132, 'char')); % int padding1
    
    series.se_datetime = fread(id, 1, 'int');  % Allocation Series Data/Time stamp
    series.se_actual_dt = fread(id, 1, 'int');  % Actual Series Data/Time stamp
    series.position = fread(id, 1, 'int');  % Patient Position
    series.entry = fread(id, 1, 'int');  % Patient Entry
    series.se_lndmrkcnt = fread(id, 1, 'int');  % Landmark Counter
    series.se_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    series.ExpType = fread(id, 1, 'int');  %
    series.TrRest = fread(id, 1, 'int');  %
    series.TrActive = fread(id, 1, 'int');  %
    series.DumAcq = fread(id, 1, 'int');  %
    series.ExptTimePts = fread(id, 1, 'int');  %    
    series.int_padding2 = modchar(fread(id, 132, 'char')); % int padding2
    
    series.se_exno = fread(id, 1, 'ushort');  % Exam Number
    series.echo1_window = fread(id, 1, 'ushort');  % Echo 1 Window Value
    series.echo2_window = fread(id, 1, 'ushort');  % Echo 2 Window Value
    series.echo3_window = fread(id, 1, 'ushort');  % Echo 3 Window Value
    series.echo4_window = fread(id, 1, 'ushort');  % Echo 4 Window Value
    series.echo5_window = fread(id, 1, 'ushort');  % Echo 5 Window Value
    series.echo6_window = fread(id, 1, 'ushort');  % Echo 6 Window Value
    series.echo7_window = fread(id, 1, 'ushort');  % Echo 7 Window Value
    series.echo8_window = fread(id, 1, 'ushort');  % Echo 8 Window Value
    series.echo8_level = fread(id, 1, 'short');  % Echo 8 Level Value
    series.echo7_level = fread(id, 1, 'short');  % Echo 7 Level Value
    series.echo6_level = fread(id, 1, 'short');  % Echo 6 Level Value
    series.echo5_level = fread(id, 1, 'short');  % Echo 5 Level Value
    series.echo4_level = fread(id, 1, 'short');  % Echo 4 Level Value
    series.echo3_level = fread(id, 1, 'short');  % Echo 3 Level Value
    series.echo2_level = fread(id, 1, 'short');  % Echo 2 Level Value
    series.echo1_level = fread(id, 1, 'short');  % Echo 1 Level Value
    series.se_no = fread(id, 1, 'short');  % Series Number
    series.se_typ = fread(id, 1, 'short');  % Series Type
    series.se_source = fread(id, 1, 'short');  % Series from which prescribed
    series.se_plane = fread(id, 1, 'short');  % Most-like Plane (for L/S)
    series.scan_type = fread(id, 1, 'short');  % Scout or Axial (for CT)
    series.se_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    series.se_contrast = fread(id, 1, 'short');  % Non-zero if > 0 image used contrast(L/S)
    series.se_pseq = fread(id, 1, 'short');  % Last Pulse Sequence Used (L/S)
    series.se_sortorder = fread(id, 1, 'short');  % Image Sort Order (L/S)
    series.se_nacq = fread(id, 1, 'short');  % Number of Acquisitions
    series.xbasest = fread(id, 1, 'short');  % Starting number for baselines
    series.xbaseend = fread(id, 1, 'short');  % Ending number for baselines
    series.xenhst = fread(id, 1, 'short');  % Starting number for enhanced scans
    series.xenhend = fread(id, 1, 'short');  % Ending number for enhanced scans
    series.table_entry = fread(id, 1, 'short');  % Table position for nMR and iMR
    series.SwingAngle = fread(id, 1, 'short');  % nMR - Swing Angle
    series.LateralOffset = fread(id, 1, 'short');  % nMR - Offset
    series.GradientCoil = fread(id, 1, 'short');  % Gradient Coil Selection
    series.se_subtype = fread(id, 1, 'short');  % supplements se_typ, see DICOM (0008,0008) //GSAge04506
    series.BWRT = fread(id, 1, 'short');  % for fMRI till ExptTimePts
    series.assetcal_serno = fread(id, 1, 'short');  % Calibration Series number
    series.assetcal_scnno = fread(id, 1, 'short');  % Calibration Scan number
    series.content_qualifn = fread(id, 1, 'short');  % PRODUCT/RESEARCH/SERVICE
    series.purecal_serno = fread(id, 1, 'short');  % Calibration Series number
    series.purecal_scnno = fread(id, 1, 'short');  % Calibration Scan number
    series.ideal = fread(id, 1, 'short');  % Ideal    
    series.short_padding = modchar(fread(id, 66, 'char'));  % Short padding
    
    series.se_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    series.se_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    series.se_suid = modchar(fread(id, 4, '*char'));  % Suite ID for this Series
    series.se_alloc_key = modchar(fread(id, 13, '*char'));  % Process that allocated this record
    series.se_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Series
    series.se_desc = modchar(fread(id, 65, '*char'));  % Series Description
    series.pr_sysid = modchar(fread(id, 9, '*char'));  % Primary Receiver Suite and Host
    series.pansysid = modchar(fread(id, 9, '*char'));  % Archiver Suite and Host
    series.anref = modchar(fread(id, 3, '*char'));  % Anatomical reference
    series.prtcl = modchar(fread(id, 25, '*char'));  % Scan Protocol Name
    series.start_ras = modchar(fread(id, 1, '*char'));  % RAS letter for first scan location (L/S)
    series.end_ras = modchar(fread(id, 1, '*char'));  % RAS letter for last scan location (L/S)
    series.series_uid = moduid(fread(id, 32, 'uint8'));  % Series Entity Unique ID
    series.landmark_uid = moduid(fread(id, 32, 'uint8'));  % Landmark Unique ID
    series.equipmnt_uid = moduid(fread(id, 32, 'uint8'));  % Equipment Unique ID
    series.refsopcuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Class UID
    series.refsopiuids = moduid(fread(id, 32, 'uint8'));  % Ref SOP Instance UID
    series.schacitval = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Value
    series.schacitdesc = modchar(fread(id, 16, '*char'));  % Sched Proc Action Item Seq - Description
    series.schacitmea = modchar(fread(id, 64, '*char'));  % Sched Proc Action Item Seq - Meaning
    series.schprocstdesc = modchar(fread(id, 65, '*char'));  % Sched Proc Step Desc
    series.schprocstid = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 1
    series.reqprocstid = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 1
    series.perprocstid = modchar(fread(id, 16, '*char'));  % PPS ID
    series.perprocstdesc = modchar(fread(id, 65, '*char'));  % PPS Description
    series.reqprocstid2 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 2
    series.reqprocstid3 = modchar(fread(id, 16, '*char'));  % Req Proc Step ID 3
    series.schprocstid2 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 2
    series.schprocstid3 = modchar(fread(id, 16, '*char'));  % Sched Proc Step ID 3
    series.refImgUID0 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image0
    series.refImgUID1 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image1
    series.refImgUID2 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image2
    series.refImgUID3 = moduid(fread(id, 32, 'uint8'));  % Dicom Reference Image3 
    series.PdgmStr = modchar(fread(id, 64, '*char'));  %
    series.PdgmDesc = modchar(fread(id, 256, '*char'));  %
    series.PdgmUID = moduid(fread(id, 64, 'uint8'));  %
    series.ApplName = modchar(fread(id, 16, '*char'));  %
    series.ApplVer = modchar(fread(id, 16, '*char'));  %
    series.asset_appl = modchar(fread(id, 12, '*char'));  % Asset application name
    series.scic_a = modchar(fread(id, 32, '*char'));  % Scic_a values from CoilConfig.cfg
    series.scic_s = modchar(fread(id, 32, '*char'));  % Scic_s values from CoilConfig.cfg
    series.scic_c = modchar(fread(id, 32, '*char'));  % Scic_c values from CoilConfig.cfg
    series.pure_cfg_params = modchar(fread(id, 64, '*char'));  % PURE Config Parameters from pure.cfg    
    series.se_padding = modchar(fread(id, 251, '*char'));  % Spare Space
end

end

function image = read_image(id, ver)

pos_image = ftell(id);

if (ver < 9.0)
    image.im_suid = modchar(fread(id, 4, '*char'));  % Suite id for this image
    image.im_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    image.im_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Image
    fseek(id, 1, 0);  % for struct byte alignment
    image.im_exno = fread(id, 1, 'ushort');  % Exam number for this image
    image.im_seno = fread(id, 1, 'short');  % Series Number for this image
    image.im_no = fread(id, 1, 'short');  % Image Number
    fseek(id, 2, 0);  % for struct byte alignment
    image.im_datetime = fread(id, 1, 'int');  % Allocation Image date/time stamp
    image.im_actual_dt = fread(id, 1, 'int');  % Actual Image date/time stamp
    image.sctime = fread(id, 1, 'float');  % Duration of scan
    image.slthick = fread(id, 1, 'float');  % Slice Thickness (mm)
    image.imatrix_X = fread(id, 1, 'short');  % Image matrix size - X
    image.imatrix_Y = fread(id, 1, 'short');  % Image matrix size - Y
    image.dfov = fread(id, 1, 'float');  % Display field of view - X (mm)
    image.dfov_rect = fread(id, 1, 'float');  % Display field of view - Y (if different)
    image.dim_X = fread(id, 1, 'float');  % Image dimension - X
    image.dim_Y = fread(id, 1, 'float');  % Image dimension - Y
    image.pixsize_X = fread(id, 1, 'float');  % Image pixel size - X
    image.pixsize_Y = fread(id, 1, 'float');  % Image pixel size - Y
    image.pdid = modchar(fread(id, 14, '*char'));  % Pixel Data ID
    image.contrastIV = modchar(fread(id, 17, '*char'));  % IV Contrast Agent
    image.contrastOral = modchar(fread(id, 17, '*char'));  % Oral Contrast Agent
    image.contmode = fread(id, 1, 'short');  % Image Contrast Mode
    image.serrx = fread(id, 1, 'short');  % Series from which prescribed
    image.imgrx = fread(id, 1, 'short');  % Image from which prescribed
    image.screenformat = fread(id, 1, 'short');  % Screen Format(8/16 bit)
    image.plane = fread(id, 1, 'short');  % Plane Type
    fseek(id, 2, 0);  % for struct byte alignment
    image.scanspacing = fread(id, 1, 'float');  % Spacing between scans (mm?)
    image.im_compress = fread(id, 1, 'short');  % Image compression type for allocation
    image.im_scouttype = fread(id, 1, 'short');  % Scout Type (AP or lateral)
    image.loc_ras = modchar(fread(id, 1, '*char'));  % RAS letter of image location
    fseek(id, 3, 0);  % for struct byte alignment
    image.loc = fread(id, 1, 'float');  % Image location
    image.ctr_R = fread(id, 1, 'float');  % Center R coord of plane image
    image.ctr_A = fread(id, 1, 'float');  % Center A coord of plane image
    image.ctr_S = fread(id, 1, 'float');  % Center S coord of plane image
    image.norm_R = fread(id, 1, 'float');  % Normal R coord
    image.norm_A = fread(id, 1, 'float');  % Normal A coord
    image.norm_S = fread(id, 1, 'float');  % Normal S coord
    image.tlhc_R = fread(id, 1, 'float');  % R Coord of Top Left Hand Corner
    image.tlhc_A = fread(id, 1, 'float');  % A Coord of Top Left Hand Corner
    image.tlhc_S = fread(id, 1, 'float');  % S Coord of Top Left Hand Corner
    image.trhc_R = fread(id, 1, 'float');  % R Coord of Top Right Hand Corner
    image.trhc_A = fread(id, 1, 'float');  % A Coord of Top Right Hand Corner
    image.trhc_S = fread(id, 1, 'float');  % S Coord of Top Right Hand Corner
    image.brhc_R = fread(id, 1, 'float');  % R Coord of Bottom Right Hand Corner
    image.brhc_A = fread(id, 1, 'float');  % A Coord of Bottom Right Hand Corner
    image.brhc_S = fread(id, 1, 'float');  % S Coord of Bottom Right Hand Corner
    image.forimgrev = modchar(fread(id, 4, '*char'));  % Foreign Image Revision
    image.tr = fread(id, 1, 'int');  % Pulse repetition time(usec)
    image.ti = fread(id, 1, 'int');  % Pulse inversion time(usec)
    image.te = fread(id, 1, 'int');  % Pulse echo time(usec)
    image.te2 = fread(id, 1, 'int');  % Second echo echo (usec)
    image.numecho = fread(id, 1, 'short');  % Number of echoes
    image.echonum = fread(id, 1, 'short');  % Echo Number
    image.tbldlta = fread(id, 1, 'float');  % Table Delta
    image.nex = fread(id, 1, 'float');  % Number of Excitations
    image.contig = fread(id, 1, 'short');  % Continuous Slices Flag
    image.hrtrate = fread(id, 1, 'short');  % Cardiac Heart Rate (bpm)
    image.tdel = fread(id, 1, 'int');  % Delay time after trigger (msec)
    image.saravg = fread(id, 1, 'float');  % Average SAR
    image.sarpeak = fread(id, 1, 'float');  % Peak SAR
    image.monsar = fread(id, 1, 'short');  % Monitor SAR flag
    image.trgwindow = fread(id, 1, 'short');  % Trigger window (% of R-R interval)
    image.reptime = fread(id, 1, 'float');  % Cardiac repetition time
    image.imgpcyc = fread(id, 1, 'short');  % Images per cardiac cycle
    image.xmtgain = fread(id, 1, 'short');  % Actual Transmit Gain (.1 db)
    image.rcvgain1 = fread(id, 1, 'short');  % Actual Receive Gain Analog (.1 db)
    image.rcvgain2 = fread(id, 1, 'short');  % Actual Receive Gain Digital (.1 db)
    image.mr_flip = fread(id, 1, 'short');  % Flip Angle for GRASS scans (deg.)
    fseek(id, 2, 0);  % for struct byte alignment
    image.mindat = fread(id, 1, 'int');  % Minimum Delay after Trigger (uSec)
    image.cphase = fread(id, 1, 'short');  % Total Cardiac Phase prescribed
    image.swappf = fread(id, 1, 'short');  % Swap Phase/Frequency Axis
    image.pauseint = fread(id, 1, 'short');  % Pause Interval (slices)
    fseek(id, 2, 0);  % for struct byte alignment
    image.pausetime = fread(id, 1, 'float');  % Pause Time
    image.obplane = fread(id, 1, 'int');  % Oblique Plane
    image.slocfov = fread(id, 1, 'int');  % Slice Offsets on Freq axis
    image.xmtfreq = fread(id, 1, 'int');  % Center Frequency (0.1 Hz)
    image.autoxmtfreq = fread(id, 1, 'int');  % Auto Center Frequency (0.1 Hz)
    image.autoxmtgain = fread(id, 1, 'short');  % Auto Transmit Gain (0.1 dB)
    image.prescan_r1 = fread(id, 1, 'short');  % PreScan R1 - Analog
    image.prescan_r2 = fread(id, 1, 'short');  % PreScan R2 - Digital
    fseek(id, 2, 0);  % for struct byte alignment
    image.user_bitmap = fread(id, 1, 'int');  % Bitmap defining user CVs
    image.cenfreq = fread(id, 1, 'short');  % Center Frequency Method
    image.imode = fread(id, 1, 'short');  % Imaging Mode
    image.iopt = fread(id, 1, 'int');  % Imaging Options
    image.pseq = fread(id, 1, 'short');  % Pulse Sequence
    image.pseqmode = fread(id, 1, 'short');  % Pulse Sequence Mode
    image.psdname = modchar(fread(id, 33, '*char'));  % Pulse Sequence Name
    fseek(id, 3, 0);  % for struct byte alignment
    image.psd_datetime = fread(id, 1, 'int');  % PSD Creation Date and Time
    image.psd_iname = modchar(fread(id, 13, '*char'));  % PSD name from inside PSD
    fseek(id, 1, 0);  % for struct byte alignment
    image.ctyp = fread(id, 1, 'short');  % Coil Type
    image.cname = modchar(fread(id, 17, '*char'));  % Coil Name
    fseek(id, 1, 0);  % for struct byte alignment
    image.surfctyp = fread(id, 1, 'short');  % Surface Coil Type
    image.surfcext = fread(id, 1, 'short');  % Extremity Coil Flag
    fseek(id, 2, 0);  % for struct byte alignment
    image.rawrunnum = fread(id, 1, 'int');  % RawData Run Number
    image.cal_fldstr = fread(id, 1, 'ulong');  % Calibrated Field Strength (x10 uGauss)
    image.supp_tech = fread(id, 1, 'short');  % SAT fat/water/none
    fseek(id, 2, 0);  % for struct byte alignment
    image.vbw = fread(id, 1, 'float');  % Variable Bandwidth (Hz)
    image.slquant = fread(id, 1, 'short');  % Number of slices in this scan group
    image.gpre = fread(id, 1, 'short');  % Graphically prescribed
    image.intr_del = fread(id, 1, 'int');  % Interimage/interloc delay (uSec)
    image.user0 = fread(id, 1, 'float');  % User Variable 0
    image.user1 = fread(id, 1, 'float');  % User Variable 1
    image.user2 = fread(id, 1, 'float');  % User Variable 2
    image.user3 = fread(id, 1, 'float');  % User Variable 3
    image.user4 = fread(id, 1, 'float');  % User Variable 4
    image.user5 = fread(id, 1, 'float');  % User Variable 5
    image.user6 = fread(id, 1, 'float');  % User Variable 6
    image.user7 = fread(id, 1, 'float');  % User Variable 7
    image.user8 = fread(id, 1, 'float');  % User Variable 8
    image.user9 = fread(id, 1, 'float');  % User Variable 9
    image.user10 = fread(id, 1, 'float');  % User Variable 10
    image.user11 = fread(id, 1, 'float');  % User Variable 11
    image.user12 = fread(id, 1, 'float');  % User Variable 12
    image.user13 = fread(id, 1, 'float');  % User Variable 13
    image.user14 = fread(id, 1, 'float');  % User Variable 14
    image.user15 = fread(id, 1, 'float');  % User Variable 15
    image.user16 = fread(id, 1, 'float');  % User Variable 16
    image.user17 = fread(id, 1, 'float');  % User Variable 17
    image.user18 = fread(id, 1, 'float');  % User Variable 18
    image.user19 = fread(id, 1, 'float');  % User Variable 19
    image.user20 = fread(id, 1, 'float');  % User Variable 20
    image.user21 = fread(id, 1, 'float');  % User Variable 21
    image.user22 = fread(id, 1, 'float');  % User Variable 22
    image.user23 = fread(id, 1, 'float');  % Projection Angle
    image.user24 = fread(id, 1, 'float');  % Concat Sat Type Flag
    image.im_alloc_key = modchar(fread(id, 13, '*char'));  %
    fseek(id, 3, 0);  % for struct byte alignment
    image.im_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    image.im_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    image.im_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    image.im_pds_a = fread(id, 1, 'int');  % PixelData size - as stored
    image.im_pds_c = fread(id, 1, 'int');  % PixelData size - Compressed
    image.im_pds_u = fread(id, 1, 'int');  % PixelData size - UnCompressed
    image.im_checksum = fread(id, 1, 'ulong');  % AcqRecon record checksum
    image.im_archived = fread(id, 1, 'long');  % Image Archive Flag
    image.im_complete = fread(id, 1, 'long');  % Image Complete Flag
    image.satbits = fread(id, 1, 'short');  % Bitmap of SAT selections
    image.scic = fread(id, 1, 'short');  % Surface Coil Intensity Correction Flag
    image.satxloc1 = fread(id, 1, 'short');  % R-side SAT pulse loc rel to lndmrk
    image.satxloc2 = fread(id, 1, 'short');  % L-side SAT pulse loc rel to lndmrk
    image.satyloc1 = fread(id, 1, 'short');  % A-side SAT pulse loc rel to lndmrk
    image.satyloc2 = fread(id, 1, 'short');  % P-side SAT pulse loc rel to lndmrk
    image.satzloc1 = fread(id, 1, 'short');  % S-side SAT pulse loc rel to lndmrk
    image.satzloc2 = fread(id, 1, 'short');  % I-side SAT pulse loc rel to lndmrk
    image.satxthick = fread(id, 1, 'short');  % Thickness of X-axis SAT pulse
    image.satythick = fread(id, 1, 'short');  % Thickness of Y-axis SAT pulse
    image.satzthick = fread(id, 1, 'short');  % Thickness of Z-axis SAT pulse
    image.flax = fread(id, 1, 'short');  % Phase contrast flow axis
    image.venc = fread(id, 1, 'short');  % Phase contrast velocity encoding
    image.thk_disclmr = fread(id, 1, 'short');  % Slice Thickness
    image.ps_flag = fread(id, 1, 'short');  % Auto/Manual Prescan flag
    image.ps_status = fread(id, 1, 'short');  % Bitmap of changed values
    image.image_type = fread(id, 1, 'short');  % Magnitude, Phase, Imaginary, or Real
    image.vas_collapse = fread(id, 1, 'short');  % Collapse Image
    image.user23n = fread(id, 1, 'float');  % User Variable 23
    image.user24n = fread(id, 1, 'float');  % User Variable 24
    image.proj_alg = fread(id, 1, 'short');  % Projection Algorithm
    image.proj_name = modchar(fread(id, 13, '*char'));  % Projection Algorithm Name
    fseek(id, 1, 0);  % for struct byte alignment
    image.x_axis_rot = fread(id, 1, 'float');  % X Axis Rotation
    image.y_axis_rot = fread(id, 1, 'float');  % Y Axis Rotation
    image.z_axis_rot = fread(id, 1, 'float');  % Z Axis Rotation
    image.thresh_min1 = fread(id, 1, 'int');  % Lower Range of Pixels 1
    image.thresh_max1 = fread(id, 1, 'int');  % Upper Range of Pixels 1
    image.thresh_min2 = fread(id, 1, 'int');  % Lower Range of Pixels 2
    image.thresh_max2 = fread(id, 1, 'int');  % Upper Range of Pixels 2
    image.echo_trn_len = fread(id, 1, 'short');  % Echo Train Length for Fast Spin Echo
    image.frac_echo = fread(id, 1, 'short');  % Fractional Echo - Effective TE Flag
    image.prep_pulse = fread(id, 1, 'short');  % Preporatory Pulse Option
    image.cphasenum = fread(id, 1, 'short');  % Cardiac Phase Number
    image.var_echo = fread(id, 1, 'short');  % Variable Echo Flag
    image.ref_img = modchar(fread(id, 1, '*char'));  % Reference Image Field
    image.sum_img = modchar(fread(id, 1, '*char'));  % Summary Image Field
    image.img_window = fread(id, 1, 'ushort');  % Window Value
    image.img_level = fread(id, 1, 'short');  % Level Value
    image.slop_int_1 = fread(id, 1, 'int');  % Number of 3D Slabs
    image.slop_int_2 = fread(id, 1, 'int');  % Slice Locs Per 3D Slab
    image.slop_int_3 = fread(id, 1, 'int');  % # of Slice Locs on Each Slab Which Overlap N eighbors
    image.slop_int_4 = fread(id, 1, 'int');  % Image Filtering 0.5/0.2T
    image.slop_int_5 = fread(id, 1, 'int');  % Integer Slop Field 5
    image.slop_float_1 = fread(id, 1, 'float');  % Float Slop Field 1
    image.slop_float_2 = fread(id, 1, 'float');  % Float Slop Field 2
    image.slop_float_3 = fread(id, 1, 'float');  % Float Slop Field 3
    image.slop_float_4 = fread(id, 1, 'float');  % Float Slop Field 4
    image.slop_float_5 = fread(id, 1, 'float');  % Float Slop Field 5
    image.slop_str_1 = modchar(fread(id, 16, '*char'));  % String Slop Field 1
    image.slop_str_2 = modchar(fread(id, 16, '*char'));  % String Slop Field 2
    image.scanactno = fread(id, 1, 'short');  % Scan Acquisition Number
    image.vasflags = fread(id, 1, 'short');  % Magnitude Weighting Flag
    image.vencscale = fread(id, 1, 'float');  % Scale Weighted Venc
    image.integrity = fread(id, 1, 'short');  % GE Image Integrity
    fseek(id, 2, 0);  % for struct byte alignment
    image.fphase = fread(id, 1, 'int');  % Number Of Phases
    image.freq_dir = fread(id, 1, 'short');  % Frequency Direction
    image.vas_mode = fread(id, 1, 'short');  % Vascular Mode
    image.image_uid = moduid(fread(id, 32, 'uint8'));  % Image Unique ID
    image.sop_uid = moduid(fread(id, 32, 'uint8'));  % Service Obj Class Unique ID
    image.dont_use_1 = fread(id, 1, 'short');  % This field is not used
    image.dont_use_2 = fread(id, 1, 'short');  % This field is not used
    image.dont_use_3 = fread(id, 1, 'short');  % This field is not used
    image.pscopts = fread(id, 1, 'short');  % bitmap of prescan options
    image.asoffsetx = fread(id, 1, 'short');  % gradient offset in X-direction
    image.asoffsety = fread(id, 1, 'short');  % gradient offset in Y-direction
    image.asoffsetz = fread(id, 1, 'short');  % gradient offset in Z-direction
    image.unoriginal = fread(id, 1, 'short');  % identifies image as original or unoriginal
    image.interleaves = fread(id, 1, 'short');  % number of EPI shots
    image.effechospace = fread(id, 1, 'short');  % effective echo spacing for EPI
    image.viewsperseg = fread(id, 1, 'short');  % views per segment
    image.rbpm = fread(id, 1, 'short');  % respiratory rate, breaths per min
    image.rtpoint = fread(id, 1, 'short');  % respiratory trigger point as percent of max.
    image.rcvrtype = fread(id, 1, 'short');  % type of receiver used
    image.dbdt = fread(id, 1, 'float');  % peak rate of change of gradient field, tesla/sec
    image.dbdtper = fread(id, 1, 'float');  % limit in units of percent of theoretical curve
    image.estdbdtper = fread(id, 1, 'float');  % PSD estimated limit in units of percent
    image.estdbdtts = fread(id, 1, 'float');  % PSD estimated limit in Teslas/sec
    image.saravghead = fread(id, 1, 'float');  % Avg head SAR
    image.neg_scanspacing = fread(id, 1, 'float');  % Negative scan spacing for overlap slices
    image.offsetfreq = fread(id, 1, 'int');  % Offset Frequency - Mag.Transfer
    image.user_usage_tag = fread(id, 1, 'ulong');  % Defines how following user CVs are to be filled in
    image.user_fill_mapMSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user_fill_mapLSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user25 = fread(id, 1, 'float');  % User Variable 25
    image.user26 = fread(id, 1, 'float');  % User Variable 26
    image.user27 = fread(id, 1, 'float');  % User Variable 27
    image.user28 = fread(id, 1, 'float');  % User Variable 28
    image.user29 = fread(id, 1, 'float');  % User Variable 29
    image.user30 = fread(id, 1, 'float');  % User Variable 30
    image.user31 = fread(id, 1, 'float');  % User Variable 31
    image.user32 = fread(id, 1, 'float');  % User Variable 32
    image.user33 = fread(id, 1, 'float');  % User Variable 33
    image.user34 = fread(id, 1, 'float');  % User Variable 34
    image.user35 = fread(id, 1, 'float');  % User Variable 35
    image.user36 = fread(id, 1, 'float');  % User Variable 36
    image.user37 = fread(id, 1, 'float');  % User Variable 37
    image.user38 = fread(id, 1, 'float');  % User Variable 38
    image.user39 = fread(id, 1, 'float');  % User Variable 39
    image.user40 = fread(id, 1, 'float');  % User Variable 40
    image.user41 = fread(id, 1, 'float');  % User Variable 41
    image.user42 = fread(id, 1, 'float');  % User Variable 42
    image.user43 = fread(id, 1, 'float');  % User Variable 43
    image.user44 = fread(id, 1, 'float');  % User Variable 44
    image.user45 = fread(id, 1, 'float');  % User Variable 45
    image.user46 = fread(id, 1, 'float');  % User Variable 46
    image.user47 = fread(id, 1, 'float');  % User Variable 47
    image.user48 = fread(id, 1, 'float');  % User Variable 48
    image.slop_int_6 = fread(id, 1, 'int');  % Integer Slop Field 6
    image.slop_int_7 = fread(id, 1, 'int');  % Integer Slop Field 7
    image.slop_int_8 = fread(id, 1, 'int');  % Integer Slop Field 8
    image.slop_int_9 = fread(id, 1, 'int');  % Integer Slop Field 9
    image.slop_int_10 = fread(id, 1, 'int');  % Integer Slop Field 10
    image.slop_int_11 = fread(id, 1, 'int');  % Integer Slop Field 11
    image.slop_int_12 = fread(id, 1, 'int');  % in use by MR-YMS
    image.slop_int_13 = fread(id, 1, 'int');  % Integer Slop Field 13
    image.slop_int_14 = fread(id, 1, 'int');  % Integer Slop Field 14
    image.slop_int_15 = fread(id, 1, 'int');  % Station Index
    image.slop_int_16 = fread(id, 1, 'int');  % Station Total
    image.slop_int_17 = fread(id, 1, 'int');  % the LAST Spare integer field!
elseif (ver >= 9.0 && ver < 11.0)
    image.im_suid = modchar(fread(id, 4, '*char'));  % Suite id for this image
    image.im_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    image.im_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Image
    fseek(id, 1, 0);  % for struct byte alignment
    image.im_exno = fread(id, 1, 'ushort');  % Exam number for this image
    image.im_seno = fread(id, 1, 'short');  % Series Number for this image
    image.im_no = fread(id, 1, 'short');  % Image Number
    fseek(id, 2, 0);  % for struct byte alignment
    image.im_datetime = fread(id, 1, 'int');  % Allocation Image date/time stamp
    image.im_actual_dt = fread(id, 1, 'int');  % Actual Image date/time stamp
    image.sctime = fread(id, 1, 'float');  % Duration of scan
    image.slthick = fread(id, 1, 'float');  % Slice Thickness (mm)
    image.imatrix_X = fread(id, 1, 'short');  % Image matrix size - X
    image.imatrix_Y = fread(id, 1, 'short');  % Image matrix size - Y
    image.dfov = fread(id, 1, 'float');  % Display field of view - X (mm)
    image.dfov_rect = fread(id, 1, 'float');  % Display field of view - Y (if different)
    image.dim_X = fread(id, 1, 'float');  % Image dimension - X
    image.dim_Y = fread(id, 1, 'float');  % Image dimension - Y
    image.pixsize_X = fread(id, 1, 'float');  % Image pixel size - X
    image.pixsize_Y = fread(id, 1, 'float');  % Image pixel size - Y
    image.pdid = modchar(fread(id, 14, '*char'));  % Pixel Data ID
    image.contrastIV = modchar(fread(id, 17, '*char'));  % IV Contrast Agent
    image.contrastOral = modchar(fread(id, 17, '*char'));  % Oral Contrast Agent
    image.contmode = fread(id, 1, 'short');  % Image Contrast Mode
    image.serrx = fread(id, 1, 'short');  % Series from which prescribed
    image.imgrx = fread(id, 1, 'short');  % Image from which prescribed
    image.screenformat = fread(id, 1, 'short');  % Screen Format(8/16 bit)
    image.plane = fread(id, 1, 'short');  % Plane Type
    fseek(id, 2, 0);  % for struct byte alignment
    image.scanspacing = fread(id, 1, 'float');  % Spacing between scans (mm?)
    image.im_compress = fread(id, 1, 'short');  % Image compression type for allocation
    image.im_scouttype = fread(id, 1, 'short');  % Scout Type (AP or lateral)
    image.loc_ras = modchar(fread(id, 1, '*char'));  % RAS letter of image location
    fseek(id, 3, 0);  % for struct byte alignment
    image.loc = fread(id, 1, 'float');  % Image location
    image.ctr_R = fread(id, 1, 'float');  % Center R coord of plane image
    image.ctr_A = fread(id, 1, 'float');  % Center A coord of plane image
    image.ctr_S = fread(id, 1, 'float');  % Center S coord of plane image
    image.norm_R = fread(id, 1, 'float');  % Normal R coord
    image.norm_A = fread(id, 1, 'float');  % Normal A coord
    image.norm_S = fread(id, 1, 'float');  % Normal S coord
    image.tlhc_R = fread(id, 1, 'float');  % R Coord of Top Left Hand Corner
    image.tlhc_A = fread(id, 1, 'float');  % A Coord of Top Left Hand Corner
    image.tlhc_S = fread(id, 1, 'float');  % S Coord of Top Left Hand Corner
    image.trhc_R = fread(id, 1, 'float');  % R Coord of Top Right Hand Corner
    image.trhc_A = fread(id, 1, 'float');  % A Coord of Top Right Hand Corner
    image.trhc_S = fread(id, 1, 'float');  % S Coord of Top Right Hand Corner
    image.brhc_R = fread(id, 1, 'float');  % R Coord of Bottom Right Hand Corner
    image.brhc_A = fread(id, 1, 'float');  % A Coord of Bottom Right Hand Corner
    image.brhc_S = fread(id, 1, 'float');  % S Coord of Bottom Right Hand Corner
    image.forimgrev = modchar(fread(id, 4, '*char'));  % Foreign Image Revision
    image.tr = fread(id, 1, 'int');  % Pulse repetition time(usec)
    image.ti = fread(id, 1, 'int');  % Pulse inversion time(usec)
    image.te = fread(id, 1, 'int');  % Pulse echo time(usec)
    image.te2 = fread(id, 1, 'int');  % Second echo echo (usec)
    image.numecho = fread(id, 1, 'short');  % Number of echoes
    image.echonum = fread(id, 1, 'short');  % Echo Number
    image.tbldlta = fread(id, 1, 'float');  % Table Delta
    image.nex = fread(id, 1, 'float');  % Number of Excitations
    image.contig = fread(id, 1, 'short');  % Continuous Slices Flag
    image.hrtrate = fread(id, 1, 'short');  % Cardiac Heart Rate (bpm)
    image.tdel = fread(id, 1, 'int');  % Delay time after trigger (msec)
    image.saravg = fread(id, 1, 'float');  % Average SAR
    image.sarpeak = fread(id, 1, 'float');  % Peak SAR
    image.monsar = fread(id, 1, 'short');  % Monitor SAR flag
    image.trgwindow = fread(id, 1, 'short');  % Trigger window (% of R-R interval)
    image.reptime = fread(id, 1, 'float');  % Cardiac repetition time
    image.imgpcyc = fread(id, 1, 'short');  % Images per cardiac cycle
    image.xmtgain = fread(id, 1, 'short');  % Actual Transmit Gain (.1 db)
    image.rcvgain1 = fread(id, 1, 'short');  % Actual Receive Gain Analog (.1 db)
    image.rcvgain2 = fread(id, 1, 'short');  % Actual Receive Gain Digital (.1 db)
    image.mr_flip = fread(id, 1, 'short');  % Flip Angle for GRASS scans (deg.)
    fseek(id, 2, 0);  % for struct byte alignment
    image.mindat = fread(id, 1, 'int');  % Minimum Delay after Trigger (uSec)
    image.cphase = fread(id, 1, 'short');  % Total Cardiac Phase prescribed
    image.swappf = fread(id, 1, 'short');  % Swap Phase/Frequency Axis
    image.pauseint = fread(id, 1, 'short');  % Pause Interval (slices)
    fseek(id, 2, 0);  % for struct byte alignment
    image.pausetime = fread(id, 1, 'float');  % Pause Time
    image.obplane = fread(id, 1, 'int');  % Oblique Plane
    image.slocfov = fread(id, 1, 'int');  % Slice Offsets on Freq axis
    image.xmtfreq = fread(id, 1, 'int');  % Center Frequency (0.1 Hz)
    image.autoxmtfreq = fread(id, 1, 'int');  % Auto Center Frequency (0.1 Hz)
    image.autoxmtgain = fread(id, 1, 'short');  % Auto Transmit Gain (0.1 dB)
    image.prescan_r1 = fread(id, 1, 'short');  % PreScan R1 - Analog
    image.prescan_r2 = fread(id, 1, 'short');  % PreScan R2 - Digital
    fseek(id, 2, 0);  % for struct byte alignment
    image.user_bitmap = fread(id, 1, 'int');  % Bitmap defining user CVs
    image.cenfreq = fread(id, 1, 'short');  % Center Frequency Method
    image.imode = fread(id, 1, 'short');  % Imaging Mode
    image.iopt = fread(id, 1, 'int');  % Imaging Options
    image.pseq = fread(id, 1, 'short');  % Pulse Sequence
    image.pseqmode = fread(id, 1, 'short');  % Pulse Sequence Mode
    image.psdname = modchar(fread(id, 33, '*char'));  % Pulse Sequence Name
    fseek(id, 3, 0);  % for struct byte alignment
    image.psd_datetime = fread(id, 1, 'int');  % PSD Creation Date and Time
    image.psd_iname = modchar(fread(id, 13, '*char'));  % PSD name from inside PSD
    fseek(id, 1, 0);  % for struct byte alignment
    image.ctyp = fread(id, 1, 'short');  % Coil Type
    image.cname = modchar(fread(id, 17, '*char'));  % Coil Name
    fseek(id, 1, 0);  % for struct byte alignment
    image.surfctyp = fread(id, 1, 'short');  % Surface Coil Type
    image.surfcext = fread(id, 1, 'short');  % Extremity Coil Flag
    fseek(id, 2, 0);  % for struct byte alignment
    image.rawrunnum = fread(id, 1, 'int');  % RawData Run Number
    image.cal_fldstr = fread(id, 1, 'ulong');  % Calibrated Field Strength (x10 uGauss)
    image.supp_tech = fread(id, 1, 'short');  % SAT fat/water/none
    fseek(id, 2, 0);  % for struct byte alignment
    image.vbw = fread(id, 1, 'float');  % Variable Bandwidth (Hz)
    image.slquant = fread(id, 1, 'short');  % Number of slices in this scan group
    image.gpre = fread(id, 1, 'short');  % Graphically prescribed
    image.intr_del = fread(id, 1, 'int');  % Interimage/interloc delay (uSec)
    image.user0 = fread(id, 1, 'float');  % User Variable 0
    image.user1 = fread(id, 1, 'float');  % User Variable 1
    image.user2 = fread(id, 1, 'float');  % User Variable 2
    image.user3 = fread(id, 1, 'float');  % User Variable 3
    image.user4 = fread(id, 1, 'float');  % User Variable 4
    image.user5 = fread(id, 1, 'float');  % User Variable 5
    image.user6 = fread(id, 1, 'float');  % User Variable 6
    image.user7 = fread(id, 1, 'float');  % User Variable 7
    image.user8 = fread(id, 1, 'float');  % User Variable 8
    image.user9 = fread(id, 1, 'float');  % User Variable 9
    image.user10 = fread(id, 1, 'float');  % User Variable 10
    image.user11 = fread(id, 1, 'float');  % User Variable 11
    image.user12 = fread(id, 1, 'float');  % User Variable 12
    image.user13 = fread(id, 1, 'float');  % User Variable 13
    image.user14 = fread(id, 1, 'float');  % User Variable 14
    image.user15 = fread(id, 1, 'float');  % User Variable 15
    image.user16 = fread(id, 1, 'float');  % User Variable 16
    image.user17 = fread(id, 1, 'float');  % User Variable 17
    image.user18 = fread(id, 1, 'float');  % User Variable 18
    image.user19 = fread(id, 1, 'float');  % User Variable 19
    image.user20 = fread(id, 1, 'float');  % User Variable 20
    image.user21 = fread(id, 1, 'float');  % User Variable 21
    image.user22 = fread(id, 1, 'float');  % User Variable 22
    image.proj_ang = fread(id, 1, 'float');  % Projection Angle
    image.concat_sat = fread(id, 1, 'float');  % Concat Sat Type Flag
    image.im_alloc_key = modchar(fread(id, 13, '*char'));  %
    fseek(id, 3, 0);  % for struct byte alignment
    image.im_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    image.im_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    image.im_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    image.im_pds_a = fread(id, 1, 'int');  % PixelData size - as stored
    image.im_pds_c = fread(id, 1, 'int');  % PixelData size - Compressed
    image.im_pds_u = fread(id, 1, 'int');  % PixelData size - UnCompressed
    image.im_checksum = fread(id, 1, 'ulong');  % AcqRecon record checksum
    image.im_archived = fread(id, 1, 'long');  % Image Archive Flag
    image.im_complete = fread(id, 1, 'long');  % Image Complete Flag
    image.satbits = fread(id, 1, 'short');  % Bitmap of SAT selections
    image.scic = fread(id, 1, 'short');  % Surface Coil Intensity Correction Flag
    image.satxloc1 = fread(id, 1, 'short');  % R-side SAT pulse loc rel to lndmrk
    image.satxloc2 = fread(id, 1, 'short');  % L-side SAT pulse loc rel to lndmrk
    image.satyloc1 = fread(id, 1, 'short');  % A-side SAT pulse loc rel to lndmrk
    image.satyloc2 = fread(id, 1, 'short');  % P-side SAT pulse loc rel to lndmrk
    image.satzloc1 = fread(id, 1, 'short');  % S-side SAT pulse loc rel to lndmrk
    image.satzloc2 = fread(id, 1, 'short');  % I-side SAT pulse loc rel to lndmrk
    image.satxthick = fread(id, 1, 'short');  % Thickness of X-axis SAT pulse
    image.satythick = fread(id, 1, 'short');  % Thickness of Y-axis SAT pulse
    image.satzthick = fread(id, 1, 'short');  % Thickness of Z-axis SAT pulse
    image.flax = fread(id, 1, 'short');  % Phase contrast flow axis
    image.venc = fread(id, 1, 'short');  % Phase contrast velocity encoding
    image.thk_disclmr = fread(id, 1, 'short');  % Slice Thickness
    image.ps_flag = fread(id, 1, 'short');  % Auto/Manual Prescan flag
    image.ps_status = fread(id, 1, 'short');  % Bitmap of changed values
    image.image_type = fread(id, 1, 'short');  % Magnitude, Phase, Imaginary, or Real
    image.vas_collapse = fread(id, 1, 'short');  % Collapse Image
    image.user23 = fread(id, 1, 'float');  % User Variable 23
    image.user24 = fread(id, 1, 'float');  % User Variable 24
    image.proj_alg = fread(id, 1, 'short');  % Projection Algorithm
    image.proj_name = modchar(fread(id, 13, '*char'));  % Projection Algorithm Name
    fseek(id, 1, 0);  % for struct byte alignment
    image.x_axis_rot = fread(id, 1, 'float');  % X Axis Rotation
    image.y_axis_rot = fread(id, 1, 'float');  % Y Axis Rotation
    image.z_axis_rot = fread(id, 1, 'float');  % Z Axis Rotation
    image.thresh_min1 = fread(id, 1, 'int');  % Lower Range of Pixels 1
    image.thresh_max1 = fread(id, 1, 'int');  % Upper Range of Pixels 1
    image.thresh_min2 = fread(id, 1, 'int');  % Lower Range of Pixels 2
    image.thresh_max2 = fread(id, 1, 'int');  % Upper Range of Pixels 2
    image.echo_trn_len = fread(id, 1, 'short');  % Echo Train Length for Fast Spin Echo
    image.frac_echo = fread(id, 1, 'short');  % Fractional Echo - Effective TE Flag
    image.prep_pulse = fread(id, 1, 'short');  % Preporatory Pulse Option
    image.cphasenum = fread(id, 1, 'short');  % Cardiac Phase Number
    image.var_echo = fread(id, 1, 'short');  % Variable Echo Flag
    image.ref_img = modchar(fread(id, 1, '*char'));  % Reference Image Field
    image.sum_img = modchar(fread(id, 1, '*char'));  % Summary Image Field
    image.img_window = fread(id, 1, 'ushort');  % Window Value
    image.img_level = fread(id, 1, 'short');  % Level Value
    image.numslabs = fread(id, 1, 'int');  % Number of 3D Slabs
    image.locsperslab = fread(id, 1, 'int');  % Slice Locs Per 3D Slab
    image.overlaps = fread(id, 1, 'int');  % # of Slice Locs on Each Slab Which Overlap N eighbors
    image.slop_int_4 = fread(id, 1, 'int');  % Image Filtering 0.5/0.2T
    image.dfax = fread(id, 1, 'int');  % Diffusion Direction for DW-EPI
    image.ihtagfa = fread(id, 1, 'float');  % Tagging Flip Angle
    image.ihtagor = fread(id, 1, 'float');  % Cardiac Tagging Orientation
    image.ihbspti = fread(id, 1, 'float');  % Blood Suppression TI
    image.rtia_timer = fread(id, 1, 'float');  % Float Slop Field 4
    image.fps = fread(id, 1, 'float');  % Float Slop Field 5
    image.filter_mode = modchar(fread(id, 16, '*char'));  % String Slop Field 1
    image.slop_str_2 = modchar(fread(id, 16, '*char'));  % String Slop Field 2
    image.scanactno = fread(id, 1, 'short');  % Scan Acquisition Number
    image.vasflags = fread(id, 1, 'short');  % Magnitude Weighting Flag
    image.vencscale = fread(id, 1, 'float');  % Scale Weighted Venc
    image.integrity = fread(id, 1, 'short');  % GE Image Integrity
    fseek(id, 2, 0);  % for struct byte alignment
    image.fphase = fread(id, 1, 'int');  % Number Of Phases
    image.freq_dir = fread(id, 1, 'short');  % Frequency Direction
    image.vas_mode = fread(id, 1, 'short');  % Vascular Mode
    image.image_uid = moduid(fread(id, 32, 'uint8'));  % Image Unique ID
    image.sop_uid = moduid(fread(id, 32, 'uint8'));  % Service Obj Class Unique ID
    image.dont_use_1 = fread(id, 1, 'short');  % This field is not used
    image.dont_use_2 = fread(id, 1, 'short');  % This field is not used
    image.dont_use_3 = fread(id, 1, 'short');  % This field is not used
    image.pscopts = fread(id, 1, 'short');  % bitmap of prescan options
    image.asoffsetx = fread(id, 1, 'short');  % gradient offset in X-direction
    image.asoffsety = fread(id, 1, 'short');  % gradient offset in Y-direction
    image.asoffsetz = fread(id, 1, 'short');  % gradient offset in Z-direction
    image.unoriginal = fread(id, 1, 'short');  % identifies image as original or unoriginal
    image.interleaves = fread(id, 1, 'short');  % number of EPI shots
    image.effechospace = fread(id, 1, 'short');  % effective echo spacing for EPI
    image.viewsperseg = fread(id, 1, 'short');  % views per segment
    image.rbpm = fread(id, 1, 'short');  % respiratory rate, breaths per min
    image.rtpoint = fread(id, 1, 'short');  % respiratory trigger point as percent of max.
    image.rcvrtype = fread(id, 1, 'short');  % type of receiver used
    image.dbdt = fread(id, 1, 'float');  % peak rate of change of gradient field, tesla/sec
    image.dbdtper = fread(id, 1, 'float');  % limit in units of percent of theoretical curve
    image.estdbdtper = fread(id, 1, 'float');  % PSD estimated limit in units of percent
    image.estdbdtts = fread(id, 1, 'float');  % PSD estimated limit in Teslas/sec
    image.saravghead = fread(id, 1, 'float');  % Avg head SAR
    image.neg_scanspacing = fread(id, 1, 'float');  % Negative scan spacing for overlap slices
    image.offsetfreq = fread(id, 1, 'int');  % Offset Frequency - Mag.Transfer
    image.user_usage_tag = fread(id, 1, 'ulong');  % Defines how following user CVs are to be filled in
    image.user_fill_mapMSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user_fill_mapLSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user25 = fread(id, 1, 'float');  % User Variable 25
    image.user26 = fread(id, 1, 'float');  % User Variable 26
    image.user27 = fread(id, 1, 'float');  % User Variable 27
    image.user28 = fread(id, 1, 'float');  % User Variable 28
    image.user29 = fread(id, 1, 'float');  % User Variable 29
    image.user30 = fread(id, 1, 'float');  % User Variable 30
    image.user31 = fread(id, 1, 'float');  % User Variable 31
    image.user32 = fread(id, 1, 'float');  % User Variable 32
    image.user33 = fread(id, 1, 'float');  % User Variable 33
    image.user34 = fread(id, 1, 'float');  % User Variable 34
    image.user35 = fread(id, 1, 'float');  % User Variable 35
    image.user36 = fread(id, 1, 'float');  % User Variable 36
    image.user37 = fread(id, 1, 'float');  % User Variable 37
    image.user38 = fread(id, 1, 'float');  % User Variable 38
    image.user39 = fread(id, 1, 'float');  % User Variable 39
    image.user40 = fread(id, 1, 'float');  % User Variable 40
    image.user41 = fread(id, 1, 'float');  % User Variable 41
    image.user42 = fread(id, 1, 'float');  % User Variable 42
    image.user43 = fread(id, 1, 'float');  % User Variable 43
    image.user44 = fread(id, 1, 'float');  % User Variable 44
    image.user45 = fread(id, 1, 'float');  % User Variable 45
    image.user46 = fread(id, 1, 'float');  % User Variable 46
    image.user47 = fread(id, 1, 'float');  % User Variable 47
    image.user48 = fread(id, 1, 'float');  % User Variable 48
    image.b_value = fread(id, 1, 'int');  % B-value for DW-EPI
    image.iopt2 = fread(id, 1, 'int');  % Imaging Option2
    image.ihtagging = fread(id, 1, 'int');  % tag type
    image.ihtagspc = fread(id, 1, 'int');  % tag space
    image.ihfcineim = fread(id, 1, 'int');  % Fast CINE interpolation method
    image.ihfcinent = fread(id, 1, 'int');  % Fast CINE normalization type
    image.num_seg = fread(id, 1, 'int');  % YMSge05074
    image.oprtarr = fread(id, 1, 'int');  % Respiratory Trigger windo
    image.averages = fread(id, 1, 'int');  % Number of averages for spectro
    image.station_index = fread(id, 1, 'int');  % Station Index
    image.station_total = fread(id, 1, 'int');  % Station Total
    image.slop_int_17 = fread(id, 1, 'int');  % the LAST Spare integer field!
    image.RegressorVal = fread(id, 1, 'float');  %
    image.delAcq = fread(id, 1, 'int');  % Delay after Acquisition (MP / fMRI screen)
    image.img_hdr_padding = modchar(fread(id, 484, '*char'));  %
elseif (ver == 11)
    image.dfov = fread(id, 1, 'float');  % Display field of view - X (mm)
    image.dfov_rect = fread(id, 1, 'float');  % Display field of view - Y (if different)
    image.sctime = fread(id, 1, 'float');  % Duration of scan
    image.slthick = fread(id, 1, 'float');  % Slice Thickness (mm)
    image.scanspacing = fread(id, 1, 'float');  % Spacing between scans (mm?)
    image.loc = fread(id, 1, 'float');  % Image location
    image.tbldlta = fread(id, 1, 'float');  % Table Delta
    image.nex = fread(id, 1, 'float');  % Number of Excitations
    image.reptime = fread(id, 1, 'float');  % Cardiac repetition time
    image.saravg = fread(id, 1, 'float');  % Average SAR
    image.sarpeak = fread(id, 1, 'float');  % Peak SAR
    image.pausetime = fread(id, 1, 'float');  % Pause Time
    image.vbw = fread(id, 1, 'float');  % Variable Bandwidth (Hz)
    image.user0 = fread(id, 1, 'float');  % User Variable 0
    image.user1 = fread(id, 1, 'float');  % User Variable 1
    image.user2 = fread(id, 1, 'float');  % User Variable 2
    image.user3 = fread(id, 1, 'float');  % User Variable 3
    image.user4 = fread(id, 1, 'float');  % User Variable 4
    image.user5 = fread(id, 1, 'float');  % User Variable 5
    image.user6 = fread(id, 1, 'float');  % User Variable 6
    image.user7 = fread(id, 1, 'float');  % User Variable 7
    image.user8 = fread(id, 1, 'float');  % User Variable 8
    image.user9 = fread(id, 1, 'float');  % User Variable 9
    image.user10 = fread(id, 1, 'float');  % User Variable 10
    image.user11 = fread(id, 1, 'float');  % User Variable 11
    image.user12 = fread(id, 1, 'float');  % User Variable 12
    image.user13 = fread(id, 1, 'float');  % User Variable 13
    image.user14 = fread(id, 1, 'float');  % User Variable 14
    image.user15 = fread(id, 1, 'float');  % User Variable 15
    image.user16 = fread(id, 1, 'float');  % User Variable 16
    image.user17 = fread(id, 1, 'float');  % User Variable 17
    image.user18 = fread(id, 1, 'float');  % User Variable 18
    image.user19 = fread(id, 1, 'float');  % User Variable 19
    image.user20 = fread(id, 1, 'float');  % User Variable 20
    image.user21 = fread(id, 1, 'float');  % User Variable 21
    image.user22 = fread(id, 1, 'float');  % User Variable 22
    image.proj_ang = fread(id, 1, 'float');  % Projection Angle
    image.concat_sat = fread(id, 1, 'float');  % Concat Sat Type Flag
    image.user23 = fread(id, 1, 'float');  % User Variable 23
    image.user24 = fread(id, 1, 'float');  % User Variable 24
    image.x_axis_rot = fread(id, 1, 'float');  % X Axis Rotation
    image.y_axis_rot = fread(id, 1, 'float');  % Y Axis Rotation
    image.z_axis_rot = fread(id, 1, 'float');  % Z Axis Rotation
    image.ihtagfa = fread(id, 1, 'float');  % Tagging Flip Angle
    image.ihtagor = fread(id, 1, 'float');  % Cardiac Tagging Orientation
    image.ihbspti = fread(id, 1, 'float');  % Blood Suppression TI
    image.rtia_timer = fread(id, 1, 'float');  % Float Slop Field 4
    image.fps = fread(id, 1, 'float');  % Float Slop Field 5
    image.vencscale = fread(id, 1, 'float');  % Scale Weighted Venc
    image.dbdt = fread(id, 1, 'float');  % peak rate of change of gradient field, tesla/sec
    image.dbdtper = fread(id, 1, 'float');  % limit in units of percent of theoretical curve
    image.estdbdtper = fread(id, 1, 'float');  % PSD estimated limit in units of percent
    image.estdbdtts = fread(id, 1, 'float');  % PSD estimated limit in Teslas/sec
    image.saravghead = fread(id, 1, 'float');  % Avg head SAR
    image.neg_scanspacing = fread(id, 1, 'float');  % Negative scan spacing for overlap slices
    image.user25 = fread(id, 1, 'float');  % User Variable 25
    image.user26 = fread(id, 1, 'float');  % User Variable 26
    image.user27 = fread(id, 1, 'float');  % User Variable 27
    image.user28 = fread(id, 1, 'float');  % User Variable 28
    image.user29 = fread(id, 1, 'float');  % User Variable 29
    image.user30 = fread(id, 1, 'float');  % User Variable 30
    image.user31 = fread(id, 1, 'float');  % User Variable 31
    image.user32 = fread(id, 1, 'float');  % User Variable 32
    image.user33 = fread(id, 1, 'float');  % User Variable 33
    image.user34 = fread(id, 1, 'float');  % User Variable 34
    image.user35 = fread(id, 1, 'float');  % User Variable 35
    image.user36 = fread(id, 1, 'float');  % User Variable 36
    image.user37 = fread(id, 1, 'float');  % User Variable 37
    image.user38 = fread(id, 1, 'float');  % User Variable 38
    image.user39 = fread(id, 1, 'float');  % User Variable 39
    image.user40 = fread(id, 1, 'float');  % User Variable 40
    image.user41 = fread(id, 1, 'float');  % User Variable 41
    image.user42 = fread(id, 1, 'float');  % User Variable 42
    image.user43 = fread(id, 1, 'float');  % User Variable 43
    image.user44 = fread(id, 1, 'float');  % User Variable 44
    image.user45 = fread(id, 1, 'float');  % User Variable 45
    image.user46 = fread(id, 1, 'float');  % User Variable 46
    image.user47 = fread(id, 1, 'float');  % User Variable 47
    image.user48 = fread(id, 1, 'float');  % User Variable 48
    image.RegressorVal = fread(id, 1, 'float');  %
    image.SliceAsset = fread(id, 1, 'float');  % Slice Asset in Asset Screen
    image.PhaseAsset = fread(id, 1, 'float');  % Phase Asset in Asset Screen
    image.sarValues = fread(id, 4, 'float');  % correspoding SAR values for defined terms
    image.shim_fov = fread(id, 2, 'float');  %
    image.shim_ctr_R = fread(id, 2, 'float');  %
    image.shim_ctr_A = fread(id, 2, 'float');  %
    image.shim_ctr_S = fread(id, 2, 'float');  %
    image.dim_X = fread(id, 1, 'float');  % Image dimension - X
    image.dim_Y = fread(id, 1, 'float');  % Image dimension - Y
    image.pixsize_X = fread(id, 1, 'float');  % Image pixel size - X
    image.pixsize_Y = fread(id, 1, 'float');  % Image pixel size - Y
    image.ctr_R = fread(id, 1, 'float');  % Center R coord of plane image
    image.ctr_A = fread(id, 1, 'float');  % Center A coord of plane image
    image.ctr_S = fread(id, 1, 'float');  % Center S coord of plane image
    image.norm_R = fread(id, 1, 'float');  % Normal R coord
    image.norm_A = fread(id, 1, 'float');  % Normal A coord
    image.norm_S = fread(id, 1, 'float');  % Normal S coord
    image.tlhc_R = fread(id, 1, 'float');  % R Coord of Top Left Hand Corner
    image.tlhc_A = fread(id, 1, 'float');  % A Coord of Top Left Hand Corner
    image.tlhc_S = fread(id, 1, 'float');  % S Coord of Top Left Hand Corner
    image.trhc_R = fread(id, 1, 'float');  % R Coord of Top Right Hand Corner
    image.trhc_A = fread(id, 1, 'float');  % A Coord of Top Right Hand Corner
    image.trhc_S = fread(id, 1, 'float');  % S Coord of Top Right Hand Corner
    image.brhc_R = fread(id, 1, 'float');  % R Coord of Bottom Right Hand Corner
    image.brhc_A = fread(id, 1, 'float');  % A Coord of Bottom Right Hand Corner
    image.brhc_S = fread(id, 1, 'float');  % S Coord of Bottom Right Hand Corner
    image.cal_fldstr = fread(id, 1, 'ulong');  % Calibrated Field Strength (x10 uGauss)
    image.im_checksum = fread(id, 1, 'ulong');  % AcqRecon record checksum
    image.user_usage_tag = fread(id, 1, 'ulong');  % Defines how following user CVs are to be filled in
    image.user_fill_mapMSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user_fill_mapLSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.im_archived = fread(id, 1, 'long');  % Image Archive Flag
    image.im_complete = fread(id, 1, 'long');  % Image Complete Flag
    image.im_datetime = fread(id, 1, 'int');  % Allocation Image date/time stamp
    image.im_actual_dt = fread(id, 1, 'int');  % Actual Image date/time stamp
    image.tr = fread(id, 1, 'int');  % Pulse repetition time(usec)
    image.ti = fread(id, 1, 'int');  % Pulse inversion time(usec)
    image.te = fread(id, 1, 'int');  % Pulse echo time(usec)
    image.te2 = fread(id, 1, 'int');  % Second echo echo (usec)
    image.tdel = fread(id, 1, 'int');  % Delay time after trigger (msec)
    image.mindat = fread(id, 1, 'int');  % Minimum Delay after Trigger (uSec)
    image.obplane = fread(id, 1, 'int');  % Oblique Plane
    image.slocfov = fread(id, 1, 'int');  % Slice Offsets on Freq axis
    image.xmtfreq = fread(id, 1, 'int');  % Center Frequency (0.1 Hz)
    image.autoxmtfreq = fread(id, 1, 'int');  % Auto Center Frequency (0.1 Hz)
    image.user_bitmap = fread(id, 1, 'int');  % Bitmap defining user CVs
    image.iopt = fread(id, 1, 'int');  % Imaging Options
    image.psd_datetime = fread(id, 1, 'int');  % PSD Creation Date and Time
    image.rawrunnum = fread(id, 1, 'int');  % RawData Run Number
    image.intr_del = fread(id, 1, 'int');  % Interimage/interloc delay (uSec)
    image.im_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    image.im_pds_a = fread(id, 1, 'int');  % PixelData size - as stored
    image.im_pds_c = fread(id, 1, 'int');  % PixelData size - Compressed
    image.im_pds_u = fread(id, 1, 'int');  % PixelData size - UnCompressed
    image.thresh_min1 = fread(id, 1, 'int');  % Lower Range of Pixels 1
    image.thresh_max1 = fread(id, 1, 'int');  % Upper Range of Pixels 1
    image.thresh_min2 = fread(id, 1, 'int');  % Lower Range of Pixels 2
    image.thresh_max2 = fread(id, 1, 'int');  % Upper Range of Pixels 2
    image.numslabs = fread(id, 1, 'int');  % Number of 3D Slabs
    image.locsperslab = fread(id, 1, 'int');  % Slice Locs Per 3D Slab
    image.overlaps = fread(id, 1, 'int');  % # of Slice Locs on Each Slab Which Overlap N eighbors
    image.slop_int_4 = fread(id, 1, 'int');  % Image Filtering 0.5/0.2T
    image.dfax = fread(id, 1, 'int');  % Diffusion Direction for DW-EPI
    image.fphase = fread(id, 1, 'int');  % Number Of Phases
    image.offsetfreq = fread(id, 1, 'int');  % Offset Frequency - Mag.Transfer
    image.b_value = fread(id, 1, 'int');  % B-value for DW-EPI
    image.iopt2 = fread(id, 1, 'int');  % Imaging Option2
    image.ihtagging = fread(id, 1, 'int');  % tag type
    image.ihtagspc = fread(id, 1, 'int');  % tag space
    image.ihfcineim = fread(id, 1, 'int');  % Fast CINE interpolation method
    image.ihfcinent = fread(id, 1, 'int');  % Fast CINE normalization type
    image.num_seg = fread(id, 1, 'int');  % YMSge05074
    image.oprtarr = fread(id, 1, 'int');  % Respiratory Trigger windo
    image.averages = fread(id, 1, 'int');  % Number of averages for spectro
    image.station_index = fread(id, 1, 'int');  % Station Index
    image.station_total = fread(id, 1, 'int');  % Station Total
    image.iopt3 = fread(id, 1, 'int');  % Imaging Option3
    image.delAcq = fread(id, 1, 'int');  % Delay after Acquisition (MP / fMRI screen)
    image.imatrix_X = fread(id, 1, 'short');  % Image matrix size - X
    image.imatrix_Y = fread(id, 1, 'short');  % Image matrix size - Y
    image.im_exno = fread(id, 1, 'ushort');  % Exam number for this image
    image.img_window = fread(id, 1, 'ushort');  % Window Value
    image.img_level = fread(id, 1, 'short');  % Level Value
    image.numecho = fread(id, 1, 'short');  % Number of echoes
    image.echonum = fread(id, 1, 'short');  % Echo Number
    image.im_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    image.im_seno = fread(id, 1, 'short');  % Series Number for this image
    image.im_no = fread(id, 1, 'short');  % Image Number
    image.contmode = fread(id, 1, 'short');  % Image Contrast Mode
    image.serrx = fread(id, 1, 'short');  % Series from which prescribed
    image.imgrx = fread(id, 1, 'short');  % Image from which prescribed
    image.screenformat = fread(id, 1, 'short');  % Screen Format(8/16 bit)
    image.plane = fread(id, 1, 'short');  % Plane Type
    image.im_compress = fread(id, 1, 'short');  % Image compression type for allocation
    image.im_scouttype = fread(id, 1, 'short');  % Scout Type (AP or lateral)
    image.contig = fread(id, 1, 'short');  % Continuous Slices Flag
    image.hrtrate = fread(id, 1, 'short');  % Cardiac Heart Rate (bpm)
    image.trgwindow = fread(id, 1, 'short');  % Trigger window (% of R-R interval)
    image.imgpcyc = fread(id, 1, 'short');  % Images per cardiac cycle
    image.xmtgain = fread(id, 1, 'short');  % Actual Transmit Gain (.1 db)
    image.rcvgain1 = fread(id, 1, 'short');  % Actual Receive Gain Analog (.1 db)
    image.rcvgain2 = fread(id, 1, 'short');  % Actual Receive Gain Digital (.1 db)
    image.mr_flip = fread(id, 1, 'short');  % Flip Angle for GRASS scans (deg.)
    image.cphase = fread(id, 1, 'short');  % Total Cardiac Phase prescribed
    image.swappf = fread(id, 1, 'short');  % Swap Phase/Frequency Axis
    image.pauseint = fread(id, 1, 'short');  % Pause Interval (slices)
    image.autoxmtgain = fread(id, 1, 'short');  % Auto Transmit Gain (0.1 dB)
    image.prescan_r1 = fread(id, 1, 'short');  % PreScan R1 - Analog
    image.prescan_r2 = fread(id, 1, 'short');  % PreScan R2 - Digital
    image.not_used_1 = fread(id, 1, 'short');  % Available for use
    image.imode = fread(id, 1, 'short');  % Imaging Mode
    image.pseq = fread(id, 1, 'short');  % Pulse Sequence
    image.pseqmode = fread(id, 1, 'short');  % Pulse Sequence Mode
    image.unused_monsar = fread(id, 1, 'short');  % Monitor SAR flag No longer is use
    image.ctyp = fread(id, 1, 'short');  % Coil Type
    image.surfctyp = fread(id, 1, 'short');  % Surface Coil Type
    image.surfcext = fread(id, 1, 'short');  % Extremity Coil Flag
    image.supp_tech = fread(id, 1, 'short');  % SAT fat/water/none
    image.slquant = fread(id, 1, 'short');  % Number of slices in this scan group
    image.gpre = fread(id, 1, 'short');  % Graphically prescribed
    image.satbits = fread(id, 1, 'short');  % Bitmap of SAT selections
    image.scic = fread(id, 1, 'short');  % Surface Coil Intensity Correction Flag
    image.satxloc1 = fread(id, 1, 'short');  % R-side SAT pulse loc rel to lndmrk
    image.satxloc2 = fread(id, 1, 'short');  % L-side SAT pulse loc rel to lndmrk
    image.satyloc1 = fread(id, 1, 'short');  % A-side SAT pulse loc rel to lndmrk
    image.satyloc2 = fread(id, 1, 'short');  % P-side SAT pulse loc rel to lndmrk
    image.satzloc1 = fread(id, 1, 'short');  % S-side SAT pulse loc rel to lndmrk
    image.satzloc2 = fread(id, 1, 'short');  % I-side SAT pulse loc rel to lndmrk
    image.satxthick = fread(id, 1, 'short');  % Thickness of X-axis SAT pulse
    image.satythick = fread(id, 1, 'short');  % Thickness of Y-axis SAT pulse
    image.satzthick = fread(id, 1, 'short');  % Thickness of Z-axis SAT pulse
    image.flax = fread(id, 1, 'short');  % Phase contrast flow axis
    image.venc = fread(id, 1, 'short');  % Phase contrast velocity encoding
    image.thk_disclmr = fread(id, 1, 'short');  % Slice Thickness
    image.ps_flag = fread(id, 1, 'short');  % Auto/Manual Prescan flag
    image.ps_status = fread(id, 1, 'short');  % Bitmap of changed values
    image.image_type = fread(id, 1, 'short');  % Magnitude, Phase, Imaginary, or Real
    image.vas_collapse = fread(id, 1, 'short');  % Collapse Image
    image.proj_alg = fread(id, 1, 'short');  % Projection Algorithm
    image.echo_trn_len = fread(id, 1, 'short');  % Echo Train Length for Fast Spin Echo
    image.frac_echo = fread(id, 1, 'short');  % Fractional Echo - Effective TE Flag
    image.prep_pulse = fread(id, 1, 'short');  % Preporatory Pulse Option
    image.cphasenum = fread(id, 1, 'short');  % Cardiac Phase Number
    image.var_echo = fread(id, 1, 'short');  % Variable Echo Flag
    image.scanactno = fread(id, 1, 'short');  % Scan Acquisition Number
    image.vasflags = fread(id, 1, 'short');  % Magnitude Weighting Flag
    image.integrity = fread(id, 1, 'short');  % GE Image Integrity
    image.freq_dir = fread(id, 1, 'short');  % Frequency Direction
    image.vas_mode = fread(id, 1, 'short');  % Vascular Mode
    image.hole = fread(id, 1, 'short');  %
    image.pscopts = fread(id, 1, 'short');  % bitmap of prescan options
    image.asoffsetx = fread(id, 1, 'short');  % gradient offset in X-direction
    image.asoffsety = fread(id, 1, 'short');  % gradient offset in Y-direction
    image.asoffsetz = fread(id, 1, 'short');  % gradient offset in Z-direction
    image.unoriginal = fread(id, 1, 'short');  % identifies image as original or unoriginal
    image.interleaves = fread(id, 1, 'short');  % number of EPI shots
    image.effechospace = fread(id, 1, 'short');  % effective echo spacing for EPI
    image.viewsperseg = fread(id, 1, 'short');  % views per segment
    image.rbpm = fread(id, 1, 'short');  % respiratory rate, breaths per min
    image.rtpoint = fread(id, 1, 'short');  % respiratory trigger point as percent of max.
    image.rcvrtype = fread(id, 1, 'short');  % type of receiver used
    image.sarMode = fread(id, 1, 'short');  % Sar Ctrl Mode (Normal, 1st or 2nd)
    image.dBdtMode = fread(id, 1, 'short');  % dBdt Ctrl Mode (Normal, 1st or 2nd)
    image.govBody = fread(id, 1, 'short');  % Governing Body MHW/IEC/FDA
    image.sarDefinition = fread(id, 1, 'short');  % Defined terms avaialble
    image.no_shimvol = fread(id, 1, 'short');  %
    image.shim_vol_type = fread(id, 1, 'short');  %
    image.psdname = modchar(fread(id, 33, '*char'));  % Pulse Sequence Name
    image.proj_name = modchar(fread(id, 13, '*char'));  % Projection Algorithm Name
    image.psd_iname = modchar(fread(id, 13, '*char'));  % PSD name from inside PSD
    image.im_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Image
    image.pdid = modchar(fread(id, 14, '*char'));  % Pixel Data ID
    image.im_suid = modchar(fread(id, 4, '*char'));  % Suite id for this image
    image.contrastIV = modchar(fread(id, 17, '*char'));  % IV Contrast Agent
    image.contrastOral = modchar(fread(id, 17, '*char'));  % Oral Contrast Agent
    image.loc_ras = modchar(fread(id, 1, '*char'));  % RAS letter of image location
    image.forimgrev = modchar(fread(id, 4, '*char'));  % Foreign Image Revision
    image.cname = modchar(fread(id, 17, '*char'));  % Coil Name
    image.im_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    image.im_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    image.im_alloc_key = modchar(fread(id, 13, '*char'));  %
    image.ref_img = modchar(fread(id, 1, '*char'));  % Reference Image Field
    image.sum_img = modchar(fread(id, 1, '*char'));  % Summary Image Field
    image.filter_mode = modchar(fread(id, 16, '*char'));  % String Slop Field 1
    image.slop_str_2 = modchar(fread(id, 16, '*char'));  % String Slop Field 2
    image.image_uid = moduid(fread(id, 32, 'uint8'));  % Image Unique ID
    image.sop_uid = moduid(fread(id, 32, 'uint8'));  % Service Obj Class Unique ID
    image.GEcname = modchar(fread(id, 24, '*char'));  % GECoilname for the cname
    image.usedCoilData = modchar(fread(id, 100, '*char'));  % Concatenated str of coilcode and chip serialID
    image.astcalseriesuid = moduid(fread(id, 32, 'uint8'));  %
    image.purecalseriesuid = moduid(fread(id, 32, 'uint8'));  %
    image.sys_type = modchar(fread(id, 64, '*char'));  %
    image.xml_psc_shm_vol = modchar(fread(id, 32, '*char'));  %
    image.img_hdr_padding = modchar(fread(id, 164, '*char'));  %
   
elseif (ver >= 14.0 & ver < 20)
    
    image.double_padding = fread(id, 64, 'char'); % double padding % by jz
    
    image.dfov = fread(id, 1, 'float');  % Display field of view - X (mm)
    image.dfov_rect = fread(id, 1, 'float');  % Display field of view - Y (if different)
    image.sctime = fread(id, 1, 'float');  % Duration of scan
    image.slthick = fread(id, 1, 'float');  % Slice Thickness (mm)
    image.scanspacing = fread(id, 1, 'float');  % Spacing between scans (mm?)
    image.loc = fread(id, 1, 'float');  % Image location
    image.tbldlta = fread(id, 1, 'float');  % Table Delta
    image.nex = fread(id, 1, 'float');  % Number of Excitations
    image.reptime = fread(id, 1, 'float');  % Cardiac repetition time
    image.saravg = fread(id, 1, 'float');  % Average SAR
    image.sarpeak = fread(id, 1, 'float');  % Peak SAR
    image.pausetime = fread(id, 1, 'float');  % Pause Time
    image.vbw = fread(id, 1, 'float');  % Variable Bandwidth (Hz)
    image.user0 = fread(id, 1, 'float');  % User Variable 0
    image.user1 = fread(id, 1, 'float');  % User Variable 1
    image.user2 = fread(id, 1, 'float');  % User Variable 2
    image.user3 = fread(id, 1, 'float');  % User Variable 3
    image.user4 = fread(id, 1, 'float');  % User Variable 4
    image.user5 = fread(id, 1, 'float');  % User Variable 5
    image.user6 = fread(id, 1, 'float');  % User Variable 6
    image.user7 = fread(id, 1, 'float');  % User Variable 7
    image.user8 = fread(id, 1, 'float');  % User Variable 8
    image.user9 = fread(id, 1, 'float');  % User Variable 9
    image.user10 = fread(id, 1, 'float');  % User Variable 10
    image.user11 = fread(id, 1, 'float');  % User Variable 11
    image.user12 = fread(id, 1, 'float');  % User Variable 12
    image.user13 = fread(id, 1, 'float');  % User Variable 13
    image.user14 = fread(id, 1, 'float');  % User Variable 14
    image.user15 = fread(id, 1, 'float');  % User Variable 15
    image.user16 = fread(id, 1, 'float');  % User Variable 16
    image.user17 = fread(id, 1, 'float');  % User Variable 17
    image.user18 = fread(id, 1, 'float');  % User Variable 18
    image.user19 = fread(id, 1, 'float');  % User Variable 19
    image.user20 = fread(id, 1, 'float');  % User Variable 20
    image.user21 = fread(id, 1, 'float');  % User Variable 21
    image.user22 = fread(id, 1, 'float');  % User Variable 22
    image.proj_ang = fread(id, 1, 'float');  % Projection Angle
    image.concat_sat = fread(id, 1, 'float');  % Concat Sat Type Flag
    image.user23 = fread(id, 1, 'float');  % User Variable 23
    image.user24 = fread(id, 1, 'float');  % User Variable 24
    image.x_axis_rot = fread(id, 1, 'float');  % X Axis Rotation
    image.y_axis_rot = fread(id, 1, 'float');  % Y Axis Rotation
    image.z_axis_rot = fread(id, 1, 'float');  % Z Axis Rotation
    image.ihtagfa = fread(id, 1, 'float');  % Tagging Flip Angle
    image.ihtagor = fread(id, 1, 'float');  % Cardiac Tagging Orientation
    image.ihbspti = fread(id, 1, 'float');  % Blood Suppression TI
    image.rtia_timer = fread(id, 1, 'float');  % Float Slop Field 4
    image.fps = fread(id, 1, 'float');  % Float Slop Field 5
    image.vencscale = fread(id, 1, 'float');  % Scale Weighted Venc
    image.dbdt = fread(id, 1, 'float');  % peak rate of change of gradient field, tesla/sec
    image.dbdtper = fread(id, 1, 'float');  % limit in units of percent of theoretical curve
    image.estdbdtper = fread(id, 1, 'float');  % PSD estimated limit in units of percent
    image.estdbdtts = fread(id, 1, 'float');  % PSD estimated limit in Teslas/sec
    image.saravghead = fread(id, 1, 'float');  % Avg head SAR
    image.neg_scanspacing = fread(id, 1, 'float');  % Negative scan spacing for overlap slices
    image.user25 = fread(id, 1, 'float');  % User Variable 25
    image.user26 = fread(id, 1, 'float');  % User Variable 26
    image.user27 = fread(id, 1, 'float');  % User Variable 27
    image.user28 = fread(id, 1, 'float');  % User Variable 28
    image.user29 = fread(id, 1, 'float');  % User Variable 29
    image.user30 = fread(id, 1, 'float');  % User Variable 30
    image.user31 = fread(id, 1, 'float');  % User Variable 31
    image.user32 = fread(id, 1, 'float');  % User Variable 32
    image.user33 = fread(id, 1, 'float');  % User Variable 33
    image.user34 = fread(id, 1, 'float');  % User Variable 34
    image.user35 = fread(id, 1, 'float');  % User Variable 35
    image.user36 = fread(id, 1, 'float');  % User Variable 36
    image.user37 = fread(id, 1, 'float');  % User Variable 37
    image.user38 = fread(id, 1, 'float');  % User Variable 38
    image.user39 = fread(id, 1, 'float');  % User Variable 39
    image.user40 = fread(id, 1, 'float');  % User Variable 40
    image.user41 = fread(id, 1, 'float');  % User Variable 41
    image.user42 = fread(id, 1, 'float');  % User Variable 42
    image.user43 = fread(id, 1, 'float');  % User Variable 43
    image.user44 = fread(id, 1, 'float');  % User Variable 44
    image.user45 = fread(id, 1, 'float');  % User Variable 45
    image.user46 = fread(id, 1, 'float');  % User Variable 46
    image.user47 = fread(id, 1, 'float');  % User Variable 47
    image.CAI_eff_res = fread(id, 1, 'float');  % 
    image.RegressorVal = fread(id, 1, 'float');  %
    image.SliceAsset = fread(id, 1, 'float');  % Slice Asset in Asset Screen
    image.PhaseAsset = fread(id, 1, 'float');  % Phase Asset in Asset Screen
    image.sarValues = fread(id, 4, 'float');  % correspoding SAR values for defined terms
    image.shim_fov = fread(id, 2, 'float');  %
    image.shim_ctr_R = fread(id, 2, 'float');  %
    image.shim_ctr_A = fread(id, 2, 'float');  %
    image.shim_ctr_S = fread(id, 2, 'float');  %
    image.dim_X = fread(id, 1, 'float');  % Image dimension - X
    image.dim_Y = fread(id, 1, 'float');  % Image dimension - Y
    image.pixsize_X = fread(id, 1, 'float');  % Image pixel size - X
    image.pixsize_Y = fread(id, 1, 'float');  % Image pixel size - Y
    image.ctr_R = fread(id, 1, 'float');  % Center R coord of plane image
    image.ctr_A = fread(id, 1, 'float');  % Center A coord of plane image
    image.ctr_S = fread(id, 1, 'float');  % Center S coord of plane image
    image.norm_R = fread(id, 1, 'float');  % Normal R coord
    image.norm_A = fread(id, 1, 'float');  % Normal A coord
    image.norm_S = fread(id, 1, 'float');  % Normal S coord
    image.tlhc_R = fread(id, 1, 'float');  % R Coord of Top Left Hand Corner
    image.tlhc_A = fread(id, 1, 'float');  % A Coord of Top Left Hand Corner
    image.tlhc_S = fread(id, 1, 'float');  % S Coord of Top Left Hand Corner
    image.trhc_R = fread(id, 1, 'float');  % R Coord of Top Right Hand Corner
    image.trhc_A = fread(id, 1, 'float');  % A Coord of Top Right Hand Corner
    image.trhc_S = fread(id, 1, 'float');  % S Coord of Top Right Hand Corner
    image.brhc_R = fread(id, 1, 'float');  % R Coord of Bottom Right Hand Corner
    image.brhc_A = fread(id, 1, 'float');  % A Coord of Bottom Right Hand Corner
    image.brhc_S = fread(id, 1, 'float');  % S Coord of Bottom Right Hand Corner
    
    image.float_padding = fread(id, 128, 'char'); % float padding
    
    image.cal_fldstr = fread(id, 1, 'long');  % Calibrated Field Strength (x10 uGauss)
    image.checksum = fread(id, 1, 'long'); 
    image.user_usage_tag = fread(id, 1, 'long');  % Defines how following user CVs are to be filled in
    
    image.user_fill_mapMSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user_fill_mapLSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.im_archived = fread(id, 1, 'long');  % Image Archive Flag
    image.im_complete = fread(id, 1, 'long');  % Image Complete Flag
    
    image.long_padding = fread(id, 32, 'char'); % long padding
    
    image.im_datetime = fread(id, 1, 'int');  % Allocation Image date/time stamp
    image.im_actual_dt = fread(id, 1, 'int');  % Actual Image date/time stamp
    image.tr = fread(id, 1, 'int');  % Pulse repetition time(usec)
    image.ti = fread(id, 1, 'int');  % Pulse inversion time(usec)
    image.te = fread(id, 1, 'int');  % Pulse echo time(usec)
    image.te2 = fread(id, 1, 'int');  % Second echo echo (usec)
    image.tdel = fread(id, 1, 'int');  % Delay time after trigger (msec)
    image.mindat = fread(id, 1, 'int');  % Minimum Delay after Trigger (uSec)
    image.obplane = fread(id, 1, 'int');  % Oblique Plane
    image.slocfov = fread(id, 1, 'int');  % Slice Offsets on Freq axis
    
    image.xmtfreq = fread(id, 1, 'int');  % Center Frequency (0.1 Hz)
    image.autoxmtfreq = fread(id, 1, 'int');  % Auto Center Frequency (0.1 Hz)
    
    image.user_bitmap = fread(id, 1, 'int');  % Bitmap defining user CVs
    image.iopt = fread(id, 1, 'int');  % Imaging Options
    image.psd_datetime = fread(id, 1, 'int');  % PSD Creation Date and Time
    image.rawrunnum = fread(id, 1, 'int');  % RawData Run Number
    image.intr_del = fread(id, 1, 'int');  % Interimage/interloc delay (uSec)
    image.im_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    image.im_pds_a = fread(id, 1, 'int');  % PixelData size - as stored
    image.im_pds_c = fread(id, 1, 'int');  % PixelData size - Compressed
    image.im_pds_u = fread(id, 1, 'int');  % PixelData size - UnCompressed
    image.thresh_min1 = fread(id, 1, 'int');  % Lower Range of Pixels 1
    image.thresh_max1 = fread(id, 1, 'int');  % Upper Range of Pixels 1
    image.thresh_min2 = fread(id, 1, 'int');  % Lower Range of Pixels 2
    image.thresh_max2 = fread(id, 1, 'int');  % Upper Range of Pixels 2
    image.numslabs = fread(id, 1, 'int');  % Number of 3D Slabs
    image.locsperslab = fread(id, 1, 'int');  % Slice Locs Per 3D Slab
    image.overlaps = fread(id, 1, 'int');  % # of Slice Locs on Each Slab Which Overlap N eighbors
    image.slop_int_4 = fread(id, 1, 'int');  % Image Filtering 0.5/0.2T
    image.dfax = fread(id, 1, 'int');  % Diffusion Direction for DW-EPI
    image.fphase = fread(id, 1, 'int');  % Number Of Phases
    image.offsetfreq = fread(id, 1, 'int');  % Offset Frequency - Mag.Transfer
    image.b_value = fread(id, 1, 'int');  % B-value for DW-EPI
    image.iopt2 = fread(id, 1, 'int');  % Imaging Option2
    image.ihtagging = fread(id, 1, 'int');  % tag type
    image.ihtagspc = fread(id, 1, 'int');  % tag space
    image.ihfcineim = fread(id, 1, 'int');  % Fast CINE interpolation method
    image.ihfcinent = fread(id, 1, 'int');  % Fast CINE normalization type
    image.num_seg = fread(id, 1, 'int');  % YMSge05074
    image.oprtarr = fread(id, 1, 'int');  % Respiratory Trigger windo
    image.averages = fread(id, 1, 'int');  % Number of averages for spectro
    image.station_index = fread(id, 1, 'int');  % Station Index
    image.station_total = fread(id, 1, 'int');  % Station Total
    image.iopt3 = fread(id, 1, 'int');  % Imaging Option3
    image.delAcq = fread(id, 1, 'int');  % Delay after Acquisition (MP / fMRI screen)
    
    image.rxmbloblen = fread(id, 1, 'int');  % rxmbloblen
    image.rxmblob = fread(id, 1, 'int');  % rxmblob_pad    
    image.int_padding = fread(id, 128, 'char'); % int padding
    
    image.imatrix_X = fread(id, 1, 'short');  % Image matrix size - X
    image.imatrix_Y = fread(id, 1, 'short');  % Image matrix size - Y
    image.im_exno = fread(id, 1, 'ushort');  % Exam number for this image
    image.img_window = fread(id, 1, 'ushort');  % Window Value
    image.img_level = fread(id, 1, 'short');  % Level Value
    image.numecho = fread(id, 1, 'short');  % Number of echoes
    image.echonum = fread(id, 1, 'short');  % Echo Number
    image.im_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    image.im_seno = fread(id, 1, 'short');  % Series Number for this image
    image.im_no = fread(id, 1, 'short');  % Image Number
    image.contmode = fread(id, 1, 'short');  % Image Contrast Mode
    image.serrx = fread(id, 1, 'short');  % Series from which prescribed
    image.imgrx = fread(id, 1, 'short');  % Image from which prescribed
    image.screenformat = fread(id, 1, 'short');  % Screen Format(8/16 bit)
    image.plane = fread(id, 1, 'short');  % Plane Type
    image.im_compress = fread(id, 1, 'short');  % Image compression type for allocation
    image.im_scouttype = fread(id, 1, 'short');  % Scout Type (AP or lateral)
    image.contig = fread(id, 1, 'short');  % Continuous Slices Flag
    image.hrtrate = fread(id, 1, 'short');  % Cardiac Heart Rate (bpm)
    image.trgwindow = fread(id, 1, 'short');  % Trigger window (% of R-R interval)
    image.imgpcyc = fread(id, 1, 'short');  % Images per cardiac cycle
    
    image.xmtgain = fread(id, 1, 'short');  % Actual Transmit Gain (.1 db)
    image.rcvgain1 = fread(id, 1, 'short');  % Actual Receive Gain Analog (.1 db)
    image.rcvgain2 = fread(id, 1, 'short');  % Actual Receive Gain Digital (.1 db)
    
    image.mr_flip = fread(id, 1, 'short');  % Flip Angle for GRASS scans (deg.)
    image.cphase = fread(id, 1, 'short');  % Total Cardiac Phase prescribed
    image.swappf = fread(id, 1, 'short');  % Swap Phase/Frequency Axis
    image.pauseint = fread(id, 1, 'short');  % Pause Interval (slices)
    
    image.autoxmtgain = fread(id, 1, 'short');  % Auto Transmit Gain (0.1 dB)
    image.prescan_r1 = fread(id, 1, 'short');  % PreScan R1 - Analog
    image.prescan_r2 = fread(id, 1, 'short');  % PreScan R2 - Digital
    
    image.not_used_1 = fread(id, 1, 'short');  % Available for use
    image.imode = fread(id, 1, 'short');  % Imaging Mode
    image.pseq = fread(id, 1, 'short');  % Pulse Sequence
    image.pseqmode = fread(id, 1, 'short');  % Pulse Sequence Mode
    image.unused_monsar = fread(id, 1, 'short');
    image.ctyp = fread(id, 1, 'short');  % Coil Type
    image.surfctyp = fread(id, 1, 'short');  % Surface Coil Type
    image.surfcext = fread(id, 1, 'short');  % Extremity Coil Flag
    image.supp_tech = fread(id, 1, 'short');  % SAT fat/water/none
    image.slquant = fread(id, 1, 'short');  % Number of slices in this scan group
    image.gpre = fread(id, 1, 'short');  % Graphically prescribed
    image.satbits = fread(id, 1, 'short');  % Bitmap of SAT selections
    image.scic = fread(id, 1, 'short');  % Surface Coil Intensity Correction Flag
    image.satxloc1 = fread(id, 1, 'short');  % R-side SAT pulse loc rel to lndmrk
    image.satxloc2 = fread(id, 1, 'short');  % L-side SAT pulse loc rel to lndmrk
    image.satyloc1 = fread(id, 1, 'short');  % A-side SAT pulse loc rel to lndmrk
    image.satyloc2 = fread(id, 1, 'short');  % P-side SAT pulse loc rel to lndmrk
    image.satzloc1 = fread(id, 1, 'short');  % S-side SAT pulse loc rel to lndmrk
    image.satzloc2 = fread(id, 1, 'short');  % I-side SAT pulse loc rel to lndmrk
    image.satxthick = fread(id, 1, 'short');  % Thickness of X-axis SAT pulse
    image.satythick = fread(id, 1, 'short');  % Thickness of Y-axis SAT pulse
    image.satzthick = fread(id, 1, 'short');  % Thickness of Z-axis SAT pulse
    image.flax = fread(id, 1, 'short');  % Phase contrast flow axis
    image.venc = fread(id, 1, 'short');  % Phase contrast velocity encoding
    image.thk_disclmr = fread(id, 1, 'short');  % Slice Thickness
    
    image.ps_flag = fread(id, 1, 'short');  % Auto/Manual Prescan flag
    image.ps_status = fread(id, 1, 'short');  % Bitmap of changed values
    
    image.image_type = fread(id, 1, 'short');  % Magnitude, Phase, Imaginary, or Real
    image.vas_collapse = fread(id, 1, 'short');  % Collapse Image
    image.proj_alg = fread(id, 1, 'short');  % Projection Algorithm
    image.echo_trn_len = fread(id, 1, 'short');  % Echo Train Length for Fast Spin Echo
    image.frac_echo = fread(id, 1, 'short');  % Fractional Echo - Effective TE Flag
    image.prep_pulse = fread(id, 1, 'short');  % Preporatory Pulse Option
    image.cphasenum = fread(id, 1, 'short');  % Cardiac Phase Number
    image.var_echo = fread(id, 1, 'short');  % Variable Echo Flag
    image.scanactno = fread(id, 1, 'short');  % Scan Acquisition Number
    image.vasflags = fread(id, 1, 'short');  % Magnitude Weighting Flag
    image.integrity = fread(id, 1, 'short');  % GE Image Integrity
    image.freq_dir = fread(id, 1, 'short');  % Frequency Direction
    image.vas_mode = fread(id, 1, 'short');  % Vascular Mode
    image.hole = fread(id, 1, 'short');  % Hole
    image.pscopts = fread(id, 1, 'short');  % bitmap of prescan options
    
    image.asoffsetx = fread(id, 1, 'short');  % gradient offset in X-direction
    image.asoffsety = fread(id, 1, 'short');  % gradient offset in Y-direction
    image.asoffsetz = fread(id, 1, 'short');  % gradient offset in Z-direction
    
    image.unoriginal = fread(id, 1, 'short');  % identifies image as original or unoriginal
    image.interleaves = fread(id, 1, 'short');  % number of EPI shots
    image.effechospace = fread(id, 1, 'short');  % effective echo spacing for EPI
    image.viewsperseg = fread(id, 1, 'short');  % views per segment
    image.rbpm = fread(id, 1, 'short');  % respiratory rate, breaths per min
    image.rtpoint = fread(id, 1, 'short');  % respiratory trigger point as percent of max.
    image.rcvrtype = fread(id, 1, 'short');  % type of receiver used
    image.sarMode = fread(id, 1, 'short');  % Sar Ctrl Mode (Normal, 1st or 2nd)
    image.dBdtMode = fread(id, 1, 'short');  % dBdt Ctrl Mode (Normal, 1st or 2nd)
    image.govBody = fread(id, 1, 'short');  % Governing Body MHW/IEC/FDA
    image.sarDefinition = fread(id, 1, 'short');  % Defined terms avaialble
    image.no_shimvol = fread(id, 1, 'short');  %
    image.shim_vol_type = fread(id, 1, 'short');  %
    
    image.current_phase = fread(id, 1, 'short');  % Current phase
    image.art_level = fread(id, 1, 'short');  % Art level
    image.slice_group_numer = fread(id, 1, 'short');  % Slice group number
    image.number_of_slice_groups = fread(id, 1, 'short');  % # of slc groups
    image.show_in_autoview = fread(id, 1, 'short');  %    
    image.short_padding = fread(id, 64, 'char'); % short padding
    
    image.psdname = modchar(fread(id, 33, '*char'));  % Pulse Sequence Name
    
    % fseek(id,pos_image+1665,-1);  % there is a bug here if no fseek, why?
    
    image.proj_name = modchar(fread(id, 13, '*char'));  % Projection Algorithm Name
    image.psd_iname = modchar(fread(id, 13, '*char'));  % PSD name from inside PSD
    image.im_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Image
    image.pdid = modchar(fread(id, 14, '*char'));  % Pixel Data ID
    image.im_suid = modchar(fread(id, 4, '*char'));  % Suite id for this image
    image.contrastIV = modchar(fread(id, 17, '*char'));  % IV Contrast Agent
    image.contrastOral = modchar(fread(id, 17, '*char'));  % Oral Contrast Agent
    image.loc_ras = modchar(fread(id, 1, '*char'));  % RAS letter of image location
    image.forimgrev = modchar(fread(id, 4, '*char'));  % Foreign Image Revision
    image.cname = modchar(fread(id, 17, '*char'));  % Coil Name
    image.im_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    image.im_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    image.im_alloc_key = modchar(fread(id, 13, '*char'));  %
    image.ref_img = modchar(fread(id, 1, '*char'));  % Reference Image Field
    image.sum_img = modchar(fread(id, 1, '*char'));  % Summary Image Field
    image.filter_mode = modchar(fread(id, 16, '*char'));  % String Slop Field 1
    image.slop_str_2 = modchar(fread(id, 16, '*char'));  % String Slop Field 2
    image.image_uid = moduid(fread(id, 32, 'uint8'));  % Image Unique ID
    image.sop_uid = moduid(fread(id, 32, 'uint8'));  % Service Obj Class Unique ID
    image.GEcname = modchar(fread(id, 24, '*char'));  % GECoilname for the cname
    image.usedCoilData = modchar(fread(id, 100, '*char'));  % Concatenated str of coilcode and chip serialID
    image.astcalseriesuid = moduid(fread(id, 32, 'uint8'));  %
    image.purecalseriesuid = moduid(fread(id, 32, 'uint8'));  %
    
    image.xml_psc_shm_vol = moduid(fread(id, 32, 'char'));  %
    image.rxmpath = moduid(fread(id, 64, 'char'));  %
    image.psdnameannot = moduid(fread(id, 33, 'char'));  %    
    image.img_hdr_padding = modchar(fread(id, 210, '*char'));  %
    
    
elseif (ver >= 20.0)
    
    image.autoSubParam0 = fread(id, 1, 'double');  % autoSubParam0
    image.autoSubParam1 = fread(id, 1, 'double');  % autoSubParam1
    image.autoSubParam2 = fread(id, 1, 'double');  % autoSubParam2
    image.autoSubParam3 = fread(id, 1, 'double');  % autoSubParam3
    image.autoSubParam4 = fread(id, 1, 'double');  % autoSubParam4
    image.autoSubParam5 = fread(id, 1, 'double');  % autoSubParam5
    
    image.double_padding = fread(id, 256, 'char'); % double padding
    
    image.dfov = fread(id, 1, 'float');  % Display field of view - X (mm)
    image.dfov_rect = fread(id, 1, 'float');  % Display field of view - Y (if different)
    image.sctime = fread(id, 1, 'float');  % Duration of scan
    image.slthick = fread(id, 1, 'float');  % Slice Thickness (mm)
    image.scanspacing = fread(id, 1, 'float');  % Spacing between scans (mm?)
    image.loc = fread(id, 1, 'float');  % Image location
    image.tbldlta = fread(id, 1, 'float');  % Table Delta
    image.nex = fread(id, 1, 'float');  % Number of Excitations
    image.reptime = fread(id, 1, 'float');  % Cardiac repetition time
    image.saravg = fread(id, 1, 'float');  % Average SAR
    image.sarpeak = fread(id, 1, 'float');  % Peak SAR
    image.pausetime = fread(id, 1, 'float');  % Pause Time
    image.vbw = fread(id, 1, 'float');  % Variable Bandwidth (Hz)
    image.user0 = fread(id, 1, 'float');  % User Variable 0
    image.user1 = fread(id, 1, 'float');  % User Variable 1
    image.user2 = fread(id, 1, 'float');  % User Variable 2
    image.user3 = fread(id, 1, 'float');  % User Variable 3
    image.user4 = fread(id, 1, 'float');  % User Variable 4
    image.user5 = fread(id, 1, 'float');  % User Variable 5
    image.user6 = fread(id, 1, 'float');  % User Variable 6
    image.user7 = fread(id, 1, 'float');  % User Variable 7
    image.user8 = fread(id, 1, 'float');  % User Variable 8
    image.user9 = fread(id, 1, 'float');  % User Variable 9
    image.user10 = fread(id, 1, 'float');  % User Variable 10
    image.user11 = fread(id, 1, 'float');  % User Variable 11
    image.user12 = fread(id, 1, 'float');  % User Variable 12
    image.user13 = fread(id, 1, 'float');  % User Variable 13
    image.user14 = fread(id, 1, 'float');  % User Variable 14
    image.user15 = fread(id, 1, 'float');  % User Variable 15
    image.user16 = fread(id, 1, 'float');  % User Variable 16
    image.user17 = fread(id, 1, 'float');  % User Variable 17
    image.user18 = fread(id, 1, 'float');  % User Variable 18
    image.user19 = fread(id, 1, 'float');  % User Variable 19
    image.user20 = fread(id, 1, 'float');  % User Variable 20
    image.user21 = fread(id, 1, 'float');  % User Variable 21
    image.user22 = fread(id, 1, 'float');  % User Variable 22
    image.proj_ang = fread(id, 1, 'float');  % Projection Angle
    image.concat_sat = fread(id, 1, 'float');  % Concat Sat Type Flag
    image.user23 = fread(id, 1, 'float');  % User Variable 23
    image.user24 = fread(id, 1, 'float');  % User Variable 24
    image.x_axis_rot = fread(id, 1, 'float');  % X Axis Rotation
    image.y_axis_rot = fread(id, 1, 'float');  % Y Axis Rotation
    image.z_axis_rot = fread(id, 1, 'float');  % Z Axis Rotation
    image.ihtagfa = fread(id, 1, 'float');  % Tagging Flip Angle
    image.ihtagor = fread(id, 1, 'float');  % Cardiac Tagging Orientation
    image.ihbspti = fread(id, 1, 'float');  % Blood Suppression TI
    image.rtia_timer = fread(id, 1, 'float');  % Float Slop Field 4
    image.fps = fread(id, 1, 'float');  % Float Slop Field 5
    image.vencscale = fread(id, 1, 'float');  % Scale Weighted Venc
    image.dbdt = fread(id, 1, 'float');  % peak rate of change of gradient field, tesla/sec
    image.dbdtper = fread(id, 1, 'float');  % limit in units of percent of theoretical curve
    image.estdbdtper = fread(id, 1, 'float');  % PSD estimated limit in units of percent
    image.estdbdtts = fread(id, 1, 'float');  % PSD estimated limit in Teslas/sec
    image.saravghead = fread(id, 1, 'float');  % Avg head SAR
    image.neg_scanspacing = fread(id, 1, 'float');  % Negative scan spacing for overlap slices
    image.user25 = fread(id, 1, 'float');  % User Variable 25
    image.user26 = fread(id, 1, 'float');  % User Variable 26
    image.user27 = fread(id, 1, 'float');  % User Variable 27
    image.user28 = fread(id, 1, 'float');  % User Variable 28
    image.user29 = fread(id, 1, 'float');  % User Variable 29
    image.user30 = fread(id, 1, 'float');  % User Variable 30
    image.user31 = fread(id, 1, 'float');  % User Variable 31
    image.user32 = fread(id, 1, 'float');  % User Variable 32
    image.user33 = fread(id, 1, 'float');  % User Variable 33
    image.user34 = fread(id, 1, 'float');  % User Variable 34
    image.user35 = fread(id, 1, 'float');  % User Variable 35
    image.user36 = fread(id, 1, 'float');  % User Variable 36
    image.user37 = fread(id, 1, 'float');  % User Variable 37
    image.user38 = fread(id, 1, 'float');  % User Variable 38
    image.user39 = fread(id, 1, 'float');  % User Variable 39
    image.user40 = fread(id, 1, 'float');  % User Variable 40
    image.user41 = fread(id, 1, 'float');  % User Variable 41
    image.user42 = fread(id, 1, 'float');  % User Variable 42
    image.user43 = fread(id, 1, 'float');  % User Variable 43
    image.user44 = fread(id, 1, 'float');  % User Variable 44
    image.user45 = fread(id, 1, 'float');  % User Variable 45
    image.user46 = fread(id, 1, 'float');  % User Variable 46
    image.user47 = fread(id, 1, 'float');  % User Variable 47
    image.CAI_eff_res = fread(id, 1, 'float');  % 
    image.RegressorVal = fread(id, 1, 'float');  %
    image.SliceAsset = fread(id, 1, 'float');  % Slice Asset in Asset Screen
    image.PhaseAsset = fread(id, 1, 'float');  % Phase Asset in Asset Screen
    image.sarValues = fread(id, 4, 'float');  % correspoding SAR values for defined terms
    image.shim_fov = fread(id, 2, 'float');  %
    image.shim_ctr_R = fread(id, 2, 'float');  %
    image.shim_ctr_A = fread(id, 2, 'float');  %
    image.shim_ctr_S = fread(id, 2, 'float');  %
    image.dim_X = fread(id, 1, 'float');  % Image dimension - X
    image.dim_Y = fread(id, 1, 'float');  % Image dimension - Y
    image.pixsize_X = fread(id, 1, 'float');  % Image pixel size - X
    image.pixsize_Y = fread(id, 1, 'float');  % Image pixel size - Y
    image.ctr_R = fread(id, 1, 'float');  % Center R coord of plane image
    image.ctr_A = fread(id, 1, 'float');  % Center A coord of plane image
    image.ctr_S = fread(id, 1, 'float');  % Center S coord of plane image
    image.norm_R = fread(id, 1, 'float');  % Normal R coord
    image.norm_A = fread(id, 1, 'float');  % Normal A coord
    image.norm_S = fread(id, 1, 'float');  % Normal S coord
    image.tlhc_R = fread(id, 1, 'float');  % R Coord of Top Left Hand Corner
    image.tlhc_A = fread(id, 1, 'float');  % A Coord of Top Left Hand Corner
    image.tlhc_S = fread(id, 1, 'float');  % S Coord of Top Left Hand Corner
    image.trhc_R = fread(id, 1, 'float');  % R Coord of Top Right Hand Corner
    image.trhc_A = fread(id, 1, 'float');  % A Coord of Top Right Hand Corner
    image.trhc_S = fread(id, 1, 'float');  % S Coord of Top Right Hand Corner
    image.brhc_R = fread(id, 1, 'float');  % R Coord of Bottom Right Hand Corner
    image.brhc_A = fread(id, 1, 'float');  % A Coord of Bottom Right Hand Corner
    image.brhc_S = fread(id, 1, 'float');  % S Coord of Bottom Right Hand Corner
    
    image.float_padding = fread(id, 132, 'char'); % float padding
    
    image.cal_fldstr = fread(id, 1, 'int');  % Calibrated Field Strength (x10 uGauss)
    image.user_usage_tag = fread(id, 1, 'int');  % Defines how following user CVs are to be filled in
    
    image.user_fill_mapMSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.user_fill_mapLSW = fread(id, 1, 'ulong');  % Define what process fills in the user CVs, ifcc or TIR
    image.im_archived = fread(id, 1, 'long');  % Image Archive Flag
    image.im_complete = fread(id, 1, 'long');  % Image Complete Flag
    
    image.int_padding1 = fread(id, 136, 'char'); % int padding 1
    
    image.im_datetime = fread(id, 1, 'int');  % Allocation Image date/time stamp
    image.im_actual_dt = fread(id, 1, 'int');  % Actual Image date/time stamp
    image.tr = fread(id, 1, 'int');  % Pulse repetition time(usec)
    image.ti = fread(id, 1, 'int');  % Pulse inversion time(usec)
    image.te = fread(id, 1, 'int');  % Pulse echo time(usec)
    image.te2 = fread(id, 1, 'int');  % Second echo echo (usec)
    image.tdel = fread(id, 1, 'int');  % Delay time after trigger (msec)
    image.mindat = fread(id, 1, 'int');  % Minimum Delay after Trigger (uSec)
    image.obplane = fread(id, 1, 'int');  % Oblique Plane
    image.slocfov = fread(id, 1, 'int');  % Slice Offsets on Freq axis
    
    image.obsolote1 = fread(id, 1, 'int');  % Center Frequency (0.1 Hz)
    image.obsolete2 = fread(id, 1, 'int');  % Auto Center Frequency (0.1 Hz)
    
    image.user_bitmap = fread(id, 1, 'int');  % Bitmap defining user CVs
    image.iopt = fread(id, 1, 'int');  % Imaging Options
    image.psd_datetime = fread(id, 1, 'int');  % PSD Creation Date and Time
    image.rawrunnum = fread(id, 1, 'int');  % RawData Run Number
    image.intr_del = fread(id, 1, 'int');  % Interimage/interloc delay (uSec)
    image.im_lastmod = fread(id, 1, 'int');  % Date/Time of Last Change
    image.im_pds_a = fread(id, 1, 'int');  % PixelData size - as stored
    image.im_pds_c = fread(id, 1, 'int');  % PixelData size - Compressed
    image.im_pds_u = fread(id, 1, 'int');  % PixelData size - UnCompressed
    image.thresh_min1 = fread(id, 1, 'int');  % Lower Range of Pixels 1
    image.thresh_max1 = fread(id, 1, 'int');  % Upper Range of Pixels 1
    image.thresh_min2 = fread(id, 1, 'int');  % Lower Range of Pixels 2
    image.thresh_max2 = fread(id, 1, 'int');  % Upper Range of Pixels 2
    image.numslabs = fread(id, 1, 'int');  % Number of 3D Slabs
    image.locsperslab = fread(id, 1, 'int');  % Slice Locs Per 3D Slab
    image.overlaps = fread(id, 1, 'int');  % # of Slice Locs on Each Slab Which Overlap N eighbors
    image.slop_int_4 = fread(id, 1, 'int');  % Image Filtering 0.5/0.2T
    image.dfax = fread(id, 1, 'int');  % Diffusion Direction for DW-EPI
    image.fphase = fread(id, 1, 'int');  % Number Of Phases
    image.offsetfreq = fread(id, 1, 'int');  % Offset Frequency - Mag.Transfer
    image.b_value = fread(id, 1, 'int');  % B-value for DW-EPI
    image.iopt2 = fread(id, 1, 'int');  % Imaging Option2
    image.ihtagging = fread(id, 1, 'int');  % tag type
    image.ihtagspc = fread(id, 1, 'int');  % tag space
    image.ihfcineim = fread(id, 1, 'int');  % Fast CINE interpolation method
    image.ihfcinent = fread(id, 1, 'int');  % Fast CINE normalization type
    image.num_seg = fread(id, 1, 'int');  % YMSge05074
    image.oprtarr = fread(id, 1, 'int');  % Respiratory Trigger windo
    image.averages = fread(id, 1, 'int');  % Number of averages for spectro
    image.station_index = fread(id, 1, 'int');  % Station Index
    image.station_total = fread(id, 1, 'int');  % Station Total
    image.iopt3 = fread(id, 1, 'int');  % Imaging Option3
    image.delAcq = fread(id, 1, 'int');  % Delay after Acquisition (MP / fMRI screen)
    
    image.rxmbloblen = fread(id, 1, 'int');  % rxmbloblen
    image.rxmblob_pad = fread(id, 1, 'int');  % rxmblob_pad
    
    image.int_padding2 = fread(id, 132, 'char'); % int padding 2
    
    image.imatrix_X = fread(id, 1, 'short');  % Image matrix size - X
    image.imatrix_Y = fread(id, 1, 'short');  % Image matrix size - Y
    image.im_exno = fread(id, 1, 'ushort');  % Exam number for this image
    image.img_window = fread(id, 1, 'ushort');  % Window Value
    image.img_level = fread(id, 1, 'short');  % Level Value
    image.numecho = fread(id, 1, 'short');  % Number of echoes
    image.echonum = fread(id, 1, 'short');  % Echo Number
    image.im_uniq = fread(id, 1, 'short');  % The Make-Unique Flag
    image.im_seno = fread(id, 1, 'short');  % Series Number for this image
    image.im_no = fread(id, 1, 'short');  % Image Number
    image.contmode = fread(id, 1, 'short');  % Image Contrast Mode
    image.serrx = fread(id, 1, 'short');  % Series from which prescribed
    image.imgrx = fread(id, 1, 'short');  % Image from which prescribed
    image.screenformat = fread(id, 1, 'short');  % Screen Format(8/16 bit)
    image.plane = fread(id, 1, 'short');  % Plane Type
    image.im_compress = fread(id, 1, 'short');  % Image compression type for allocation
    image.im_scouttype = fread(id, 1, 'short');  % Scout Type (AP or lateral)
    image.contig = fread(id, 1, 'short');  % Continuous Slices Flag
    image.hrtrate = fread(id, 1, 'short');  % Cardiac Heart Rate (bpm)
    image.trgwindow = fread(id, 1, 'short');  % Trigger window (% of R-R interval)
    image.imgpcyc = fread(id, 1, 'short');  % Images per cardiac cycle
    
    image.obsolete3 = fread(id, 1, 'short');  % Actual Transmit Gain (.1 db)
    image.obsolete4 = fread(id, 1, 'short');  % Actual Receive Gain Analog (.1 db)
    image.obsolete5 = fread(id, 1, 'short');  % Actual Receive Gain Digital (.1 db)
    
    image.mr_flip = fread(id, 1, 'short');  % Flip Angle for GRASS scans (deg.)
    image.cphase = fread(id, 1, 'short');  % Total Cardiac Phase prescribed
    image.swappf = fread(id, 1, 'short');  % Swap Phase/Frequency Axis
    image.pauseint = fread(id, 1, 'short');  % Pause Interval (slices)
    
    image.obsolete6 = fread(id, 1, 'short');  % Auto Transmit Gain (0.1 dB)
    image.obsolete7 = fread(id, 1, 'short');  % PreScan R1 - Analog
    image.obsolete8 = fread(id, 1, 'short');  % PreScan R2 - Digital
    
    image.not_used_1 = fread(id, 1, 'short');  % Available for use
    image.imode = fread(id, 1, 'short');  % Imaging Mode
    image.pseq = fread(id, 1, 'short');  % Pulse Sequence
    image.pseqmode = fread(id, 1, 'short');  % Pulse Sequence Mode
    
    image.ctyp = fread(id, 1, 'short');  % Coil Type
    image.surfctyp = fread(id, 1, 'short');  % Surface Coil Type
    image.surfcext = fread(id, 1, 'short');  % Extremity Coil Flag
    image.supp_tech = fread(id, 1, 'short');  % SAT fat/water/none
    image.slquant = fread(id, 1, 'short');  % Number of slices in this scan group
    image.gpre = fread(id, 1, 'short');  % Graphically prescribed
    image.satbits = fread(id, 1, 'short');  % Bitmap of SAT selections
    image.scic = fread(id, 1, 'short');  % Surface Coil Intensity Correction Flag
    image.satxloc1 = fread(id, 1, 'short');  % R-side SAT pulse loc rel to lndmrk
    image.satxloc2 = fread(id, 1, 'short');  % L-side SAT pulse loc rel to lndmrk
    image.satyloc1 = fread(id, 1, 'short');  % A-side SAT pulse loc rel to lndmrk
    image.satyloc2 = fread(id, 1, 'short');  % P-side SAT pulse loc rel to lndmrk
    image.satzloc1 = fread(id, 1, 'short');  % S-side SAT pulse loc rel to lndmrk
    image.satzloc2 = fread(id, 1, 'short');  % I-side SAT pulse loc rel to lndmrk
    image.satxthick = fread(id, 1, 'short');  % Thickness of X-axis SAT pulse
    image.satythick = fread(id, 1, 'short');  % Thickness of Y-axis SAT pulse
    image.satzthick = fread(id, 1, 'short');  % Thickness of Z-axis SAT pulse
    image.flax = fread(id, 1, 'short');  % Phase contrast flow axis
    image.venc = fread(id, 1, 'short');  % Phase contrast velocity encoding
    image.thk_disclmr = fread(id, 1, 'short');  % Slice Thickness
    
    image.obsolete9 = fread(id, 1, 'short');  % Auto/Manual Prescan flag
    image.obsolete10 = fread(id, 1, 'short');  % Bitmap of changed values
    
    image.image_type = fread(id, 1, 'short');  % Magnitude, Phase, Imaginary, or Real
    image.vas_collapse = fread(id, 1, 'short');  % Collapse Image
    image.proj_alg = fread(id, 1, 'short');  % Projection Algorithm
    image.echo_trn_len = fread(id, 1, 'short');  % Echo Train Length for Fast Spin Echo
    image.frac_echo = fread(id, 1, 'short');  % Fractional Echo - Effective TE Flag
    image.prep_pulse = fread(id, 1, 'short');  % Preporatory Pulse Option
    image.cphasenum = fread(id, 1, 'short');  % Cardiac Phase Number
    image.var_echo = fread(id, 1, 'short');  % Variable Echo Flag
    image.scanactno = fread(id, 1, 'short');  % Scan Acquisition Number
    image.vasflags = fread(id, 1, 'short');  % Magnitude Weighting Flag
    image.integrity = fread(id, 1, 'short');  % GE Image Integrity
    image.freq_dir = fread(id, 1, 'short');  % Frequency Direction
    image.vas_mode = fread(id, 1, 'short');  % Vascular Mode
    
    image.pscopts = fread(id, 1, 'short');  % bitmap of prescan options
    
    image.obsolete11 = fread(id, 1, 'short');  % gradient offset in X-direction
    image.obsolete12 = fread(id, 1, 'short');  % gradient offset in Y-direction
    image.obsolete13 = fread(id, 1, 'short');  % gradient offset in Z-direction
    
    image.unoriginal = fread(id, 1, 'short');  % identifies image as original or unoriginal
    image.interleaves = fread(id, 1, 'short');  % number of EPI shots
    image.effechospace = fread(id, 1, 'short');  % effective echo spacing for EPI
    image.viewsperseg = fread(id, 1, 'short');  % views per segment
    image.rbpm = fread(id, 1, 'short');  % respiratory rate, breaths per min
    image.rtpoint = fread(id, 1, 'short');  % respiratory trigger point as percent of max.
    image.rcvrtype = fread(id, 1, 'short');  % type of receiver used
    image.sarMode = fread(id, 1, 'short');  % Sar Ctrl Mode (Normal, 1st or 2nd)
    image.dBdtMode = fread(id, 1, 'short');  % dBdt Ctrl Mode (Normal, 1st or 2nd)
    image.govBody = fread(id, 1, 'short');  % Governing Body MHW/IEC/FDA
    image.sarDefinition = fread(id, 1, 'short');  % Defined terms avaialble
    image.no_shimvol = fread(id, 1, 'short');  %
    image.shim_vol_type = fread(id, 1, 'short');  %
    
    image.current_phase = fread(id, 1, 'short');  % Current phase
    image.art_level = fread(id, 1, 'short');  % Art level
    image.slice_group_numer = fread(id, 1, 'short');  % Slice group number
    image.number_of_slice_groups = fread(id, 1, 'short');  % # of slc groups
    image.show_in_autoview = fread(id, 1, 'short');  %
    image.slice_number_inGroup = fread(id, 1, 'short');  %
    
    image.short_padding = fread(id, 78, 'char'); % short padding
    
    image.psdname = modchar(fread(id, 33, '*char'));  % Pulse Sequence Name
    
    fseek(id,pos_image+1665,-1);  % there is a bug here if no fseek, why?
    
    image.proj_name = modchar(fread(id, 13, '*char'));  % Projection Algorithm Name
    image.psd_iname = modchar(fread(id, 13, '*char'));  % PSD name from inside PSD
    image.im_diskid = modchar(fread(id, 1, '*char'));  % Disk ID for this Image
    image.pdid = modchar(fread(id, 14, '*char'));  % Pixel Data ID
    image.im_suid = modchar(fread(id, 4, '*char'));  % Suite id for this image
    image.contrastIV = modchar(fread(id, 17, '*char'));  % IV Contrast Agent
    image.contrastOral = modchar(fread(id, 17, '*char'));  % Oral Contrast Agent
    image.loc_ras = modchar(fread(id, 1, '*char'));  % RAS letter of image location
    image.forimgrev = modchar(fread(id, 4, '*char'));  % Foreign Image Revision
    image.cname = modchar(fread(id, 17, '*char'));  % Coil Name
    image.im_verscre = modchar(fread(id, 2, '*char'));  % Genesis Version - Created
    image.im_verscur = modchar(fread(id, 2, '*char'));  % Genesis Version - Now
    image.im_alloc_key = modchar(fread(id, 13, '*char'));  %
    image.ref_img = modchar(fread(id, 1, '*char'));  % Reference Image Field
    image.sum_img = modchar(fread(id, 1, '*char'));  % Summary Image Field
    image.filter_mode = modchar(fread(id, 16, '*char'));  % String Slop Field 1
    image.slop_str_2 = modchar(fread(id, 16, '*char'));  % String Slop Field 2
    image.image_uid = moduid(fread(id, 32, 'uint8'));  % Image Unique ID
    image.sop_uid = moduid(fread(id, 32, 'uint8'));  % Service Obj Class Unique ID
    image.GEcname = modchar(fread(id, 24, '*char'));  % GECoilname for the cname
    image.usedCoilData = modchar(fread(id, 100, '*char'));  % Concatenated str of coilcode and chip serialID
    image.astcalseriesuid = moduid(fread(id, 32, 'uint8'));  %
    image.purecalseriesuid = moduid(fread(id, 32, 'uint8'));  %
    
    image.xml_psc_shm_vol = moduid(fread(id, 32, 'char'));  %
    image.rxmpath = moduid(fread(id, 64, 'char'));  %
    image.psdnameannot = moduid(fread(id, 33, 'char'));  %
    
    image.img_hdr_padding = modchar(fread(id, 250, '*char'));  %
end
end

function prescan = read_prescan(id)
prescan.command = fread(id, 1, 'long');
prescan.mps_r1 = fread(id, 1, 'int');
prescan.mps_r2 = fread(id, 1, 'int');
prescan.mps_tg = fread(id, 1, 'int');
prescan.mps_freq = fread(id, 1, 'uint32');
prescan.aps_r1 = fread(id, 1, 'int');
prescan.aps_r2 = fread(id, 1, 'int');
prescan.aps_tg = fread(id, 1, 'int');
prescan.aps_freq = fread(id, 1, 'uint32');
prescan.scalei = fread(id, 1, 'float');
prescan.scaleq = fread(id, 1, 'float');
prescan.snr_warning= fread(id, 1, 'int');
prescan.aps_or_mps= fread(id, 1, 'int');
prescan.maps_bitmap= fread(id, 1, 'int');
prescan.powerspec = modchar(fread(id, 256, '*char'));
prescan.filter1 = fread(id, 1, 'int');
prescan.filter2 = fread(id, 1, 'int');
prescan.xshim = fread(id, 1, 'short');
prescan.yshim = fread(id, 1, 'short');
prescan.zshim = fread(id, 1, 'short');
prescan.recon_enable = fread(id, 1, 'short');
prescan.autoshim_status = fread(id, 1, 'int');

prescan.rec_std0 = fread(id, 1, 'float');
prescan.rec_std1 = fread(id, 1, 'float');
prescan.rec_std2 = fread(id, 1, 'float');
prescan.rec_std3 = fread(id, 1, 'float');
prescan.rec_std4 = fread(id, 1, 'float');
prescan.rec_std5 = fread(id, 1, 'float');
prescan.rec_std6 = fread(id, 1, 'float');
prescan.rec_std7 = fread(id, 1, 'float');
prescan.rec_std8 = fread(id, 1, 'float');
prescan.rec_std9 = fread(id, 1, 'float');
prescan.rec_std10 = fread(id, 1, 'float');
prescan.rec_std11 = fread(id, 1, 'float');
prescan.rec_std12 = fread(id, 1, 'float');
prescan.rec_std13 = fread(id, 1, 'float');
prescan.rec_std14 = fread(id, 1, 'float');
prescan.rec_std15 = fread(id, 1, 'float');
prescan.rec_std16 = fread(id, 1, 'float');
prescan.rec_std17 = fread(id, 1, 'float');
prescan.rec_std18 = fread(id, 1, 'float');
prescan.rec_std19 = fread(id, 1, 'float');
prescan.rec_std20 = fread(id, 1, 'float');
prescan.rec_std21 = fread(id, 1, 'float');
prescan.rec_std22 = fread(id, 1, 'float');
prescan.rec_std23 = fread(id, 1, 'float');
prescan.rec_std24 = fread(id, 1, 'float');
prescan.rec_std25 = fread(id, 1, 'float');
prescan.rec_std26 = fread(id, 1, 'float');
prescan.rec_std27 = fread(id, 1, 'float');
prescan.rec_std28 = fread(id, 1, 'float');
prescan.rec_std29 = fread(id, 1, 'float');
prescan.rec_std30 = fread(id, 1, 'float');
prescan.rec_std31 = fread(id, 1, 'float');
prescan.rec_std32 = fread(id, 1, 'float');
prescan.rec_std33 = fread(id, 1, 'float');
prescan.rec_std34 = fread(id, 1, 'float');
prescan.rec_std35 = fread(id, 1, 'float');
prescan.rec_std36 = fread(id, 1, 'float');
prescan.rec_std37 = fread(id, 1, 'float');
prescan.rec_std38 = fread(id, 1, 'float');
prescan.rec_std39 = fread(id, 1, 'float');
prescan.rec_std40 = fread(id, 1, 'float');
prescan.rec_std41 = fread(id, 1, 'float');
prescan.rec_std42 = fread(id, 1, 'float');
prescan.rec_std43 = fread(id, 1, 'float');
prescan.rec_std44 = fread(id, 1, 'float');
prescan.rec_std45 = fread(id, 1, 'float');
prescan.rec_std46 = fread(id, 1, 'float');
prescan.rec_std47 = fread(id, 1, 'float');
prescan.rec_std48 = fread(id, 1, 'float');
prescan.rec_std49 = fread(id, 1, 'float');
prescan.rec_std50 = fread(id, 1, 'float');
prescan.rec_std51 = fread(id, 1, 'float');
prescan.rec_std52 = fread(id, 1, 'float');
prescan.rec_std53 = fread(id, 1, 'float');
prescan.rec_std54 = fread(id, 1, 'float');
prescan.rec_std55 = fread(id, 1, 'float');
prescan.rec_std56 = fread(id, 1, 'float');
prescan.rec_std57 = fread(id, 1, 'float');
prescan.rec_std58 = fread(id, 1, 'float');
prescan.rec_std59 = fread(id, 1, 'float');
prescan.rec_std60 = fread(id, 1, 'float');
prescan.rec_std61 = fread(id, 1, 'float');
prescan.rec_std62 = fread(id, 1, 'float');
prescan.rec_std63 = fread(id, 1, 'float');
prescan.rec_std64 = fread(id, 1, 'float');
prescan.rec_std65 = fread(id, 1, 'float');
prescan.rec_std66 = fread(id, 1, 'float');
prescan.rec_std67 = fread(id, 1, 'float');
prescan.rec_std68 = fread(id, 1, 'float');
prescan.rec_std69 = fread(id, 1, 'float');
prescan.rec_std70 = fread(id, 1, 'float');
prescan.rec_std71 = fread(id, 1, 'float');
prescan.rec_std72 = fread(id, 1, 'float');
prescan.rec_std73 = fread(id, 1, 'float');
prescan.rec_std74 = fread(id, 1, 'float');
prescan.rec_std75 = fread(id, 1, 'float');
prescan.rec_std76 = fread(id, 1, 'float');
prescan.rec_std77 = fread(id, 1, 'float');
prescan.rec_std78 = fread(id, 1, 'float');
prescan.rec_std79 = fread(id, 1, 'float');
prescan.rec_std80 = fread(id, 1, 'float');
prescan.rec_std81 = fread(id, 1, 'float');
prescan.rec_std82 = fread(id, 1, 'float');
prescan.rec_std83 = fread(id, 1, 'float');
prescan.rec_std84 = fread(id, 1, 'float');
prescan.rec_std85 = fread(id, 1, 'float');
prescan.rec_std86 = fread(id, 1, 'float');
prescan.rec_std87 = fread(id, 1, 'float');
prescan.rec_std88 = fread(id, 1, 'float');
prescan.rec_std89 = fread(id, 1, 'float');
prescan.rec_std90 = fread(id, 1, 'float');
prescan.rec_std91 = fread(id, 1, 'float');
prescan.rec_std92 = fread(id, 1, 'float');
prescan.rec_std93 = fread(id, 1, 'float');
prescan.rec_std94 = fread(id, 1, 'float');
prescan.rec_std95 = fread(id, 1, 'float');
prescan.rec_std96 = fread(id, 1, 'float');
prescan.rec_std97 = fread(id, 1, 'float');
prescan.rec_std98 = fread(id, 1, 'float');
prescan.rec_std99 = fread(id, 1, 'float');
prescan.rec_std100 = fread(id, 1, 'float');
prescan.rec_std101 = fread(id, 1, 'float');
prescan.rec_std102 = fread(id, 1, 'float');
prescan.rec_std103 = fread(id, 1, 'float');
prescan.rec_std104 = fread(id, 1, 'float');
prescan.rec_std105 = fread(id, 1, 'float');
prescan.rec_std106 = fread(id, 1, 'float');
prescan.rec_std107 = fread(id, 1, 'float');
prescan.rec_std108 = fread(id, 1, 'float');
prescan.rec_std109 = fread(id, 1, 'float');
prescan.rec_std110 = fread(id, 1, 'float');
prescan.rec_std111 = fread(id, 1, 'float');
prescan.rec_std112 = fread(id, 1, 'float');
prescan.rec_std113 = fread(id, 1, 'float');
prescan.rec_std114 = fread(id, 1, 'float');
prescan.rec_std115 = fread(id, 1, 'float');
prescan.rec_std116 = fread(id, 1, 'float');
prescan.rec_std117 = fread(id, 1, 'float');
prescan.rec_std118 = fread(id, 1, 'float');
prescan.rec_std119 = fread(id, 1, 'float');
prescan.rec_std120 = fread(id, 1, 'float');
prescan.rec_std121 = fread(id, 1, 'float');
prescan.rec_std122 = fread(id, 1, 'float');
prescan.rec_std123 = fread(id, 1, 'float');
prescan.rec_std124 = fread(id, 1, 'float');
prescan.rec_std125 = fread(id, 1, 'float');
prescan.rec_std126 = fread(id, 1, 'float');
prescan.rec_std127 = fread(id, 1, 'float');

prescan.rec_mean0 = fread(id, 1, 'float');
prescan.rec_mean1 = fread(id, 1, 'float');
prescan.rec_mean2 = fread(id, 1, 'float');
prescan.rec_mean3 = fread(id, 1, 'float');
prescan.rec_mean4 = fread(id, 1, 'float');
prescan.rec_mean5 = fread(id, 1, 'float');
prescan.rec_mean6 = fread(id, 1, 'float');
prescan.rec_mean7 = fread(id, 1, 'float');
prescan.rec_mean8 = fread(id, 1, 'float');
prescan.rec_mean9 = fread(id, 1, 'float');
prescan.rec_mean10 = fread(id, 1, 'float');
prescan.rec_mean11 = fread(id, 1, 'float');
prescan.rec_mean12 = fread(id, 1, 'float');
prescan.rec_mean13 = fread(id, 1, 'float');
prescan.rec_mean14 = fread(id, 1, 'float');
prescan.rec_mean15 = fread(id, 1, 'float');
prescan.rec_mean16 = fread(id, 1, 'float');
prescan.rec_mean17 = fread(id, 1, 'float');
prescan.rec_mean18 = fread(id, 1, 'float');
prescan.rec_mean19 = fread(id, 1, 'float');
prescan.rec_mean20 = fread(id, 1, 'float');
prescan.rec_mean21 = fread(id, 1, 'float');
prescan.rec_mean22 = fread(id, 1, 'float');
prescan.rec_mean23 = fread(id, 1, 'float');
prescan.rec_mean24 = fread(id, 1, 'float');
prescan.rec_mean25 = fread(id, 1, 'float');
prescan.rec_mean26 = fread(id, 1, 'float');
prescan.rec_mean27 = fread(id, 1, 'float');
prescan.rec_mean28 = fread(id, 1, 'float');
prescan.rec_mean29 = fread(id, 1, 'float');
prescan.rec_mean30 = fread(id, 1, 'float');
prescan.rec_mean31 = fread(id, 1, 'float');
prescan.rec_mean32 = fread(id, 1, 'float');
prescan.rec_mean33 = fread(id, 1, 'float');
prescan.rec_mean34 = fread(id, 1, 'float');
prescan.rec_mean35 = fread(id, 1, 'float');
prescan.rec_mean36 = fread(id, 1, 'float');
prescan.rec_mean37 = fread(id, 1, 'float');
prescan.rec_mean38 = fread(id, 1, 'float');
prescan.rec_mean39 = fread(id, 1, 'float');
prescan.rec_mean40 = fread(id, 1, 'float');
prescan.rec_mean41 = fread(id, 1, 'float');
prescan.rec_mean42 = fread(id, 1, 'float');
prescan.rec_mean43 = fread(id, 1, 'float');
prescan.rec_mean44 = fread(id, 1, 'float');
prescan.rec_mean45 = fread(id, 1, 'float');
prescan.rec_mean46 = fread(id, 1, 'float');
prescan.rec_mean47 = fread(id, 1, 'float');
prescan.rec_mean48 = fread(id, 1, 'float');
prescan.rec_mean49 = fread(id, 1, 'float');
prescan.rec_mean50 = fread(id, 1, 'float');
prescan.rec_mean51 = fread(id, 1, 'float');
prescan.rec_mean52 = fread(id, 1, 'float');
prescan.rec_mean53 = fread(id, 1, 'float');
prescan.rec_mean54 = fread(id, 1, 'float');
prescan.rec_mean55 = fread(id, 1, 'float');
prescan.rec_mean56 = fread(id, 1, 'float');
prescan.rec_mean57 = fread(id, 1, 'float');
prescan.rec_mean58 = fread(id, 1, 'float');
prescan.rec_mean59 = fread(id, 1, 'float');
prescan.rec_mean60 = fread(id, 1, 'float');
prescan.rec_mean61 = fread(id, 1, 'float');
prescan.rec_mean62 = fread(id, 1, 'float');
prescan.rec_mean63 = fread(id, 1, 'float');
prescan.rec_mean64 = fread(id, 1, 'float');
prescan.rec_mean65 = fread(id, 1, 'float');
prescan.rec_mean66 = fread(id, 1, 'float');
prescan.rec_mean67 = fread(id, 1, 'float');
prescan.rec_mean68 = fread(id, 1, 'float');
prescan.rec_mean69 = fread(id, 1, 'float');
prescan.rec_mean70 = fread(id, 1, 'float');
prescan.rec_mean71 = fread(id, 1, 'float');
prescan.rec_mean72 = fread(id, 1, 'float');
prescan.rec_mean73 = fread(id, 1, 'float');
prescan.rec_mean74 = fread(id, 1, 'float');
prescan.rec_mean75 = fread(id, 1, 'float');
prescan.rec_mean76 = fread(id, 1, 'float');
prescan.rec_mean77 = fread(id, 1, 'float');
prescan.rec_mean78 = fread(id, 1, 'float');
prescan.rec_mean79 = fread(id, 1, 'float');
prescan.rec_mean80 = fread(id, 1, 'float');
prescan.rec_mean81 = fread(id, 1, 'float');
prescan.rec_mean82 = fread(id, 1, 'float');
prescan.rec_mean83 = fread(id, 1, 'float');
prescan.rec_mean84 = fread(id, 1, 'float');
prescan.rec_mean85 = fread(id, 1, 'float');
prescan.rec_mean86 = fread(id, 1, 'float');
prescan.rec_mean87 = fread(id, 1, 'float');
prescan.rec_mean88 = fread(id, 1, 'float');
prescan.rec_mean89 = fread(id, 1, 'float');
prescan.rec_mean90 = fread(id, 1, 'float');
prescan.rec_mean91 = fread(id, 1, 'float');
prescan.rec_mean92 = fread(id, 1, 'float');
prescan.rec_mean93 = fread(id, 1, 'float');
prescan.rec_mean94 = fread(id, 1, 'float');
prescan.rec_mean95 = fread(id, 1, 'float');
prescan.rec_mean96 = fread(id, 1, 'float');
prescan.rec_mean97 = fread(id, 1, 'float');
prescan.rec_mean98 = fread(id, 1, 'float');
prescan.rec_mean99 = fread(id, 1, 'float');
prescan.rec_mean100 = fread(id, 1, 'float');
prescan.rec_mean101 = fread(id, 1, 'float');
prescan.rec_mean102 = fread(id, 1, 'float');
prescan.rec_mean103 = fread(id, 1, 'float');
prescan.rec_mean104 = fread(id, 1, 'float');
prescan.rec_mean105 = fread(id, 1, 'float');
prescan.rec_mean106 = fread(id, 1, 'float');
prescan.rec_mean107 = fread(id, 1, 'float');
prescan.rec_mean108 = fread(id, 1, 'float');
prescan.rec_mean109 = fread(id, 1, 'float');
prescan.rec_mean110 = fread(id, 1, 'float');
prescan.rec_mean111 = fread(id, 1, 'float');
prescan.rec_mean112 = fread(id, 1, 'float');
prescan.rec_mean113 = fread(id, 1, 'float');
prescan.rec_mean114 = fread(id, 1, 'float');
prescan.rec_mean115 = fread(id, 1, 'float');
prescan.rec_mean116 = fread(id, 1, 'float');
prescan.rec_mean117 = fread(id, 1, 'float');
prescan.rec_mean118 = fread(id, 1, 'float');
prescan.rec_mean119 = fread(id, 1, 'float');
prescan.rec_mean120 = fread(id, 1, 'float');
prescan.rec_mean121 = fread(id, 1, 'float');
prescan.rec_mean122 = fread(id, 1, 'float');
prescan.rec_mean123 = fread(id, 1, 'float');
prescan.rec_mean124 = fread(id, 1, 'float');
prescan.rec_mean125 = fread(id, 1, 'float');
prescan.rec_mean126 = fread(id, 1, 'float');
prescan.rec_mean127 = fread(id, 1, 'float');

prescan.line_width = fread(id, 1, 'int');
prescan.ws_flip = fread(id, 1, 'int');
prescan.supp_lvl = fread(id, 1, 'int');
prescan.psc_reuse = fread(id, 1, 'int');
prescan.psc_reuse_string= modchar(fread(id, 52, 'char'));
prescan.psc_ta = fread(id, 1, 'int');
prescan.phase_correction_status = fread(id, 1, 'int');
prescan.broad_band_select = fread(id, 1, 'int');
prescan.buffer = modchar(fread(id, 64, 'char'));
end

function daqtab = read_daqtab(id, ver, n)

if (ver >= 14)    
    for k=1:n
        daqtab(k).pass_number = fread(id, 1, 'short');
        daqtab(k).slice_in_pass = fread(id, 1, 'short');
        daqtab(k).gw_point = fread(id, [3,3], 'float');
        daqtab(k).transpose = fread(id, 1, 'short');
        daqtab(k).rotate = fread(id, 1, 'short');
        daqtab(k).swiftcoilinfo = fread(id, 1, 'short');
        fseek(id,2,'cof');
    end
elseif (ver == 11)
    for k=1:n
        daqtab(k).pass_number = fread(id, 1, 'short');
        daqtab(k).slice_in_pass = fread(id, 1, 'short');
        daqtab(k).gw_point = fread(id, [3,3], 'float');
        daqtab(k).transpose = fread(id, 1, 'short');
        daqtab(k).rotate = fread(id, 1, 'short');
    end
else    
    for k=1:n
        daqtab(k).pass_number = fread(id, 1, 'short');
        daqtab(k).slice_in_pass = fread(id, 1, 'short');
        daqtab(k).gw_point = fread(id, [3,3], 'float');      
    end   
end

end


function rdb = read_multi_rcv(id, rdb)
buff = fread(id, 8, 'short');
rdb.dab_start_rcv = buff(1:2:end);
rdb.dab_stop_rcv = buff(2:2:end);
end


function v = read_vartype(id);
v = fread(id, 1, 'ulong');
fseek(id, 4, 0);
end


function a = moduid(x)
a = [];
for k=1:length(x)  % loop over each 8 bit value...
    n = bitshift(bitand(x(k),240),-4) - 1;
    if (n < 0), break; end
    if (n > 9), a=[a,'.']; else a=[a,num2str(n)]; end
    n = bitand(x(k),15) - 1;
    if (n < 0), break; end
    if (n > 9), a=[a,'.']; else a=[a,num2str(n)]; end
end
a = sprintf('%s', a);
end


function a = modchar(x)
warning off;
a = [];
for k=1:length(x)  % loop over each character...
    if (int32(x(k)) < 32 || int32(x(k)) > 126), break; end  % not a printing char
    a = [a, x(k)];
end
a = sprintf('%s', a);
end

