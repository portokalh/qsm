function img = RussReconADMask_ME(varargin)
%%function raw = Pseries(varargin)
% Before calling the function, put the pfiles (single- or multi-echo) in a
% separate folder. Call the function in the created folder. 
%
%   
% varargin is a list of options in 'single quotes' separated by commas
%   (no input)  - output is complex, raw, frequency space data
%   'raw'       - output is complex data in k-space
%   'echofill'  - replaces faulty data in last echo when phase1*phase2 is
%                 not divisible by baseline_spacing (usually 2048)
%   'shift'     - linearly corrects even echo shift to match odd echoes
%   [-1 -1 -1]  - flips readout, phase1, and phase2 array orientations
%                 (good for getting data into WHS if you scanned it weird)
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
rawswitch = false;
cryomask = false;
echoflag = false;
shift9flag = false;
maskinflag = false;
fermi_flag = false;
% if matlabpool('size') > 0
%     pswitch = true;
% end
for narg = 1:nargin
    if ischar(varargin{narg})
        if strcmp(varargin{narg},'theworks')
            imgswitch = true;
            shiftswitch = true;
            phaseswitch = true;
            maskswitch = true;
        elseif strcmp(varargin{narg},'bruker')
            bruker_flag = true;
        elseif strcmp(varargin{narg},'cryomask')
            cryomask_flag = true;
        elseif strcmp(varargin{narg},'multimask')
            multimask_flag = true;
            maskswitch = true;
        elseif strcmp(varargin{narg},'Ashift9')
            shift9flag = true;
        elseif strcmp(varargin{narg},'agilent')
            agilent_postrecon = true;
        elseif strcmp(varargin{narg},'nomask')
            imgswitch = true;
            shiftswitch = true;
            phaseswitch = true;
        elseif strcmp(varargin{narg},'raw')
            rawswitch = true;
        elseif strcmp(varargin{narg},'X')
            Xswitch = true;
        elseif strcmp(varargin{narg},'post')
            postprocessingonly = true;
        elseif strcmp(varargin{narg},'zcorrect')
            zcorrect_flag = true;
        elseif strcmp(varargin{narg},'fermi')
            fermi_flag = true;
        elseif strcmp(varargin{narg},'parallel')
            pswitch = true;
            if matlabpool('size') == 0
                matlabpool open
            end
        elseif strcmp(varargin{narg},'rp')
            rp_flag = true;
            rp_file = dir('*.rp');
        elseif strcmp(varargin{narg},'T2')
            T2switch = true;
        elseif strcmp(varargin{narg},'shift')
            shiftswitch = true;
        elseif strcmp(varargin{narg},'phaseshift')
            phaseswitch = true;
        elseif strcmp(varargin{narg},'mask')
            maskswitch = true;
        elseif strcmp(varargin{narg},'center')
            centerswitch = true;
            imgswitch = true;
        elseif strcmp(varargin{narg},'int')
            pfile_type = 'int';
        elseif strcmp(varargin{narg},'save')
            matsave = true;
        elseif strcmp(varargin{narg},'short')
            pfile_type = 'short';    
        elseif strcmp(varargin{narg},'convert')
            SHORT_SIZE = 2; % in bytes
            convertflag = true;
        elseif strcmp(varargin{narg}(1:4),'mask')
            maskinflag = true;
            mask = open_nii(varargin{narg});
        elseif length(varargin{narg}) < 5
            scan_name = varargin{narg};
        elseif strcmp(varargin{narg}(1:4),'echo')
            echoflag = true;
            iecho = str2num(varargin{narg}(5:end));
        elseif strcmp(varargin{narg}(1:4),'BIAC')
            BIACswitch = true;
            maskswitch = true;
%             iBIAC = str2num(varargin{narg}(5:end));
            runno = str2num(varargin{narg}(5:end));
        elseif ~strcmp(varargin{narg}(1:4),'BIAC')
            scan_name = varargin{narg};
        end
    else
        if length(varargin{narg}) == 3    
            if length(varargin{narg}(abs(varargin{narg}) > 1)) > 1
                permuteswitch = true;
                permute_vector = varargin{narg};
                save permute_vector permute_vector
                permute_vector = [permute_vector 4];
            else
                orientswitch = true;
                orientation = varargin{narg};
            end
        end
    end   
end

%% Recon summary



%% Convert GE p-file series by removing baselines and creating echo files

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
%     baseline_spacing = hdr.rdb.nframes;
    SHORT_SIZE = 2;  
%     flag='l';  % default byte swap
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
    echo_spacing = hdr.rdb.user33/1e3; % Only valid in GE multi-gradient echo
    TE = (0:dims(4)-1)*echo_spacing+hdr.rdb.te/1e6; save TE TE
    save hdr hdr header_bytes

    raw=zeros(dims(1),dims(2),dims(3));
    back = 0; time = 0;
    fid = fopen(rp_file(1).name,'r','l');
    [~]=fread(fid,[1 header_bytes/SHORT_SIZE],'short');  % skip header
    for z=1:dims(3)
        [back,time] = progress(z,dims(3),'Reading slice',back,time);
            for y=1:dims(2)    
                dump=fread(fid,[2 dims(1)], pfile_type);
                raw(:,y,z)=squeeze(dump(2,1:end)+1i*dump(1,:));
            end
    end
    fclose(fid);   
    
end


%% Read in raw data, (Agilent, Bruker or GE)

if exist('bruker_flag','var')   %%% Bruker Style
    
    B0 = 7;
    
    if ~exist('postprocessingonly','var')
        [raw,hdr] = BrukerReconRawMGEMS;
    else
        hdr.acqp = readBrukerHeader('acqp');
        hdr.method = readBrukerHeader('method');
%         raw = readBrukerFID('',hdr.method);
        hdr.subject = readBrukerHeader('../subject');
        load img2channels
        raw = permute(img,[1 2 3 5 4]); clear img;
     
    end
    
    if strcmp(hdr.method.PVM_SPackArrSliceOrient,'coronal')
        B0dir = [1 0 0]; 
        BIACpermute = [2 3 1];
    elseif strcmp(hdr.method.PVM_SPackArrSliceOrient,'axial')
        B0dir = [0 0 1];
        BIACpermute = [1 2 3];
    elseif strcmp(hdr.method.PVM_SPackArrSliceOrient,'sagittal')
        B0dir = [1 0 0];
        BIACpermute = [3 2 1];
    end
    scan_orientation = hdr.method.PVM_SPackArrSliceOrient;
    save B0dir B0dir BIACpermute scan_orientation
    
    
     
    FOV = hdr.method.PVM_Fov;
    raw = single(raw);

%     vox = hdr.method.PVM_SpatResol; 
     
    %%% Pad the array here if desired, e.g.: raw = padarray(raw,[0 0 88 0],'both');
%     dims = [hdr.method.PVM_Matrix hdr.method.PVM_NEchoImages];
    dims = size(raw); 
    dims = dims(1:4);
    vox = FOV./dims(1:3);
    if length(vox) < 3
        if strcmp(hdr.method.Method,'MDEFT');
            FOV(3) = hdr.method.PVM_SliceThick*hdr.method.PVM_SPackArrNSlices;
            TE = hdr.method.PVM_EchoTime1/1000; 
        else
            FOV(3) = hdr.method.PVM_SliceThick;
            TE = hdr.method.EffectiveTE/1000;
        end
        vox(3) = hdr.method.PVM_SliceThick;
    else
        TE = hdr.method.EffectiveTE/1000;
    end
    save vox vox
    % Corrects Bruker's shift in phase encode 2 due to bad pulses.
    if dims(4) > 1 
        if exist('zcorrect_flag','var')
            raw = shiftMEbruker(raw); 
        end
    end  
    
    if ~exist('postprocessingonly','var')
        if strcmp(hdr.method.PVM_SpatDimEnum,'2D') 
            fw = fermiW([dims(1:2) 1]);
            for s = 1:dims(3)
                for k = 1:dims(4)  
                    for icoil = 1:size(raw,5) 
                        if fermi_flag
                            raw(:,:,s,k,icoil) = ifftnc(fw.*raw(:,:,s,k,icoil));
                        else
                            raw(:,:,s,k,icoil) = ifftnc(raw(:,:,s,k,icoil));
                        end
                    end
                end
            end 
        else
            fw = fermiW(dims(1:3));
            back = 0; time = 0;
            for k = 1:dims(4) 
                [back,time] = progress(k,dims(4),'Performing IFFT on echo',back,time);
                for icoil = 1:size(raw,5)
                    if fermi_flag
                        raw(:,:,:,k,icoil) = ifftnc(fw.*raw(:,:,:,k,icoil));
                    else
                        raw(:,:,:,k,icoil) = ifftnc(raw(:,:,:,k,icoil));
                    end
                end
            end
        end
    end
    
    % Centers image object in FOV
    if exist('centerswitch','var')
        if strcmp(hdr.method.PVM_SpatDimEnum,'3D')
            fprintf('   Centering image data...\n');
%             raw = centerME(raw); % manual centering is lame/slow!
            if exist(fullfile(pwd,'shift_vector.mat'),'file')
%                 load shift_vector
%                 [raw,shift_vector] = autoC(raw,shift_vector);
                [raw,shift_vector] = autoC(raw);

            else
                [raw,shift_vector] = autoC(raw);
            end
            save shift_vector shift_vector
        end
        
    else
        shift_vector = [0 0 0 0 0];
    end
    
%     save raw_coil_images raw
    
    % Combine data from multiple coils (like the cryoprobe setup)
%     back = 0; time = 0;
    img = zeros(dims);
    
    if strcmp(hdr.method.Method,'RARE')
        acqmethod = 'SE';
    else
        acqmethod = 'GRE'; 
    end
    
    %%% masking
%calculate noise variance to threshold phase info
    phivar = var7(unwrap_phase(angle(raw(:,:,:,1,1)),[1 1 1],[0 0 0]));% raw is complex x,y,z, echo, channel; 111 is voxel size, assumed isotropic; 0 0 0 is padding
    mag1 = mean(abs(raw(:,:,:,1,:)),5);
    [~,NoiseParams,SignalParams] = autoSNRpipe(mag1);
    noise_floor = NoiseParams.Mean;
    noise_thresh = noise_floor*6;%empirical threshold, below this the phase ino is not reliable
    threshmask = (phivar < .2)|(mag1 > noise_thresh);
    phimask = threshmask;
    n_ops = 4;
    for k = 1:n_ops
        phimask = imerode(phimask,strel3d_odd(3));
    end
    phimask = scrub(phimask,6,1:2);% keep the two larges compmenents is this case when you have a tube, otherwise just keep 1 not 1:2
    for k = 1:n_ops-1
        phimask = imdilate(phimask,strel3d_odd(3));
    end  

    maskME = false(dims);
    for k = 1:dims(4) % use 3*noise_floor for non BOMUS, 2 for BOMUS
         maskME(:,:,:,k) = imerode(single(phimask).*mean(abs(raw(:,:,:,k,:)),5) > 3*noise_floor,strel3d_odd(3));
        maskME(:,:,:,k) = single(phimask).*mean(abs(raw(:,:,:,k,:)),5) > 3*noise_floor; 
       maskME(:,:,:,k) = single(phimask).*mean(abs(raw(:,:,:,k,:)),5) > 3*noise_floor; %alex reduced to 2* noise_floor

    end
     save_nii(make_nii(uint8(maskME),vox,[0 0 0],2),['maskME.nii']);
    maskinflag = true;
   
    if strcmp(acqmethod,'GRE');
        if hdr.method.PVM_EncNReceivers > 1
            for k = 1:dims(4)
                    fprintf('   Combining coil data for echo %d/%d...\n',k,dims(4));

%                 [back,time] = progress(k,dims(4),'Combining coil data for echo',back,time);
                if dims(3) == 1
                    img(:,:,1,k) = coil_combine_complex(permute(raw(:,:,1,k,:),[1 2 3 5 4])); 
                else
                    if maskinflag
                        if k <= size(maskME,4)
                            mask2 = maskME(:,:,:,k);
                        else
                            mask2 = maskME(:,:,:,1);
                        end
                        img(:,:,:,k) = coil_combine_complex(squeeze(raw(:,:,:,k,:)),mask2);% this does the phase processing, unwrapping, bacgroud removal
                    else
                        img(:,:,:,k) = coil_combine_complex(squeeze(raw(:,:,:,k,:)));
                    end
                end
            end
        else
            img = raw;% this needs to do the phase processing, unwrapping, bacgroud removal - add from coil_combine
        end
    else
        img = squeeze(mean(abs(raw),5));% this needs to do the phase processing, unwrapping, bacgroud removal - add from coil_combine
    end
    autoSNR(squeeze(mean(abs(raw),5)));
    load meanSNR meanSNR
     [m1,m2,m3,m4] = size(raw);
     disp(m1)
     disp(m2)
     disp(m3)
     disp(m4)
     T2 = T2mapMaskedWeighted(mean(abs(raw),5),TE,phimask(:,:,:,1),meanSNR);
    mycoeffs = Xcoeff(B0,TE);
    Freq = zeros(size(img));
    for k = 1:dims(4)
        Freq(:,:,:,k) = angle(img(:,:,:,k))*mycoeffs.Freq(k);
    end
     idx = 1:6;%alex
    %idx = 1;
   [W,FreqC] = nanMEW(TE(idx),T2,Freq(:,:,:,idx));%nanMEW - do not use empty voxels from Freq; average meanaingul voxels only
    save_nii(make_nii(Freq,vox,[0 0 0],16),'Freq.nii');
    save_nii(make_nii(FreqC,vox,[0 0 0],16),'FreqC.nii');
    save_nii(make_nii(uint8(maskME),vox,[0 0 0],2),'maskME.nii');

    save TE TE
    save FOV FOV 
    save dims dims
    save hdr hdr
    
elseif agilent_postrecon
    load dims
    load recon_data
    load TE
    FOV = [procpar.lro procpar.lpe procpar.lpe2]*10;
    save FOV FOV
    vox = FOV./dims(1:3);
    save vox vox

    if ~echoflag
       iecho = 1:dims(4);
    end
    nechoes = length(iecho);
    
    
    raw = readMEraw('k','','single',iecho);
    raw = phase_ramp_remove1(raw);
    
%     if shift9flag == 1
%         display('   Shifting echoes from delayed Agilent acquisition...');
%         raw = circshift(raw,[-9 0 0 0]);
%     end
    
    if fermi_flag
        raw = raw.*fermi_filterV(dims(1:3),.15,.75);
    end
    
    

    
    if ~rawswitch
        % Perform IFFT
        back = 0; time = 0;
        for k = 1:nechoes
            [back,time] = progress(k,nechoes,'Performing ifft on echo',back,time);
            raw(:,:,:,k) = ifftnc(raw(:,:,:,k));
        end
    end
    img = raw; clear raw;
    
else    %%% GE Style
      
    load dims;
    load hdr hdr;

    echofiles = dir([scan_name '*']);
    nechofiles = length(echofiles);
    if length(dims) < 4
        dims(4) = 1;
        save dims dims
    end

    if dims(4) > nechofiles
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
        if exist('pswitch','var')
            raw = readMErawP(scan_name,'',pfile_type,1:length(dir([scan_name '*'])));
        else 
            raw = readMEraw(scan_name,'',pfile_type,1:length(dir([scan_name '*'])));
        end
    end
    
    % Multi-echo odd-even readout corrections
    if dims(4) > 1
        % Corrects odd-even phase shift
        if exist('phaseswitch','var')
            fprintf('  Performing odd/even phase correction in readout direction...\n');    
            raw = phaseME(raw); 
        end

        % Attempts to correct odd-even spatial shift/transform due to
        % misaligned TE. ACHTUNG! Uses fuzzy math and changes image intensity
        if exist('shiftswitch','var')
            fprintf('   Correcting even echo shift...\n'); 
            raw = shiftME(raw); 
        end
    end
    
    if ~rawswitch
        % Perform IFFT
        back = 0; time = 0;
        for k = 1:dims(4)
            [back,time] = progress(k,dims(4),'Performing ifft on echo',back,time);
            raw(:,:,:,k) = ifftnc(raw(:,:,:,k));
        end
    end
    img = raw; clear raw;
    
    if ~rawswitch
        figure(71);  
        subplot(211); montageME(abs(img),ceil(dims(3)/2));
        subplot(212); montageME(angle(img),ceil(dims(3)/2));

        % Centers image object in FOV
        if exist('centerswitch','var')
            img = centerME(img);
            figure(1); montageME(abs(img),ceil(dims(3)/2),'1D',dims(1:3));
        end
    end

end

if numel(dims) < 4
    dims(4) = 1;
    save dims dims;
end

clear raw

% Flip x, y, and/or z direction as specified
if exist('orientswitch','var')
    img = orientME(img,orientation);
    save orientation orientation
end

% Create image mask (works best for rodent brain)
if exist('maskswitch','var') 
    
    if exist('multimask_flag','var')
        [~] = brainextM(img);
        
    else
    
        if echoflag
           maskidx = 1;
        else 
           maskidx = TE < .020;
        end
        if exist('bruker_flag','var')
            [~] = brainext(img(:,:,:,maskidx),2,0.6);
    %         [~] = brainext(img(:,:,:,maskidx),2,0.4);
        else
            [~] = brainext(img(:,:,:,maskidx));
        end
    end 
end

% Calculate filtered phase and susceptibility maps
if exist('BIACswitch','var')
    fprintf('   Saving mask and phase information for BIAC computations...\n');
    load B0dir
    if dims(4) < 10
        iBIAC = runno*10+1;
    elseif dims(4) >= 10
        iBIAC = runno*100+1;
    end   
    iBIACs = iBIAC:iBIAC+dims(4)-1;
    save iBIACs iBIACs B0 runno
    BIACprep2(permute(img,[BIACpermute 4]),iBIACs(1), ...
              permute(open_nii('msk.nii'),BIACpermute));
    
    pts = prod(dims(1:3));
    if pts < 5e6
        vmem = 8;
    elseif pts < 20e6
        vmem = 16;
    elseif pts < 150e6
        vmem = 31;
    elseif pts < 500e6
        vmem = 46;
    else
        vmem = 92;
    end
    XBIAC3(iBIACs,vmem,30);
end

% Put image in proper orientation as specified
if exist('permuteswitch','var')
    fprintf('   Reorienting image according to your preference...\n');
    img = permute(img,permute_vector); 
    dims = size(img);
    save dims dims
    vox = vox(permute_vector(1:3));
    save vox vox
    if exist('maskswitch','var')
        msk = open_nii('msk.nii');
        msk = permute(msk,permute_vector);
        mat2nii(msk);
    end
    
    if numel(dims) < 4
        dims(4) = 1;
        save dims dims;
    end
end

% Save data as mat file 
if exist('matsave','var')
    fprintf('   Saving raw image data as img.mat...\n');
    img = single(img);
    save -V7.3 img img 
end

% Automatic SNR calculation
if exist('msk.nii','file') == 2
    if exist('bruker_flag','var')
        SNR = autoSNR(img,open_nii('msk.nii'),3,'gaussian','mean'); % mouse brains in surface coil only
    else
        SNR = autoSNR(img,open_nii('msk.nii'),3);
    end
end

% T2 map calculation
if dims(4) > 1   
    if exist('T2switch','var')
        if exist('msk.nii','file') == 2
%             T2 = T2mapMaskedWeighted(abs(img),TE,open_nii('msk.nii'),SNR(1,:));
            load meanSNR
            T2 = T2mapMaskedWeighted(abs(img),TE,open_nii('msk.nii'),meanSNR);
            writeNii(T2,'T2map.nii');
        end
    end
end
 
% retrieves filtered phase and susceptibility maps if they are complete
if exist('BIACswitch','var')
    if vmem == 8
        display('   Retrieving your susceptibility maps...please wait.');
        BIACpullwait 
    else
        display('   Susceptibility maps are likely still processing.');
    end
end

fprintf('   Data reconstruction complete!\n');
end
