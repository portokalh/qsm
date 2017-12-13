function raw = Pseries0(pfile,ee)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[hdr, header_bytes] = WrapGEheader(pfile(1).name);
dim1 = hdr.rdb.da_xres; % get raw dimensions from pfile header
dim2 = hdr.rdb.user7;
dim3 = hdr.rdb.user8;
dim4 = hdr.rdb.nechoes;
bp_half_element = hdr.rdb.point_size;
baseline_spacing = hdr.rdb.nframes;
SHORT_SIZE = 2;  % in bytes
flag='l';  % default byte swap
if bp_half_element == SHORT_SIZE
    pfile_type = 'short';
else
    if bp_half_element ~= 4
        error('\nPfile rdb.point_size unknown.\n');
    end
    pfile_type = 'int';
end

num_pfile = length(pfile);


v = 0;
raw = zeros(dim1,dim2,dim3);
for ipfile = 1:num_pfile

% Get p-file parameters
filename = pfile(ipfile).name;
[hdr, header_bytes] = WrapGEheader(filename);
fprintf('   %d/%d: %s',ipfile, num_pfile, filename);


fprintf(' 2x%d byte %s Pfile, nframes %d\n', bp_half_element, pfile_type, baseline_spacing);
% TE(ipfile) = hdr.rdb.te*1.0e-6;
filesize = dir(filename);
datapoints = (filesize.bytes-header_bytes)/bp_half_element;
sets_pe = datapoints/(baseline_spacing+1)/dim1/2/dim4; % Sets per echo
views_pe = sets_pe*(baseline_spacing);                 % Views per echo in this pfile

% Read p-file data into raw data array
fd=fopen(filename,'r',flag);
skip=fread(fd,[1 header_bytes/SHORT_SIZE],'short');  % header                  
skip=fread(fd,[2 dim1], pfile_type);                 % first baseline
back = 0; time = 0;
for s = 1:ceil(sets_pe)
    [back,time] = progress(s,ceil(sets_pe),'   Reading view set',back,time);
        vstart = v+1;
    for k = 1:dim4
        if s <= sets_pe
            for v = vstart:vstart+baseline_spacing-1               
                z = ceil(v/dim2);
                y = v-(z-1)*dim2;
                dump = fread(fd,[2 dim1], pfile_type);
                if ee == k
                    raw(:,y,z) = squeeze(dump(2,1:end)+1i*dump(1,:));
%                     eval(sprintf('raw%d(:,y,z) = squeeze(dump(2,1:end)+1i*dump(1,:));',ee));
                end
            end
        else
            for v = vstart:vstart+(views_pe-baseline_spacing*(s-1))-1; 
                z = ceil(v/dim2);
                y = v-(z-1)*dim2;
%                 if k ~= dim4
                dump=fread(fd,[2 dim1], pfile_type);
                if ee == k
                    raw(:,y,z) = squeeze(dump(2,1:end)+1i*dump(1,:));
%                     eval(sprintf('raw%d(:,y,z) = squeeze(dump(2,1:end)+1i*dump(1,:));',ee));
                end
%                 else
%                     dump=fread(fd,[2 dim1], pfile_type);
%                     raw(:,y,z,k)=squeeze(dump(2,1:end)+1i*dump(1,:));
%                 end
            end       
        end
    skip = fread(fd,[2 dim1], pfile_type);  % skip other baselines    
    end
end
% whatsleft = fread(fd, pfile_type);  % see if any data is left behind   

fclose(fd);
end %end pfile loop






end

