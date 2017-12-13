function data_buffer = load_fid_ME(fidpath,max_blocks,ntraces,npoints,bitdepth,cyclenum,voldims,nechoes)

try
    fid = fopen(fidpath,'r','ieee-be');
catch ME
    disp(ME)
end


%find out how many bytes per point
if strcmp(bitdepth,'int16');
    bytes_per_point=2;
else
    bytes_per_point=4;
end

%preallocate complex array
display('preallocating complex array');
data_buffer=zeros((npoints/2)*ntraces,max_blocks,cyclenum,'single');
data_buffer=complex(data_buffer,data_buffer);

%fseek to the right place, skip 60 byte header and 28 byte block header, then data
byteskip=60+max_blocks*npoints*ntraces*bytes_per_point*(cyclenum-1)+28*(cyclenum-1)*(max_blocks);
fseek(fid,byteskip,'bof');

display('reading blocks');

inx=1;%index pointer
for b = 1:max_blocks    
    data = fread(fid,npoints*ntraces,[bitdepth '=>single']); 
    data_buffer(:,inx)=complex(data(1:2:end),data(2:2:end));
    fseek(fid,28,'cof'); %skip block header
    inx=inx+1;
end  % done reading one block  

fclose(fid);
 
data_buffer= permute(reshape(data_buffer,[voldims(1) nechoes voldims(2:3)]), ...
             [1 3 4 2]);%reshape into 4d array    

