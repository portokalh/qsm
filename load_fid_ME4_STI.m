% function load_fid_ME4_STI(fidpath,max_blocks,ntraces,npoints,bitdepth,cyclenum,dims)

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
% data_buffer=zeros((npoints/2)*ntraces,max_blocks,'single');
% data_buffer=complex(data_buffer,data_buffer);

%fseek to the right place, skip 60 byte header and 28 byte block header, then data
byteskip=60+max_blocks*npoints*ntraces*bytes_per_point*(cyclenum-1)+28*(cyclenum-1)*(max_blocks);
fseek(fid,byteskip,'bof');

display('reading blocks');


nechoes = dims(4);
echodata = zeros(prod(dims(1:3)),nechoes,'single');

inx=1;%index pointer
for b = 1:max_blocks    
    
    data = fread(fid,npoints*ntraces,[bitdepth '=>single']); 
    fseek(fid,28,'cof'); %skip block header
    inx=inx+1;
    
    for t = 1:ntraces
        q = ntraces/nechoes;
        k = 1+mod(t-1,nechoes);
%         fwrite(fidvec(k),data((t-1)*npoints+1:t*npoints),'single');
        data2 = (t-1)*npoints+1:t*npoints;
        data2 = complex(data2(1:2:end),data2(2:2:end));
        echodata(1+(b-1)*ntraces+(t-1)*npoints/2:(b-1)*ntraces+t*npoints/2,k)
    end
    
    
end  % done reading one block  

fclose(fid);


% size(data_buffer)
% data_buffer= permute(reshape(data_buffer,[voldims(1) nechoes voldims(2) ...
%                      voldims(3)/nechoes]),[1 3 4 2]);%reshape into 4d array    

