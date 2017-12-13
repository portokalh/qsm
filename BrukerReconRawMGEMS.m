function [raw,bruker] = BrukerReconRawMGEMS
% "raw" is the raw data, "bruker" is the headfile
% Original file by Evan Calbrese, 6/8/11
% rmd - added multi-echo read, 6/8/11
% rmd - adapted for newest data and the PseriesP reconstruction, 8/6/12

% Read in Bruker parameters
bruker.acqp = readBrukerHeader('acqp');
bruker.method = readBrukerHeader('method');
raw = readBrukerFID('',bruker.method);
% try
%     bruker.subject = readBrukerHeader('../subject');
% catch
%     bruker.subject = readBrukerHeader('subject');
% end

if length(bruker.method.PVM_Matrix) < 3
 

        bruker.method.PVM_Matrix(1,3) = 1;
%     bruker.method.PVM_EncMatrix(1,3) = 1;

end

% Correct Bruker raw data with zero filling if array size is not 2^n
if prod([bruker.method.PVM_Matrix bruker.method.PVM_NEchoImages ...
         bruker.method.PVM_EncNReceivers]) ~= length(raw)
    [bruker,raw] = zerorid(bruker,raw);  
    
end

nechoes = bruker.method.PVM_NEchoImages;
dims = bruker.method.PVM_Matrix;


  
% else % works for 1-slice images only?

if strcmp(bruker.method.Method,'MDEFT')
    z = bruker.method.PVM_SPackArrNSlices;
    if prod([bruker.method.PVM_Matrix z ...
        bruker.method.PVM_EncNReceivers]) == length(raw)
        dims = bruker.method.PVM_Matrix;
    elseif prod([bruker.method.PVM_EncMatrix z ...
        bruker.method.PVM_EncNReceivers]) == length(raw)
        dims = bruker.method.PVM_EncMatrix;
    end
    dims(3) = z;
    
    
    img_idx = 1:bruker.acqp.ACQ_phase_factor*bruker.method.SegmNumber;
    img_idx = reshape(img_idx,bruker.method.SegmNumber,bruker.acqp.ACQ_phase_factor)';
    img_idx = img_idx(:);
    
    
    
        raw2 = reshape(raw,bruker.method.PVM_EncMatrix(1,1), ...
                           bruker.method.PVM_EncNReceivers,...
                           bruker.method.SegmNumber, ...
                               bruker.method.PVM_EncMatrix(1,2)/bruker.method.SegmNumber,...
                           z*nechoes ...
                           );
                       raw2 = reshape(permute(raw2,[1 3 4 5 2]),[220 160 56 2]);
                       raw2 = raw2(:,img_idx,:,:);
    
    
    
    
    
    
    raw2 = reshape(raw,bruker.method.PVM_EncMatrix(1,1), ...
                  bruker.method.PVM_EncNReceivers, ...
                  bruker.method.SegmNumber,z*nechoes,bruker.method.PVM_EncMatrix(1,2)/bruker.method.SegmNumber);
%   raw2 = permute(raw2,[1 5 3 4 2]); 
      raw2 = reshape(raw2,bruker.method.PVM_EncMatrix(1,1),...
                  bruker.method.PVM_EncMatrix(1,2), ...
                  z*nechoes,1,bruker.method.PVM_EncNReceivers);
  
  
  
%     raw2 = reshape(raw,bruker.method.PVM_EncMatrix(1,1),bruker.method.PVM_EncNReceivers, ...
%         bruker.method.PVM_EncMatrix(1,2), z);
% %     raw2 = raw2(:,:,img_idx,:);
%         raw2 = permute(raw2,[1 3 4 5 2]);

    raw2 = reshape(raw,bruker.method.PVM_EncMatrix(1,1), ...
        bruker.method.PVM_EncMatrix(1,2),bruker.method.PVM_EncNReceivers,z);
%     raw2 = raw2(:,:,img_idx,:);
        raw2 = permute(raw2,[1 2 4 3 5]);


elseif (strcmp(bruker.method.Method,'RARE') || (bruker.method.PVM_SPackArrNSlices > 1)) % multi-slice reconstruction, supports multi-echo
    z = bruker.method.PVM_Matrix(1,3)*bruker.acqp.NSLICES;
    rare = bruker.acqp.ACQ_rare_factor;

    if strcmp(bruker.method.PVM_SpatDimEnum,'2D')||strcmp(bruker.method.PVM_SpatDimEnum,'<2D>')
    raw = reshape(raw,bruker.method.PVM_EncMatrix(1,1), ...
                      bruker.method.PVM_EncNReceivers, ...
                      rare,z*nechoes,bruker.method.PVM_EncMatrix(1,2)/rare);
    raw = permute(raw,[1 5 3 4 2]); 
    raw = reshape(raw,bruker.method.PVM_EncMatrix(1,1),...
                  bruker.method.PVM_EncMatrix(1,2), ...
                  z*nechoes,1,bruker.method.PVM_EncNReceivers);
    raw(:,:,bruker.acqp.ACQ_obj_order+1,:,:) = raw;
    dims(2) = bruker.method.PVM_EncMatrix(1,2);
    raw = reshape(raw,[dims(1:2) nechoes z bruker.method.PVM_EncNReceivers]);
    raw = permute(raw,[1 2 4 3 5]);

    elseif strcmp(bruker.method.PVM_SpatDimEnum,'3D') % not sure if this block works
    raw = reshape(raw,bruker.method.PVM_EncMatrix(1,1),...
                      bruker.method.PVM_EncNReceivers, ...
                      rare,bruker.method.PVM_EncMatrix(1,2)/rare,z);
    raw = permute(raw,[1 4 3 5 2]);
    raw = reshape(raw,bruker.method.PVM_EncMatrix(1,1),... 
                  bruker.method.PVM_EncMatrix(1,2), ...
                  z,1,bruker.method.PVM_EncNReceivers);
    end
elseif bruker.method.PVM_NEchoImages > 1 % multi-echo recon
    if bruker.method.PVM_AntiAlias(1) == 1
        raw=reshape(raw,bruker.method.PVM_Matrix(1,1),  ...
        bruker.method.PVM_EncNReceivers,nechoes,  ...
        bruker.method.PVM_Matrix(1,2), bruker.method.PVM_Matrix(1,3));
        raw = permute(raw,[1 4 5 3 2]);
    %     raw = raw(:,:,:,1:2:end)+raw(:,:,:,2:2:end); 
        if strcmp(bruker.method.EchoAcqMode,'allEchoes') == 1
            raw(:,:,:,2:2:end,:) = raw(end:-1:1,:,:,2:2:end,:);   
        end
    else % reconstruct MGE with AntiAlias ~= 1
        raw=reshape(raw,bruker.method.PVM_Matrix(1,1)*bruker.method.PVM_AntiAlias(1),  ...
        bruker.method.PVM_EncNReceivers,nechoes,  ...
        bruker.method.PVM_Matrix(1,2), bruker.method.PVM_Matrix(1,3));
        raw = permute(raw,[1 4 5 3 2]);
        
        raw = ifftc(raw,[],1);
        raw = fftc(raw((size(raw,1)-dims(1))/2+1:(size(raw,1)+dims(1))/2,:,:,:));
                
    end 
 
else % Reconstruction for presumably all normally-acquired 2- or 3-D cartesian data
        raw = reshape(raw,bruker.method.PVM_EncMatrix(1,1),bruker.method.PVM_EncNReceivers, ...
        bruker.method.PVM_EncMatrix(1,2), bruker.method.PVM_Matrix(1,3));
        raw = permute(raw,[1 3 4 5 2]);
        
       

        

end

end





