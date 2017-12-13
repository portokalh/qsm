function [phs,scl] = CleanFreq(phase_img,mag_img,varargin)
%CleanFreq.m: Cleans up that messy noise in MapFreq images
%
% CleanFreq(phase_img,mag_img,matrix_size) if phase_img is string.  If you
% are lazy and your data is in the form of m x m x m (m^3), don't enter
% matrix size.  If you are even lazier, and phase_img and mag_img have the
% same base name, you only have to enter one input.
%
% CleanFreq(phase_img,mag_img,'echo#') if phase_img and mag_img are already  
% variable matrices.  Include 'echo#' as a string if working with a single 
% echo so that it is properly labeled.

if ischar(phase_img) == 1

    fid = fopen(phase_img,'r','l');
    freq = fread(fid,inf,'float');
    fclose(fid);
    echo_only = phase_img(1:end-5);
    
    if nargin < 3
        
        if nargin < 2
            mag_img = [echo_only '.mag'];
        end      
        
        if ischar(mag_img) == 1       
            res = (round((length(freq))^(1/3)));
            if (res^3 ~= length(freq))
                res = sqrt(length(freq));
                matrix_size = [res res 1];
                disp(['Assuming 2D isotropic data, ' num2str(res) '^2']);        
            else
                matrix_size = [res res res];
                disp(['Assuming 3D isotropic data, ' num2str(res) '^3']);
            end
        else
            matrix_size = mag_img;
            mag_img = [echo_only '.mag'];
        end

    else
        matrix_size = varargin{1};                              
    end
      
    fid = fopen(mag_img,'r','l');
    mag = fread(fid,inf,'float');
    fclose(fid);
    
    voxels = 1;
    for k = 1:length(matrix_size(:))
        voxels = voxels*matrix_size(k);
    end
    
    freq = reshape(freq(1:voxels),matrix_size);
    mag = reshape(mag(1:voxels),matrix_size);
else
    mag = mag_img;
    freq = phase_img;
    if size(mag,4) == 1
        echo_only = ['echo' varargin{1}];
    end
end

noise_frame = mag(1:size(mag,1),1:ceil(size(mag,2)/64),1:ceil(size(mag,3)/64));
noise_mean = mean(noise_frame(:));
noise_std = std(noise_frame(:));

% Pick threshold criteria
% magmax = max(mag(:));
% thresh = .05*magmax               % cutoff fraction of magmax
thresh = noise_mean+6*noise_std;    % 6-sigma quality!
Required_SNR = thresh/noise_std;
disp(['  SNR cutoff used for clean up: ' num2str(Required_SNR)]);

% Creates threshold max from magnitude data
mask = mag > thresh;

% Creates phase image with zero (gray) background (phase +/-)
phs = freq.*mask;

% Creates phase image mask for SWI
scl = phs;
scl(scl >= 0) = 1;
scl(scl < 0) = scl(scl < 0)/180+1; % pi radians in degrees
scl = scl.^4; % 4th power mapping

center_img = ceil(size(mag,3)/2);
figure(85)
subplot(221); show(mag,center_img); title('Magnitude');
subplot(222); show(freq,center_img); title('Frequency Map');
subplot(223); show(phs,center_img); title('Cleaned');
subplot(224); show(scl,center_img); title('Scaled');

for k = 1:size(mag,4) 
    writeME(phs,k);
    writeME(scl,k);
end

end

function [] = show(image,varargin)
%show.m: easier image display

if nargin > 1
    imagesc(image(:,:,varargin{1})); colormap gray; axis image; axis off;
else
    imagesc(image); colormap gray; axis image; axis off;
end

end

