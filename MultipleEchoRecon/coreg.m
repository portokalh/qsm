function [ output_args ] = coreg(img,xform)

if ~exist(xform,'var')
    if ~exist(xform,'file')
        w = 5; count = 0;
        while ~exist(xform,'file')
            wstr = ['   Script has waited ' num2str(count,4) ...
                ' seconds for FSL transform matrix, ' xform '.'];
            fprintf('%s',wstr);
            pause(w)
            count = count + w;
            for k = 1:length(wstr)
                fprintf('\b');
            end
        end
        fprintf('   FSL transform matrix complete!\n');
    end
end

if exist(xform,'var')
    img(:,:,:,ipfile) = registerM(raw(:,:,:,ipfile),xform,vox);
else
    img(:,:,:,ipfile) = circshift(raw(:,:,:,ipfile),-echoshift);
end


end

