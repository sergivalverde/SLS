% *************************************************************************
% SAVE_COMPRESSED_NII --> compress and save a nifti image
% 
% uses the niftitoolbox - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
% uses standard unix tools (gzio, gunzip)
% 
% usage:    save_compressed_nii(NII_struct, path_to_saved_image);
%
%           NII_struct: nifti struct 
%           path_to_saved_image: path to save the image withouth extension
%
%
% Sergi Valverde 2015
% sergi.valverde@udg.edu
%
% *************************************************************************

function save_compressed_nii(nii_image, image_path)

    save_untouch_nii(nii_image, [image_path,'_tmp_.nii']); 
    system(['gzip ', image_path,'_tmp_.nii']);
    system(['mv ', image_path,'_tmp_.nii.gz ', image_path,'.nii.gz']);
    % if current image exists delete it
    if exist([image_path,'.nii']) == 2
        delete([image_path,'.nii']);
    end
end