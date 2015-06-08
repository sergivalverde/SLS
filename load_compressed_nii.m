function [nifti_out] = load_compressed_nii(image_path)
% *************************************************************************
% 
% Simple and stupid way to read compressed NIFTI scans in MATLAB. 
%
% Uses the niftitoolbox - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
% uses standard unix tools (gzio, gunzip)
% 
% usage:    [NII_struct] = load _compressed_nii(imag_path);
%
%           NII_struct: nifti struct 
%          
%
% Sergi Valverde 2015
% sergi.valverde@udg.edu
%
% *************************************************************************

%
%    
%  April 2015 Sergi Valverde
%  sergi.valverde@udg.edu
% ------------------------------------------------------------------------
      
     % generate a new identifier for the current image
     s = clock;
     tmp_id = num2str(round(rand*1000*s(6)));
     system(['cp ', image_path,'.nii.gz ' image_path,'_',tmp_id,'.nii.gz']);
     system(['gunzip ', image_path,'_',tmp_id,'.nii.gz']);
     nifti_out = load_untouch_nii([image_path,'_',tmp_id,'.nii']);
     system(['rm ', image_path,'_',tmp_id,'.nii']);
     
end