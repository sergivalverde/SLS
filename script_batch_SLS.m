% ************************************************************************
% SCRIPT for SLS BATCH (MATLAB implementation)
% June 2015 Sergi Valverde
% sergi.valverde@udg.edu
% 
% NOTES: 
%
%************************************************************************


%% options  
clear all;

HOME='/media/s/external/lesion_data'; %/sergi external 
image_folder = [HOME,'/VH/FLAIR'];

% SLS segmentation parameters
alpha = 2;          % standard deviation multiplier (t = \mu + \alpha*\sigma)
omega_t = 0.8;      % percentage of lesion voxels in either GM or WM 
omega_n = 0.6;      % percentage of neighbor lesion voxels in WM
min_size = 10;      % minimum number of voxels per lesion

%% read scans

dir_names = dir(image_folder);
for i=3:3%size(dir_names,1)
    
    % ********************************************************************
    % 1) load data 
    % - intensity normalized FLAIR image
    % - Brain Mask 
    % - SPM8 tissue segmentation
    % ********************************************************************
    current_scan = dir_names(i).name;
    flair_name = [image_folder,'/',current_scan,'/rmFLAIR'];
    bm_name = [image_folder,'/',current_scan,'/brainMask.betPope'];
    spm_name = [image_folder,'/',current_scan,'/spm8tissueseg.nii'];
    flair_scan = load_compressed_nii(flair_name);
    bms_scan = load_compressed_nii(bm_name);
    spm8_seg_scan = load_untouch_nii(spm_name);
    
    flair_img = flair_scan.img;
    brainmask = bms_scan.img;
    spm8_img = spm8_seg_scan.img;
    
    flair_img(~brainmask) = 0;
    spm8_img(~brainmask) = 0;
    
    disp(' ');
    disp('************************************************');
    disp(['Processing scan: ', current_scan ]);
    disp('************************************************');
    disp(' ');
    
   
    % ********************************************************************
    % 2) LESION SEGMENTATION 
    % ********************************************************************
    tic;
    lesion_mask = lesionSegmentationTool(flair_img, spm8_img, alpha, omega_t, omega_n, min_size);
    toc;
        
    % save the result
    flair_scan.img = lesion_mask;  
    save_compressed_nii(flair_scan, [image_folder,'/',current_scan,'/final_mask']);
    
end
    
    
    
    
        
