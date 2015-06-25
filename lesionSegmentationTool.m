function [lesion_mask] = lesionSegmentationTool(input_image, tissueProb, alpha, omega_t, omega_n, min_size)
% ------------------------------------------------------------------------
% [lesion_mask] = lesionSegmentationTool(input_image, tissueProb, alpha, omega_t, omega_n, min_size) 
% 
% WM lesion segmentation of MS lesions on FLAIR images. 
%
%
%  -input_image  --> Pre-processed FLAIR (skull-stripped + intensity corrected) 
%  -tissueProb   --> labeled T1-w tissue segmentation (1 = CSF, 2 = GM, 3 = WM)
%  -alpha:       --> Parameter for controlling the weight of the sigma computed from 
%                    [mu, sigma] = FWHM function.
%  - omega_t     --> Percentage of either GM / WM for candidate lesion regions
%  - omega_n     --> Percentage of WM voxels for 6-connectivity 3D boundary
%                    voxels of candidate_regions.
%   -min_size    --> minimum area of candidate lesion regions.
%
%  - lesion_mask  --> Binary mask of the final candidate lesion regions.
%
%  June 2015 Sergi Valverde / Eloy Roura 
%  sergi.valverde@udg.edu
% ------------------------------------------------------------------------

    % vars:
    lesion_mask = zeros(size(input_image)); % output lesionmask
    gm_mask = tissueProb == 2;              % GM TISSUE
    wm_mask = tissueProb == 3;              % WM TISSUE

    
    % ********************************************************************
    % 1. Threshold evaluation
    %
    % - by default: 512 bins are computed. The final threshold is computed
    %   following T = mu + alpha * sigma.
    % - So far, the FWHM methods are different between implementations :(
    %
    % ********************************************************************
    
    gm_flair = input_image(gm_mask == 1);
    [mu, sigma] = compute_fwhm(gm_flair,512);
    T = mu + (alpha * sigma);
    thresholded_mask = input_image >= T;
    disp(['Threshold: ', num2str(T), ' (',num2str(mu),' + ',num2str(alpha),' * ',num2str(sigma),')']);
    
    
    % ********************************************************************
    % 2. Filtering step
    %
    % - Connected components are computed using 6-neighbor connectivity in
    %   3D. (followin c++ implementation)
    % ********************************************************************
    
    % Find connected components in 3D
    CC = bwconncomp(thresholded_mask,6);
    candidate_regions = labelmatrix(CC);
    candidate_labels = unique(nonzeros(candidate_regions));
    
    disp(['Number of lesions before refinement: ', num2str(CC.NumObjects)]);
    
    % RULE 1). 3d lesion size has to be higher than min_size parameter
    labels_filter_1 = cellfun(@(x) numel(x) < min_size, CC.PixelIdxList);
    candidate_labels(labels_filter_1) = 0;
    
    % RULE 2). omega_t percentage of lesion voxels either GM and WM 
    labels_filter_2 = cellfun(@(x) ((sum(gm_mask(x)) + sum(wm_mask(x))) / numel(x)) < omega_t, CC.PixelIdxList);
    candidate_labels(labels_filter_2) = 0;  
  
    
    % RULE 3) Omega_n Percentage of neighbors as WM
    % so far, have to visit each lesion candidate with a for loop :(
    candidate_labels = nonzeros(candidate_labels);
    num_lesions = 0;
    for l=candidate_labels'
            % current label 
            region_voxels = find(candidate_regions == l);
            
            % default config for neighbors: 24 connectivity in 3D (radius 1)
            % Neighbor voxels for each region are found by substracting the
            % total neighbors of the candidate region and the same
            % candidate region. Not shure if efficient :/
            boundary_region_voxels = compute_neighborhoods(region_voxels, size(candidate_regions), 1,3);
            candidate_regions(boundary_region_voxels) = l;
            candidate_regions(region_voxels) = 0;
            boundary_voxels = find(candidate_regions == l);
            
            %boundary_voxels = boundary_region_voxels{1};
            perc_wm_elem = sum(wm_mask(boundary_voxels)) / numel(boundary_voxels);
            if perc_wm_elem >= omega_n
                lesion_mask(region_voxels) = 1;
                num_lesions = num_lesions +1;
            end
    end
    disp(['Final number of lesions: ', num2str(num_lesions)]);
    
end