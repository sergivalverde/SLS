
function [neighbor_voxels] = get_3D_neighbors(image_size, current_region, levels)
% ------------------------------------------------------------------------
% [neighbor_voxels] = get_3D_neighbors(image_size, current_region, levels) 
% 
% Compute the n-level neighbor voxels with 6-neighbor connectivity of a 
% given 3D  region from a nifti like matrix with shape similar 
% to (row,cols,slices). 
%
%  -image size: vector containing the dimension of the 3D matrix.
%  -current_region: Region voxels expressed as linear indices (as given by
%   the FIND function (region = FIND(3dmatrix == label))
%  -levels: number of neighbors: 1) first level, 2) first + and second level
%   (neighbors of the first level), and so on...
%
%  - neighbor_voxels: A matrix with indices as levels containing 
%    each one the neighbor voxels at each level as linear indices
%
% NOTES: 
%
% Strategy: Linear indices of region voxels are converted to subscripts 
%(x,y,z) and then diplaced +1,-1 position in x,y,z (6 connectivity 
% components). Then, all subscripts are converted back to linear indices 
% and repeated ones are remove generating a dilated mask of the input region. 
% Neighbor voxels are computed by substracting the dilated voxels and 
% the input region. 
%    
%  June 2015 Sergi Valverde
%  sergi.valverde@udg.edu
% ------------------------------------------------------------------------
  
    
    rows = image_size(1);
    cols = image_size(2); 
    slices = image_size(3);
    
    % iterate throught the neighbor levels. The dilated mask generated at
    % each particular level is used as an input for the next level. Only
    % unique voxel position as stored for each level. 
    % With levels > 1, works but slow, but makes sense if we want to obtain
    % information splitted by distance with respect to the lesion region.
    
    next_region = current_region;
    boundary_space = zeros(image_size);
    for l=1:levels
        %disp(['level :', num2str(l)]);
        [region_row, region_col, region_s] = ind2sub(image_size, next_region);

        rleft_cs = [max(region_row-1,1) , region_col, region_s];
        rright_cs = [min(region_row+1,rows), region_col, region_s];
        r_cleft_s = [region_row , max(region_col-1,1), region_s];
        r_cright_s = [region_row , min(region_col+1,cols), region_s];
        rc_sleft = [region_row , region_col, max(region_s-1,1)];
        rc_sright = [region_row , region_col, min(region_s+1,slices)];

        total_voxels = [rleft_cs; rright_cs; r_cleft_s; r_cright_s; rc_sleft; rc_sright];
        dilated_voxels = unique(sub2ind(image_size, total_voxels(:,1), total_voxels(:,2), total_voxels(:,3)));
    
        % For each level,remove region voxels and return just the boundary voxels enclosing
        % them. Maybe not so elegant :(
        
        boundary_space(dilated_voxels) = 1;
        current_region = next_region;
        next_region = find(boundary_space == 1);
        boundary_space(current_region) = 0;
        neighbor_voxels{l} = find(boundary_space == 1);
    end
end

