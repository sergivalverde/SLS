# SLS (Salem lesion segmentation)

## About:

An implementation of the SLS toolbox originally developed in C++ into Matlab. 


## Performance tests:

We have evaluated the performance of the Matlab implementation against the original C++ code. Comparision tests have been run using the set of 70 images (3T) provided by the Vall d'Hebron hospital center. All the images have been skull-stripped, intensity-corrected and denoised before. The same image masks are employed in either the Matlab and C++ SLS implementations.

First, we evaluated the Dice coefficient (DSC)  between the final lesion masks obtained with both SLS implementations. The DSC values between implementations were higher than 0.83 in all the images. Observed differences were due to small (*<0.01*) changes in the *sigma* estimation.  

![SLS-MSLS Dice coefficients](DSC_between_masks.png)


However, differences between masks were not relevant in terms of the accuracy of each implementation when compared with the provided ground-truth annotations. We compared the DSC coefficient of each SLS implementation against the manually annotated masks (GT). The performance of both SLS implementation was virtually identical.


![SLS-MSLS Dice coefficients](DSC_between_gt.png)


## Usage:

Find attached the ```test_SLS.m``` script to test the toolbox with one example image containing a MS patient image with the necessary *FLAIR* and the *T1* 3-class tissue segmentation.


```
clear all;

image_folder = 'test_data';

% SLS segmentation parameters
alpha = 2;          % standard deviation multiplier (t = \mu + \alpha*\sigma)
omega_t = 0.8;      % percentage of lesion voxels in either GM or WM
omega_n = 0.6;      % percentage of neighbor lesion voxels in WM
min_size = 10;      % minimum number of voxels per lesion

%% read scans

% ********************************************************************
% 1) load data
% - intensity normalized FLAIR image
% - Brain Mask
% - SPM8 tissue segmentation
% ********************************************************************
flair_name = [image_folder,'/FLAIR'];
t1_segmentation = [image_folder,'/T1segmentation'];

flair_scan = load_compressed_nii(flair_name);
seg_scan = load_compressed_nii(t1_segmentation);

flair_img = flair_scan.img;
brainmask = flair_img > 0;
t1seg_img = seg_scan.img;

flair_img(~brainmask) = 0;
t1seg_img(~brainmask) = 0;


% ********************************************************************
% 2) LESION SEGMENTATION
% ********************************************************************
tic;
[lesion_mask] = lesionSegmentationTool(flair_img, t1seg_img, alpha, omega_t, omega_n, min_size);
toc;

% store the result
flair_scan.img = lesion_mask;
save_compressed_nii(flair_scan, [image_folder,'/SLS_mask']);

```

## TODO:
+ Run a series of tests to show the performance between the C++ and Matlab implementations. 

