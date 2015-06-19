# SLS (Salem lesion segmentation)

## About:

We implement the SLS toolbox originally developed in C++ into Matlab. This is a work in progress.

## Differences:

+ So far, the histogram used in the __compute_fwhm.m__ function is different, and we obtain a different sigma threshold.
+ When fixing the same threshold, the number of estimated lesions after connected-components, and after applying the rules 1 and 2 (*min_size* and *omega_t* parameters) is identical.


## TODO:
+ test the differences in the histogram (compute_fwmh.m).
+ Rule 3 (omega_n) appears to behave differently.

