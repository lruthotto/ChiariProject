# Atlas-Based Segmentation

This folder contains the MATLAB files for running the atlas-based segmentations.

## Pre-Requisites

The histogram equalization uses the MATLAB image processing toolbox and the registration component requires [FAIR.m](https://github.com/C4IR/FAIR.m) (we used commit 975edeb).

## Overview

The main function is contained in the file 'chiari_atlas.m'. 

This function makes use of the following helper routines

- `biomarker.m` - Given a mask and DENSE image, produce the spatial-average-temporal-peak biomarker
- `chiari_example.m` - perform pairwise registration
- `chiari_example_average.m` - average many pairwise registrations
- `dice_jaccard.m` - quantify the segmentation accuracy
- `generate_avg.m` - 
- `histnormal.m` - histogram normalization
- `viewContour2D.m` - visualizes masks as contour plot
