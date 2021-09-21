# Neural-Network-Based Semantic Segmentation

This directory contains our implementation of a U-Net for segmenting the chiari data. 

## Requirements

This code requires recent versions of Python and pytorch as provided, e.g., on Google Colab. 

## Overview

[`MRI_segmentation.ipynb`]() - Loads the optimal model and segments given DENSE MR images outputting segmentation masks, DICE similarities, and loss.

[`U-Net_Emory_REU.ipynb`]() - Sets up the U-Net, trains a model (saving at lowest validation loss), and outputs segmentation masks, DICE similarities, and loss.
