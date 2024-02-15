# A Dataset Containing Sentinel-2 Time Series of German National Forest Inventory Trees

[Download dataset (FIX ME)](https://zenodo.org)

This repo contains the code that was used to create the dataset published [here]. The dataset contains Sentinel-2 bottom of atmosphere reflectance time series for a subset of the trees measured in the 2012 and 2022 national forest inventory in Germany. It can be used to train machine learning and deep learning models to classify tree species from S2 time series.

The scripts provided here can only be used if you have access to the confidential NFI sampling locations. Consequently, they are mainly here to allow the public to check the validity of the code. If you need access to the NFI data, please contact the Thünen Institute directly and not the authors.

|ID|species|unix time|reflectance|quality flags|
|---|---|---|---|---|
|123|100|4568720|10 16 bit signed integers|16 bit unsigned int|

