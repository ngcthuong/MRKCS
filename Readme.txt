------------------------------------------------------------------------------------------

    Demo software for MULTI-SCALE/MULTI-RESOLUTION KRONECKER COMPRESSIVE IMAGING
	
              Public release ver. 1.1 (27 Jan. 2016)

------------------------------------------------------------------------------------------

The software reproduces the experiments published in the paper

 T. N. Canh, K. Q. Dinh and B. Jeon, “multi-scale/multi-resolution Kronecker compressive 
 imaging,�? Proc. of IEEE Inter. Conf. Image Process., pp. 2700 – 2704, Sep. 2015. 
 
 DOI http://dx.doi.org/10.1109/ICIP.2015.7351293

------------------------------------------------------------------------------------------

authors:               Thuong Nguyen Canh

web page:              https://sites.google.com/site/ngcthuong/

contact:               ngcthuong @ skku dot vn

------------------------------------------------------------------------------------------
Copyright (c) 2016 Sungkyunkwan University.
All rights reserved. 
This work should be used for nonprofit purposes only.
------------------------------------------------------------------------------------------

Contents
--------
demo_MRKCS_woPrior.m       	  	main script
Sensing/wavelet_matrix.m 		generate separable wavelet matrix
Sensing/measAlloc.m 			measurement allocation for MRKCS
Sensing/MRKCS_Sensing 			Generate MRKCS sensing matrix 
Utilities						Tools: tools, write output, IQA 
								
Requirements
------------
This demo is designed for Matlab for Windows (ver. 7.4 and above)
KCS_Rec_DTVNL_v1.1 is required, can be downloaded from the author website 

Description
-----------

The demo software reproduces the experiments published in the paper

 T. N. Canh, K. Q. Dinh and B. Jeon, Multi-scale/multi-resolution Kronecker compressive 
 imaging, Proc. of IEEE Inter. Conf. Image Process., pp. 2700 – 2704, Sep. 2015. 
 
 DOI http://dx.doi.org/10.1109/ICIP.2015.7351293

Runing
----------
1. Download and unpack KCS_Rec_DTVNL_v1.1 packet 
2. Run the scrip demo_MRKCS_woPior.m
3. Enjoy.  

Note
---------
1. Running time is depend on the size of test image
2. Performance migh vary due to the recovery

------------------------------------------------------------------------------------------

Disclaimer
----------

Any unauthorized use of these routines for industrial or profit-oriented activities is
expressively prohibited. By downloading and/or using any of these files, you implicitly
agree to all the terms of the TUT limited license (included in the file Legal_Notice.txt).
------------------------------------------------------------------------------------------