------------------------------------------------------------------

Demo software for image denoising and deblurring examples in 

Y. Lou, T. Zeng, S. Osher, and J. Xin, "A Weighted Difference of 
Anisotropic and Isotropic Total Variation Model for Image Processing," 
submitted to SIAM J. Imaging Sci, September 2014.

------------------------------------------------------------------

Copyright (c) Yifei Lou 
https://sites.google.com/site/louyifei/
This work should be used only for nonprofit purposes.


------------------------------------------------------------------
Contents
------------------------------------------------------------------

The package comprises these functions

*) denoiseTV.m    	: TV denoising by split Bregman
*) denoiseL12ap.m      	: Ani-iso TV denoising by DCA+split Bregman
*) deblurTV.m		: TV deblurring by split Bregman
*) deblurL12ap.m   	: Ani-iso TV deblurring by DCA+split Bregman
*) MRreconTV.m		: MRI reconstruction by split Bregman
*) MRreconL1L2ap.m	: Ani-iso MRI reconstruction by DCA+split Bregman

*) testDenoise.m	: a demo code for denoising methods
*) testDeblur.m		: a demo code for deblurring methods
*) testMRI.m		: a demo code for MRI reconstruction methods

*) ssim_index.m		: SSIM to measure the performance
*) myconv.m 		: convolution code
*) MRImask.m		: generate radial lines for MRI reconstruction


------------------------------------------------------------------
Feedback
------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact Yifei Lou at: first.last@utdallas.edu

