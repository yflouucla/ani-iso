close all; clear; clc

I = im2double(imread('shape.png'));

g = fspecial('motion', 15,1);

sig = 0.1;
f = myconv(I,g)+sig*randn(size(I)); 


%% deblur
pm.mu = 5;
pm.lambda = 1;
uTV = deblurTV(f,g,pm);
uL1L2ap5 =deblurL1L2ap(f,g, .5,pm);
pm.maxDCA =2;
uL1L2ap =deblurL1L2ap(f,g, 1, pm);



% plot results, which are quantatively measured by SSIM
K = [0.05 0.05];
window = ones(8);
L = 1;


figure; 
subplot(222)
imshow(uTV); 
title(['L_1,SSIM=', num2str(ssim_index(uTV,I,K,window,L),3)])
subplot(221); imshow(f); title('input')
subplot(224); imshow(uL1L2ap5);
title(['L_1-0.5L_2, SSIM=',num2str(ssim_index(uL1L2ap5,I,K,window,L),3)])
subplot(223); imshow(uL1L2ap);
title(['L_1-L_2, SSIM=',num2str(ssim_index(uL1L2ap,I,K,window,L),3)])


