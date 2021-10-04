clear; close all; 
clc

I = im2double(imread('lena512.png'));
I = imresize(I,0.5);
sig = 0.05;
f = imnoise(I, 'gaussian',0,sig^2);
 
pm.mu = 20;

% Denoise
uTV = denoiseTV(f,pm);
uL1L2ap5 =denoiseL1L2ap(f, 0.5, pm);
pm.maxDCA = 2;
uL1L2ap =denoiseL1L2ap(f, 1, pm); % in the paper, we stop at two outer iterations


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


