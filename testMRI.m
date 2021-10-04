clear; close all; clc

N = 256; 
L = 8; 

F = phantom(N);

Mask = fftshift(double(MRImask(N, L)));


data = Mask.*fft2(F)/N;

pm.mu = 1e3;
pm.lambda = 10;
pm.maxit = 1000;
pm.u_orig = F;

pmTV = pm;
pmTV.maxit = pm.maxit * 10;

tic
uTV = MRreconTV(Mask, data, pmTV);
pm.alpha = 1;
uL12ap1 = MRreconL1L2ap(Mask, data, pm);
pm.alpha = 0.5;
uL12ap5 = MRreconL1L2ap(Mask, data, pm);
toc


uFBP = abs(ifft2(data));

figure;
subplot(2,2,1);
imshow(uFBP,[]); colormap('gray');
title(['FBP, ' num2str(norm(uFBP-F, 'fro')/norm(F, 'fro'),4)]);
subplot(2,2,2);
imshow(uTV,[]); colormap('gray');
title(['TV, ', num2str(norm(abs(uTV)-F, 'fro')/norm(F, 'fro'),4)]);
subplot(2,2,3)
imshow(abs(uL12ap5),[]); colormap('gray');
title(['L_1-0.5L_2, ' num2str(norm(abs(uL12ap5)-F, 'fro')/norm(F, 'fro'),4)]);
subplot(2,2,4);
imshow(abs(uL12ap1),[]); colormap('gray');
title(['L_1-L_2, ' num2str(norm(abs(uL12ap1)-F, 'fro')/norm(F, 'fro'),4)]);
