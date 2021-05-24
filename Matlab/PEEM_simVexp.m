clear all

OG = imread('SCO700eV_processed_HS.png');

randImg = imread('randSpinLatt.png');
zerosImg = imread('zeros.png');
onesImg = imread('ones.png');
meanOG = imread('mean.png');

J47 = imread('210422_J47K_T6.3191K_1_0.7071_0.5_S83.9_squeeze6LG.png');
J47pt5 = imread('210421_J47.5K_T6.2526K_1_0.7071_0.5_S83.9_squeeze6LG.png');
J48 = imread('210421_J48K_T6.1875K_1_0.7071_0.5_S83.9_squeeze6LG.png');
J48pt5 = imread('210421_J48.5K_T6.1237K_1_0.7071_0.5_S83.9_squeeze6LG.png');

%% Do Cross Correlation

cself = xcorr2(OG, OG);
cselfmax = max(cself(:))

czeros = xcorr2(OG, zerosImg);
czerosmax = max(czeros(:))

cones = xcorr2(OG, onesImg);
conesmax = max(cones(:))

cmean = xcorr2(OG, meanOG);
cmeanmax = max(cmean(:))

crand = xcorr2(OG, randImg);
crandmax = max(crand(:))

c47 = xcorr2(OG, J47);
c47pt5max = max(c47(:))

c47pt5 = xcorr2(OG, J47pt5);
c47pt5max = max(c47pt5(:))

c48 = xcorr2(OG, J48);
c48max = max(c48(:))

c48pt5 = xcorr2(OG, J48pt5);
c48pt5max = max(c48pt5(:))



%% Show cross correlations
figure;
subplot(2, 2, 1)
imshow(OG)
title('A = Original Image')
axis square

subplot(2, 2, 2)
imshow(OG)
title('B = Original Image')
axis square

subplot(2, 2, 3)
imagesc(cself)
axis square
title('Cross Correlation of A and B')

subplot(2, 2, 4)
plot(cself(:))
title('Cross Correlation of A and B: Max')
%%
figure;
subplot(2, 2, 1)
imshow(OG)
title('A = Original Image')
axis square

subplot(2, 2, 2)
imshow(zerosImg)
title('B = Zeros')
axis square

subplot(2, 2, 3)
imagesc(czeros)
axis equal
title('Cross Correlation of A & B')

subplot(2, 2, 4)
plot(czeros(:))
title('Cross Correlation of A & B: Max')
%%
figure;
subplot(2, 2, 1)
imshow(OG)
title('A = Original Image')
axis square

subplot(2, 2, 2)
imshow(onesImg)
title('B = ones')
axis square

subplot(2, 2, 3)
imagesc(cones)
axis square
title('Cross Correlation of A & B')

subplot(2, 2, 4)
plot(cones(:))
title('Cross Correlation of A & B')
axis square
%%
figure;
subplot(2, 2, 1)
imshow(OG)
title('A = Original Image')
axis square

subplot(2, 2, 2)
imshow(meanOG)
title('B = Original Image Mean')
axis square

subplot(2, 2, 3)
imagesc(cmean)
axis square
title('Cross Correlation of A & B')

subplot(2, 2, 4)
plot(cmean(:))
axis square
title('Cross Correlation of A & B: Max')

%%
figure;

subplot(2, 2, 1)
imshow(OG)
axis square
title('A = Original Image')

subplot(2, 2, 2)
imshow(randImg)
axis square
title('B = Random Image')

subplot(2, 2, 3)
imagesc(crand)
axis square
title('Cross Correlation of A & B')

subplot(2, 2, 4)
plot(crand(:))
axis square
title('Cross Correlation of A & B')

%%

figure;

subplot(2, 2, 1)
imshow(OG)
axis square
title('A = Original Image')

subplot(2, 2, 2)
imshow(J47)
axis square
title('B = J47 K')

subplot(2, 2, 3)
imagesc(c47)
axis square
title('Cross Correlation A & B')

subplot(2, 2, 4)
plot(c47(:))
axis square
title('Cross Correlation A & B: Max')

%%

figure;

subplot(2, 2, 1)
imshow(OG)
axis square
title('A = Original Image')

subplot(2, 2, 2)
imshow(J47pt5)
axis square
title('B = J47.5 K')

subplot(2, 2, 3)
imagesc(c47pt5)
axis square
title('Cross Correlation A & B')

subplot(2, 2, 4)
plot(c47pt5(:))
axis square
title('Cross Correlation A & B: Max')


%%
figure;
subplot(2, 2, 1)
imshow(OG)
axis square
title('A = Original Image')

subplot(2, 2, 2)
imshow(J48)
axis square
title('B = J48 K')

subplot(2, 2, 3)
imagesc(c48)
axis square
title('Cross Correlation A & B')

subplot(2, 2, 4)
plot(c48(:))
axis square
title('Cross Correlation A & B: Max')

%%

figure;

subplot(2, 2, 1)
imshow(OG)
axis square
title('A = Original Image')

subplot(2, 2, 2)
imshow(J48pt5)
axis square
title('B = J48.5 K')

subplot(2, 2, 3)
imagesc(c48pt5)
axis square
title('Cross Correlation A & B')

subplot(2, 2, 4)
plot(c48pt5(:))
axis square
title('Cross Correlation A & B: Max')



%% Do SSIM

[svalself, smapself] = ssim(OG, OG);
[svalzeros, smapzeros] = ssim(zerosImg, OG);
[svalones, smapones] = ssim(onesImg, OG);
[svalmean, smapmean] = ssim(meanOG, OG);
[svalrand, smaprand] = ssim(randImg, OG);

[sval47pt5, smap47pt5] = ssim(J47pt5, OG);
[sval48, smap48] = ssim(J48, OG);
[sval48pt5, smap48pt5] = ssim(J48pt5, OG);

%% Show SSIM

figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(OG)
axis square
title('B: Original Image')

subplot(1, 3, 3)
imagesc(smapself)
axis square
title('SSIM A & B')
axis off

%%
figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(zeros)
axis square
title('B: Zeros')

subplot(1, 3, 3)
imagesc(smapzeros)
axis square
title('SSIM A & B')
axis off


%%

figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(onesImg)
axis square
title('B: Ones')

subplot(1, 3, 3)
imagesc(smapones)
axis square
title('SSIM A & B')
axis off

%%
figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(meanOG)
axis square
title('B: Original Image Mean')

subplot(1, 3, 3)
imagesc(smapmean)
axis square
axis off
title('SSIM A & B')
%%
figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(randImg)
axis square
title('B: Random')

subplot(1, 3, 3)
imagesc(smaprand)
axis square
axis off
title('SSIM A & B')

%%
figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(J47pt5)
axis square
title('B: J = 47.5 K')

subplot(1, 3, 3)
imagesc(smap47pt5)
title('SSIM A & B')
axis square
axis off

%%

figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(J48)
axis square
title('B: J = 48 K')

subplot(1, 3, 3)
imagesc(smap48)
title('SSIM A & B')
axis square
axis off

%%

figure;
subplot(1, 3, 1)
imshow(OG)
axis square
title('A: Original Image')

subplot(1, 3, 2)
imshow(J48pt5)
axis square
title('B: J = 48.5 K')

subplot(1, 3, 3)
imagesc(smap48pt5)
title('SSIM A & B')
axis square
axis off


%% Do 2DFFT

fftOG = fft2(OG);
fft47pt5 = fft2(J47pt5);
fft48 = fft2(J48);
fft48pt5 = fft2(J48pt5);

%% Show 2D FFT
figure
imagesc(real(fftOG))
set(gca, 'ColorScale', 'log')
title('2D FFT of Original Image')
axis equal

figure
imagesc(abs(fft47pt5))
set(gca, 'ColorScale', 'log')
title('2D FFT J = 47.5 K')
axis equal

figure
imagesc(abs(fft48))
set(gca, 'ColorScale', 'log')
title('2D FFT J = 48 K')
axis equal

figure
imagesc(abs(fft48pt5))
set(gca, 'ColorScale', 'log')
title('2D FFT J = 48.5 K')
axis equal

%% Do Cross Correlation on the 2D FFT imgs



