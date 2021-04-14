function img_bin = binarizeSpins(img, levels, show, save, trial_dir)

%blur image to remove thermal noise
img8b = uint8(255 * mat2gray(img));
img = filter2(fspecial('average',20),img8b/255);
figure
imshow(img)

if levels > 1
    thresholds = multithresh(img, levels);
    img_bin = imquantize(img, thresholds);
else
    img_bin = imbinarize(img, 'adaptive', 'Sensitivity', 0.7);
end

%remove noise


if show == true
    figure;
    imshow(img_bin, [])
    axis square
    colorbar
    colormap gray
    title('Thresholded Spin Array')
end

if save == true
    imgName = strcat(trial_dir, 'flattenSpins.png');
    img8b = uint8(255 * mat2gray(img));
    imwrite(img8b, imgName);
end

end