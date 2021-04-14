%{
takes a 3D spin array (padded with zeros) and returns the sum of all the
layers added together scaled from 0 to 1.  Options to save img and show
img.
%}
function img = flattenSpins(spins, show, save, trial_dir)

[M, N, P] = size(spins);
img = zeros(M, N);

for i = 1:P
    img = spins(:,:,i) + img;
end

img_max = max(max(img));
img_min = min(min(img));
img = (img - img_min)./(img_max - img_min);

if show == true
    figure;
    imagesc(img)
    axis square
    colorbar
    colormap gray
    title('Flattened Spin Array')
end

if save == true
    imgName = strcat(trial_dir, 'flattenSpins.png');
    img8b = uint8(255 * mat2gray(img));
    imwrite(img8b, imgName);
end

end