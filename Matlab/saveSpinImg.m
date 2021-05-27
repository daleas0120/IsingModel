function saveSpinImg(img, imgName)
fig = figure;
ah = axes('Units', 'Normalize', 'Position', [0 0 1 1]);
imagesc(img)
axis square
axis off
caxis([-1 1])
fig.Position(3) = fig.Position(4);
saveas(gcf, imgName)


end