he = imread('tissue_hires_image.png');
imshow(he), title('H&E image');
text(size(he, 2), size(he, 1)+15, ...
    'Image courtesy of Alan Partin, Johns Hopkins University', ...
    'FontSize', 7, 'HorizontalAlignment', 'right');

lab_he = rgb2lab(he);

ab = lab_he(:, :, 2:3);
ab = im2single(ab);
nColors = 3;
% repeat the clustering 3 times to avoid local minima
pixel_labels = imsegkmeans(ab, nColors, 'NumAttempts', 3);

imshow(pixel_labels, [])
title('Image Labeled by Cluster Index');

%%
figure;
for k = 1:3
    mask1 = pixel_labels == k;
    cluster1 = he .* uint8(mask1);
    subplot(2, 2, k)
    imshow(cluster1)
    title(sprintf('Objects in Cluster %d', k));
end

%%

L = lab_he(:, :, 1);
L_blue = L .* double(mask3);
L_blue = rescale(L_blue);
idx_light_blue = imbinarize(nonzeros(L_blue));