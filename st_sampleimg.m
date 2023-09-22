function [imgarray] = st_sampleimg(img, xy, sz)

if nargin < 3, sz = 10; end
% [n,m]=size(img,[1 2]);
N = size(xy, 1);
imgarray = cell(N, 1);
for k = 1:N
    imgarray{k} = imcrop(img, [xy(k, :) - round(sz/2), sz, sz]);
end
end
