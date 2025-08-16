function M = readFeatureSlice(h5file, geneName)
if nargin<1, h5file='F:\Visium\test_data\Visium_HD_3prime_Human_Ovarian_Cancer_feature_slice.h5'; end
if nargin<2, geneName='MALAT1'; end

% readFeatureSlice  Load Visium HD feature slice for a specific gene.
%
%   M = readFeatureSlice(h5file, geneName)
%   - h5file: path to feature_slice.h5
%   - geneName: string, gene name to extract
%   - M: sparse matrix (spatial expression grid for that gene)
%
% Example:
%   M = readFeatureSlice('feature_slice.h5', 'MALAT1');
%   imagesc(log1p(full(M))); axis equal tight; colormap hot;

    %--- read gene list ---
    allgenes = h5read(h5file, '/features/name');
    if ischar(allgenes) % char array, convert to cellstr
        allgenes = cellstr(allgenes);
    end
    allgenes=strip(allgenes,char(0));
    idx = find(strcmp(allgenes, geneName));
    if isempty(idx)
        error('Gene %s not found in file.', geneName);
    end

    %--- construct dataset path ---    
    dspath = sprintf('/feature_slices/%d', idx-1); % some files are 0-based
    info = h5info(h5file, dspath);
    
    % Check structure: triplets (x,y,value) or dense matrix
    hasTriplets = any(strcmp({info.Datasets.Name}, 'row')) && ...
                  any(strcmp({info.Datasets.Name}, 'col')) && ...
                  any(strcmp({info.Datasets.Name}, 'data'));

    if hasTriplets
        x = h5read(h5file, [dspath '/row']);
        y = h5read(h5file, [dspath '/col']);
        v = h5read(h5file, [dspath '/data']);
        try
         M = sparse(double(y+1), double(x+1), double(v));
        catch ME
            disp(ME.message)
        end
    else
        % fallback: read entire dataset (dense grid)
        M = h5read(h5file, [dspath '/data']);
    end
end
