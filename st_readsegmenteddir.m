function [ste, c] = st_readsegmenteddir(selpath)
%Read 10x folder
%Read files from a 10x Visium cellranger output folder
%readVisium

if nargin < 1, selpath = uigetdir; end
fprintf('Processing %s...\n', selpath);
[~, aff] = i_guessmtxfile(selpath);

if ~isempty(aff)
    h5fname = fullfile(selpath, sprintf('%sraw_feature_cell_matrix.h5', aff));
    zh5fname = fullfile(selpath, sprintf('%sraw_feature_cell_matrix.h5.gz', aff));
else
    h5fname = fullfile(selpath, 'raw_feature_cell_matrix.h5');
    zh5fname = fullfile(selpath, 'raw_feature_cell_matrix.h5.gz');
end


% h5fname=fullfile(selpath,'filtered_feature_bc_matrix.h5');
% h5fname=fullfile(selpath,'raw_feature_bc_matrix.h5');
if ~exist(h5fname, 'file')
    if ~exist(zh5fname, 'file')
        error('No raw_feature_cell_matrix.h5 file.');
    else
        [~, nametxt] = fileparts(zh5fname);
        fprintf('Unzipping %s.gz...\n', nametxt);
        gunzip(zh5fname);
    end
end

imgfolder = fullfile(selpath, 'spatial');
if ~isempty(aff)
    image_file = fullfile(imgfolder, sprintf('%stissue_hires_image.png', aff));
else
    image_file = fullfile(imgfolder, 'tissue_hires_image.png');
end
img = imread(image_file);

json_paths = fullfile(sprintf('%sgraphclust_annotated_cell_segmentations.geojson', aff));
jsonText = fileread(json_paths);
geoData = jsondecode(jsonText);

n = length(geoData.features);
xy = zeros(n, 2);
c = string(ones(n, 1));
for i = 1:n
    geometry = geoData.features(i).geometry;
    if strcmp(geometry.type, 'Polygon')
        s = squeeze(geometry.coordinates);
        [xy(i,1),xy(i,2)]=st.pkg.polycentroid(s);        
        if isfield(geoData.features(i).properties,'classification')
            c(i) = string(geoData.features(i).properties.classification.name);
        end
    end
end

% assignin("base","c",c);
[X, genelist, celllist] = sc_readhdf5file(h5fname);

assert(size(X, 2) == n);

sce = SingleCellExperiment(X, genelist);
sce.c_cell_id = celllist;
sce.c_cluster_id = c;
sce.c = findgroups(c);

metainfo = sprintf("Source: %s", selpath);
sce = sce.appendmetainfo(metainfo);


ste = SpatialTranscriptomicsExperiment(sce, xy, img);
% ste.tissue_positions_list = T;
% ste.scalefactors_json = scalef;
metainfo = sprintf("Source: %s", selpath);
ste = ste.appendmetainfo(metainfo);

    function [out, aff] = i_guessmtxfile(selpath)
        out = [];
        aff = [];
        a = dir(selpath);
        for k = 1:length(a)
            if contains(a(k).name, 'raw_feature_cell_matrix.h5')
                out = a(k).name;
                aff = extractBefore(out, 'raw_feature_cell_matrix.h5');
                continue;
            end
        end
end %  guess_mtx

end