function [ste] = st_read10xdir(selpath)
%Read 10x folder
%Read files from a 10x Visium cellranger output folder
%readVisium

if nargin < 1, selpath = uigetdir; end
fprintf('Processing %s...\n', selpath);
[~, aff] = i_guessmtxfile(selpath);

if ~isempty(aff)
    h5fname = fullfile(selpath, sprintf('%sfiltered_feature_bc_matrix.h5', aff));
    zh5fname = fullfile(selpath, sprintf('%sfiltered_feature_bc_matrix.h5.gz', aff));
else
    h5fname = fullfile(selpath, 'filtered_feature_bc_matrix.h5');
    zh5fname = fullfile(selpath, 'filtered_feature_bc_matrix.h5.gz');
end


% h5fname=fullfile(selpath,'filtered_feature_bc_matrix.h5');
% h5fname=fullfile(selpath,'raw_feature_bc_matrix.h5');
if ~exist(h5fname, 'file')
    if ~exist(zh5fname, 'file')
        error('[st_read10xdir] No matrix.h5 file.');
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

json_paths = fullfile(imgfolder, sprintf('%sscalefactors_json.json', aff));
txt = fileread(json_paths);
scalef = jsondecode(txt);
%     scalef.tissue_hires_scalef=value.tissue_hires_scalef;
%     scalef.spot_diameter_fullres=value.spot_diameter_fullres;
%     scalef.fiducial_diameter_fullres=value.fiducial_diameter_fullres;
%     scalef.tissue_lowres_scalef=value.tissue_lowres_scalef;


position_file = fullfile(imgfolder, sprintf('%stissue_positions_list.csv', aff));
if ~exist(position_file, 'file')
    position_file = fullfile(imgfolder, sprintf('%stissue_positions.csv', aff));
    T = readtable(position_file, 'ReadVariableNames', true);
else
    T = readtable(position_file, 'ReadVariableNames', false);
end


[X, genelist, celllist] = sc_readhdf5file(h5fname);
sce = SingleCellExperiment(X, genelist);
sce.c_cell_id = celllist;
metainfo = sprintf("Source: %s", selpath);
sce = sce.appendmetainfo(metainfo);

if ismember('barcode', T.Properties.VariableNames)
    assert(all(ismember(sce.c_cell_id, string(T.barcode))))
    [~, idx] = ismember(sce.c_cell_id, string(T.barcode));
else
    assert(all(ismember(sce.c_cell_id, string(T.Var1))))
    [~, idx] = ismember(sce.c_cell_id, string(T.Var1));
end


t = T(idx, :);
if ismember('pxl_row_in_fullres', T.Properties.VariableNames) && ...
        ismember('pxl_col_in_fullres', T.Properties.VariableNames)
    x_pixel = t.pxl_row_in_fullres;
    y_pixel = t.pxl_col_in_fullres;
else
    %s=[t.Var3,t.Var4];
    %ste.s=s;
    x_pixel = t.Var5;
    y_pixel = t.Var6;
end

xy = [x_pixel, y_pixel];
% sce.s=[x_pixel y_pixel];


xy = xy - mean(xy);
r = [range(xy(:, 1)), range(xy(:, 2))];
a = size(img, 1:2);
xy = xy .* (0.65 * a ./ r);
xy = xy + a / 2;

ste = SpatialTranscriptomicsExperiment(sce, xy, img);
ste.tissue_positions_list = T;
ste.scalefactors_json = scalef;
metainfo = sprintf("Source: %s", selpath);
ste = ste.appendmetainfo(metainfo);

    function [out, aff] = i_guessmtxfile(selpath)
        out = [];
        aff = [];
        a = dir(selpath);
        for k = 1:length(a)
            if contains(a(k).name, 'filtered_feature_bc_matrix.h5')
                out = a(k).name;
                aff = extractBefore(out, 'filtered_feature_bc_matrix.h5');
                continue;
            end
        end
end %  guess_mtx

end