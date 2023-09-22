function [ste] = st_sce2ste(sce, xyinfo, imginfo)

narginchk(2, 3);
validateattributes(sce, {'SingleCellExperiment'}, {});
if nargin < 3, imginfo = []; end
ste = [];

%switch typeid
%    case 1    % direct xy
%        ste=SpatialTranscriptomicsExperiment(sce,xyinfo);
%    case 2    % xy location file
if ~exist(xyinfo, 'file')
    [fname, pathname] = uigetfile( ...
        {'*.csv', 'MTX Format Files (*.csv)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a xy location file');
    if isequal(fname, 0), return; end
    xyinfo = fullfile(pathname, fname);
end
if exist(xyinfo, 'file')
    t = readtable(xyinfo, 'ReadVariableNames', false);
    s = string(t.(t.Properties.VariableNames{1}));
    xy = [t.(t.Properties.VariableNames{2}), ...
        t.(t.Properties.VariableNames{3})];
    [a, xi, xj] = intersect(sce.c_cell_id, s);
    sce = sce.selectcells(xi);
    assert(isequal(a, sce.c_cell_id));
    xy = xy(xj, :);
    ste = SpatialTranscriptomicsExperiment(sce, xy);
end
if exist(imginfo, 'file')
    img = imread(imginfo);
    ste.img = img;
end
%    case 3
%end
