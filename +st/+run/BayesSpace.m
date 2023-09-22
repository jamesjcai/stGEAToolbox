function [T, X, ste1] = BayesSpace(ste, markers)

if nargin < 2, markers = []; end
if ~isempty(markers)
    [y, idx] = ismember(markers, ste.sce.g);
    assert(all(y));
else
    idx = [];
end

T = [];
X = [];
isdebug = false;

rscriptdir = 'R_BayesSpace';
oldpth = pwd();
[isok, msg] = commoncheck_R(rscriptdir);
if ~isok, error(msg); end

tmpfilelist = {'input.txt', 'input.mat', 'positions.csv', ...
    'output.h5', 'positions_enhanced.csv', 'Rplots.pdf'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

% writematrix(ste.sce.X,'input.csv');
% sc_writefile('input.txt',ste.sce.X, ste.sce.g);

X = ste.sce.X;
save('input.mat', 'X', 'idx', '-v7.3')
writepositions(ste, 'positions.csv');

Rpath = getpref('scgeatoolbox', 'rexecutablepath');
pkg.RunRcode('script.R', Rpath);

if ~exist('output.h5', 'file'), return; end
web('Rplots.pdf', '-browser');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
T = readtable("positions_enhanced.csv");
warning('on');

if ~isempty(idx)
    X = h5read('output.h5', '/X');
end
if size(X, 1) ~= length(idx)
    X = [];
end

if nargout > 2 && ~isempty(X) && ~isempty(markers)
    sce1 = SingleCellExperiment(X, markers);
    sce1.struct_cell_clusterings.bayesspace = T.spatial_cluster;
    ste1 = SpatialTranscriptomicsExperiment(sce1, [T.row, T.col]);
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
