function [c, T, S] = i_ctypescore(X, g)

stag = 'hs';
markerfile = sprintf('marker_%s.mat', stag);

pw1 = fileparts(mfilename('fullpath'));
markerfile = fullfile(pw1, '..', 'scGEAToolbox', ...
    '+run', 'thirdparty', 'alona_panglaodb', markerfile);

if exist(markerfile, 'file')
    load(markerfile, 'Tw', 'Tm');
else
    error('xxx');
end
disp('Marker database loaded.');

tic;
S = zeros(size(X, 2), 10);
ctlist = string(Tm.Var1);
for idx = 1:length(ctlist)
    ctmarkers = Tm.Var2{idx};
    posg = string(strsplit(ctmarkers, ','));
    posg(strlength(posg) == 0) = [];
    S(:, idx) = sc_cellscore(X, g, posg);
end
toc;
[~, id] = max(S, [], 2);
c = ctlist(id);
t = tabulate(c);
T = cell2table(t, 'VariableNames', ...
    {'Value', 'Count', 'Percent'});
T = sortrows(T, 'Count', 'descend');
