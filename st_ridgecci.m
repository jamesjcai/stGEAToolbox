function [T, cx] = st_ridgecci(ste, group, idx1, idx2, K, methodid)

if nargin < 6, methodid = 1; end
if nargin < 5, K = 5; end

T = [];
if nargin < 4, idx2 = 2; end % target cell group id (receptor gene expressed)
if nargin < 3, idx1 = 1; end % source cell group id (ligand gene expressed)
% if nargin<2, [group]=st_cluster_img(ste,'method','resnet18kmeans','k',5); end


[a, b] = cdgea;
load(fullfile(a, 'resources', 'Ligand_Receptor.mat'), ...
    'ligand', 'receptor');
cd(b);
valididx = ismember(ligand, ste.sce.g) & ismember(receptor, ste.sce.g);
lrdatabase = [ligand(valididx), receptor(valididx)];

% K=gui.i_inputnumk(5);
% if isempty(K)
%     error('K=5');
% end

[A] = logical(full(sc_knngraph(ste.xy, K)));
%X1=ste.sce.X(:,group==idx1);
%X2=ste.sce.X(:,group==idx2);
MASK = logical(group == idx1*(group == idx2)');
MASK1 = MASK & A;
MASK2 = MASK & (~MASK1);

if nargout > 1
    cx1 = zeros(size(group));
    [i1, i2] = find(MASK1);
    cx1(i1) = 1;
    cx1(i2) = 2;

    cx2 = zeros(size(group));
    [i1, i2] = find(MASK2);
    cx2(i1) = 1;
    cx2(i2) = 2;
    cx = {cx1, cx2};
end

N = size(lrdatabase, 1);
% N=480;

pval = ones(N, 1);
kstat = zeros(N, 1);
numneigb = zeros(N, 1);
numouter = zeros(N, 1);
lgene = lrdatabase(1:N, 1);
rgene = lrdatabase(1:N, 2);
avgneigb = zeros(N, 1);
avgouter = zeros(N, 1);
stdneigb = zeros(N, 1);
stdouter = zeros(N, 1);

%RES=cell(N,1);

[i1, i2] = find(MASK1);

fw = gui.gui_waitbar_adv;
for xx = 1:N
    gui.gui_waitbar_adv(fw, xx/N);
    [~, lidx] = ismember(lrdatabase(xx, 1), ste.sce.g);
    [~, ridx] = ismember(lrdatabase(xx, 2), ste.sce.g);
    x1 = ste.sce.X(lidx, :); % ligand
    x2 = ste.sce.X(ridx, :); % receptor

    switch methodid
        case 1
            X12 = x1' * x2;
            nv = X12(MASK1);
            vv = X12(MASK2);
        case 2
            x1 = x1(i1)';
            x2 = x2(i2)';
            nv = x1 .* x2;
            x2 = x2(randperm(size(x2, 1)));
            vv = x1 .* x2;
    end

    numneigb(xx) = length(nv);
    numouter(xx) = length(vv);
    avgneigb(xx) = mean(nv);
    avgouter(xx) = mean(vv);
    stdneigb(xx) = std(nv);
    stdouter(xx) = std(vv);
    [~, pval(xx), kstat(xx)] = kstest2(nv, vv);
    %RES{xx}={nv,vv};
end
gui.gui_waitbar_adv(fw);

[~, ~, ~, fdr] = pkg.fdr_bh(pval);
T = table(lgene, rgene, numneigb, numouter, kstat, pval, fdr, ...
    avgneigb, avgouter, stdneigb, stdouter);
[T] = sortrows(T, 'kstat', 'descend');
end
