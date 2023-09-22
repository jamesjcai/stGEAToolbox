function [T1, T2] = st_ridgenet(ste, group, idx1, idx2, K)

if nargin < 5, K = 10; end
if nargin < 4, idx2 = 2; end % target cell group id (receptor gene expressed)
if nargin < 3, idx1 = 1; end % source cell group id (ligand gene expressed)

[A] = logical(full(sc_knngraph(ste.xy, K)));
MASK = logical(group == idx1*(group == idx2)');
MASK1 = MASK & A;
MASK2 = MASK & (~MASK1);

cx1 = zeros(size(group));
[i1, i2] = find(MASK1);
cx1(i1) = 1; % in neighbor cell type 1   e.g., immune cells
cx1(i2) = 2; % in neighbor cell type 2   e.g., tumor cells

cx2 = zeros(size(group));
[i1, i2] = find(MASK2);
cx2(i1) = 1; % out neighbor cell type 1  e.g., immune cells
cx2(i2) = 2; % out neighbor cell type 2  e.g., tumor cells
% cx={cx1,cx2};

X11 = ste.sce.X(:, cx1 == 1); % immune cells neighbor
X12 = ste.sce.X(:, cx1 == 2); % tumor cells neighbor
X21 = ste.sce.X(:, cx2 == 1); % immune cells nonneighbor
X22 = ste.sce.X(:, cx2 == 2); % tumor cells nonneighbor

n = min([size(X11, 2), size(X21, 2)]);
X11 = X11(:, 1:n); % immune cells neighbor
X21 = X21(:, 1:n); % immune cells nonneighbor

n = min([size(X12, 2), size(X22, 2)]);
X12 = X12(:, 1:n); % tumor cells neighbor
X22 = X22(:, 1:n); % tumor cells nonneighbor


[T1] = ten.sctenifoldnet(X11, X21, ste.sce.g); % immune
[T2] = ten.sctenifoldnet(X12, X22, ste.sce.g); % tumor

end
