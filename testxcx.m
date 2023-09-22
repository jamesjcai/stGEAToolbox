load example_data\DuctualCarcinoma.mat
c = ste.sce.struct_cell_clusterings.sc3;
[T, RES, cx] = st_ridgecci(ste, c, 3, 4, 10);

%%
% [T1_tumor,T2_immu]=st_ridgenet(ste,c,3,4,10);

%%


[A] = logical(full(sc_knngraph(ste.xy, 50)));

x1 = 4;
x2 = 3;

MASK0 = logical(c == x1*(c == 3)'); % 1-2
MASK1 = MASK0 & A; % neib = 1-2 & knnlink
MASK2 = MASK0 & (~MASK1); % nonneib = 1-2 & not neib


MASKZ = logical(c == x1*(c ~= x1 & c ~= x2)') | logical(c == x2*(c ~= x1 & c ~= x2)');
% MASK=logical(c==x1*(c==2)');   % 1-2
% MASK3=MASK&A;           % neib but not 1-2
MASK4 = MASK2 & ~MASKZ; % nonneib = 1-2 & not neib

%%
[i1, i2] = find(MASK4);
newc = zeros(size(c));
newc(i1) = 1;
newc(i2) = 2;
figure;
scatter(ste.xy(:, 1), ste.xy(:, 2), 15, newc, 'filled');
axis ij;
colormap lines(3)

return;

%%

[i1, i2] = find(MASK1);
newc = zeros(size(c));
newc(i1) = 1;
newc(i2) = 2;

figure;
scatter3(ste.xy(:, 1), ste.xy(:, 2), ones(size(ste.xy(:, 2))), 5, newc);
axis ij;
colormap lines(3)
hold on

[i1, i2] = find(MASK3);
newc = zeros(size(c));
newc(i1) = 1;
newc(i2) = 2;
scatter3(ste.xy(:, 1), ste.xy(:, 2), zeros(size(ste.xy(:, 2))), 5, newc, 'filled');
axis ij;
colormap lines(3)
