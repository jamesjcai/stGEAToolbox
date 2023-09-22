function st_showccires2(T, cx, ste)

answer = questdlg('Show cell-cell interaction scores of top 5 gene pairs?');
if ~strcmp(answer, 'Yes'), return; end

M = logical((cx{1} == 1)*(cx{1} == 2)');
A = logical(full(sc_knngraph(ste.xy, 10)));
M = M & A;
[x1, x2] = find(M);
a = ste.xy(x1, :);
b = ste.xy(x2, :);
ab = [mean([a(:, 1), b(:, 1)], 2), mean([a(:, 2), b(:, 2)], 2)];

% I2=sub2ind(size(M),x1,x2);
% I=M;
% isequal(I,I2)

for kx = 1:5
    x = ste.sce.X(ste.sce.g == T.lgene(kx), :)';
    y = ste.sce.X(ste.sce.g == T.rgene(kx), :)';

    xy = x * y';
    xy = xy .* M;


    % figure;
    %a=ste.xy(x1,:);
    % scatter(a(:,1),a(:,2));
    %b=ste.xy(x2,:);
    % hold on
    % scatter(b(:,1),b(:,2));
    % axis ij
    %xy=log(xy+1);

    cnew = xy(M);

    %         figure;
    %         scatter(ab(:,1),ab(:,2),40,cnew);
    %         axis ij
    %         colormap jet
    %         box on


    figure
    scatter3(a(:, 1), a(:, 2), zeros(size(a(:, 1))), 10, ...
        [0, 0.4470, 0.7410], 'filled');
    hold on;
    scatter3(b(:, 1), b(:, 2), zeros(size(b(:, 1))), 10, ...
        [0.8500, 0.3250, 0.0980], 'filled');
    stem3(ab(:, 1), ab(:, 2), cnew, 'marker', 'none', 'color', 'g');
    title(sprintf('%s - %s', T.lgene(kx), T.rgene(kx)));
    axis ij
    box on
end

end