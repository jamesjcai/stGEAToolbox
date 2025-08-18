function [cx, cy, A] = polycentroid(x, y)
% POLYCENTROID Compute the centroid of a simple polygon.
% [cx, cy, A] = polycentroid(X, Y) returns the centroid (cx, cy) and the
% signed area A of the polygon with vertices (X(i), Y(i)), i=1..N.
% Vertices should be ordered along the boundary (clockwise or counter-clockwise).
%
% [cx, cy, A] = polycentroid(P) accepts an N-by-2 array P = [X Y].
%
% Notes:
% - If the polygon is closed (first vertex equals last), the function
% removes the duplicated last point internally.
% - A is signed: counter-clockwise polygons yield A > 0, clockwise A < 0.
% - If the polygon is degenerate (|A| ~ 0), the centroid is the mean of vertices.

%{

x = [0 2 2 0];
y = [0 0 1 1];
[cx, cy, A] = st.pkg.polycentroid(x, y)
P = [0 0; 2 0; 2 1; 0 1];
[cx, cy, A] = st.pkg.polycentroid(P)
%}

% Handle input formats
if nargin == 1
P = x;
if size(P, 2) ~= 2
error('Single input must be an N-by-2 array of [x y] coordinates.');
end
x = P(:,1);
y = P(:,2);
elseif nargin ~= 2
error('Usage: [cx, cy, A] = polycentroid(x, y) or polycentroid([x y])');
end

x = x(:);
y = y(:);
n = numel(x);

if n < 3
error('Polygon must have at least 3 vertices.');
end

% Remove duplicated last point if polygon is already closed
if x(1) == x(end) && y(1) == y(end)
x(end) = [];
y(end) = [];
n = n - 1;
end

% Close the polygon by appending the first vertex at the end
x2 = [x; x(1)];
y2 = [y; y(1)];

% Shoelace terms
cross = x2(1:end-1).*y2(2:end) - x2(2:end).*y2(1:end-1);
A = 0.5 * sum(cross); % signed area

if abs(A) < eps(max(abs([x; y])))*n
% Degenerate polygon: fall back to mean of vertices
cx = mean(x);
cy = mean(y);
A = 0; % area effectively zero
return;
end

cx = (1/(6*A)) * sum((x2(1:end-1) + x2(2:end)) .* cross);
cy = (1/(6*A)) * sum((y2(1:end-1) + y2(2:end)) .* cross);
end