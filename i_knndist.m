function [D]=i_knndist(s,k)
if nargin<2, k=4; end
[A]=sc_knngraph(s,k);
G=graph(A);
D = distances(G);
% D=zeros(size(s,1));
% for k=1:size(s,1)-1
%     for l=k+1:size(s,1)
% 
%     end
% end

