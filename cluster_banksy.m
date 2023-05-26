function [idx]=cluster_banksy(X,xy,numclass,nnk,lamda)

%REF: https://doi.org/10.1101/2022.04.14.488259

if nargin<3, numclass=4; end
if nargin<4, nnk=10; end
if nargin<5, lamda=0.3; end

[mIdx] = knnsearch(xy,xy,'K',nnk+1);
Y=X;
for l=1:size(mIdx,1)
    Y(:,l)=mean(X(:,mIdx(l,2:end)),2);
end
Data=[lamda*X; (1-lamda)*Y];
idx=run.SC3(Data,numclass);

