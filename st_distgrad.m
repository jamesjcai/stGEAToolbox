function [T,sortid]=st_distgrad(ste,group,idx)
if nargin<3, idx=1; end
if nargin<2, [group]=st_cluster_img(ste,'method','resnet18kmeans','k',5); end
I=group==idx;
xy=ste.xy;
% [D]=pdist2(xy,xy);
[D]=i_knndist(xy);
dt=zeros(size(xy,1),1);
for k=1:size(xy,1)
    d=D(k,:);
    dt(k)=min(d(I));
end
gene=ste.sce.g;
X=sc_norm(ste.sce.X);
pval=zeros(ste.sce.NumGenes,1);
r=zeros(ste.sce.NumGenes,1);
for k=1:ste.sce.NumGenes    
    [r(k),pval(k)]=corr(dt,X(k,:)','type','Spearman');
end
pval(isnan(pval))=1.0;
r(isnan(r))=0.0;
[~,~,~,fdr]=pkg.fdr_bh(pval);
T=table(gene,r,pval,fdr);
[T,sortid]=sortrows(T,'fdr','ascend');


