function [T,sortid]=st_anova(ste,group)
%ST_ANOVA 
    % 
    % group : 
    %
    % SEE ALSO: 

if nargin<2
    [group]=st_cluster_img(ste,'method','resnet18kmeans','k',5);
end
X=sc_norm(ste.sce.X);
pval=zeros(ste.sce.NumGenes,1);
F=zeros(ste.sce.NumGenes,1);
for k=1:ste.sce.NumGenes
    [pval(k),t]=anova1(X(k,:),group,'off');
    F(k)=t{2,5};
end
pval(isnan(pval))=1.0;
F(isnan(F))=0.0;
gene=ste.sce.g;
[~,~,~,fdr]=pkg.fdr_bh(pval);
T=table(gene,F,pval,fdr);
[T,sortid]=sortrows(T,'F','descend');

