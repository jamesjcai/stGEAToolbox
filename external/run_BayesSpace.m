function [T,X]=run_BayesSpace(ste)

rscriptdir='R_BayesSpace';
oldpth=pwd();
[isok,msg]=commoncheck_R(rscriptdir);
if ~isok, error(msg); end

pkg.i_deletefiles({'input.txt','input.mat','positions.csv', ...
    'output.h5','positions_enhanced.csv','Rplots.pdf'});

% writematrix(ste.sce.X,'input.csv');
% sc_writefile('input.txt',ste.sce.X, ste.sce.g);

X=ste.sce.X;
save('input.mat','X','-v7.3')
writepositions(ste,'positions.csv');
pkg.RunRcode('script.R');

web('Rplots.pdf','-browser')
T=readtable("positions_enhanced.csv");
scatter(T.row,T.col,[],T.spatial_cluster)
X=h5read('output.h5','/X');

cd(oldpth);
end

