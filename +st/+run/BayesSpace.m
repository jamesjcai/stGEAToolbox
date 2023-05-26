function [T,X]=BayesSpace(ste)

T=[];
X=[];
isdebug=false;
rscriptdir='R_BayesSpace';
oldpth=pwd();
[isok,msg]=commoncheck_R(rscriptdir);
if ~isok, error(msg); end

tmpfilelist={'input.txt','input.mat','positions.csv', ...
    'output.h5','positions_enhanced.csv','Rplots.pdf'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

% writematrix(ste.sce.X,'input.csv');
% sc_writefile('input.txt',ste.sce.X, ste.sce.g);

X=ste.sce.X;
save('input.mat','X','-v7.3')
writepositions(ste,'positions.csv');

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if ~exist('output.h5','file'), return; end
web('Rplots.pdf','-browser')
T=readtable("positions_enhanced.csv");
scatter(T.row,T.col,[],T.spatial_cluster)
X=h5read('output.h5','/X');

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

