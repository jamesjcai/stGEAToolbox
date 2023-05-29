function [T,X]=BayesSpace(ste,markers)

if nargin<2, markers=[]; end

% assert(all(ismember(markers,ste.sce.g)))
T=[]; X=[];
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

figure;
scatter(T.row,T.col,[],T.spatial_cluster)
X=h5read('output.h5','/X');
figure;
subplot(2,2,1)
scatter(T.row,T.col,[],X(1,:))
subplot(2,2,2)
scatter(T.row,T.col,[],X(2,:))
subplot(2,2,3)
scatter(T.row,T.col,[],X(3,:))
subplot(2,2,4)
scatter(T.row,T.col,[],X(4,:))


if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

