function [sce, img, xy] = st_readgeoaccessionST(acc)
if nargin < 1, acc = 'GSM4284323'; end
% https://www.biorxiv.org/content/10.1101/2022.01.25.477389v1.full.pdf
% r GSE111672, GSE103322, and GSE144240
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4100721

% if length(strsplit(acc,{',',';',' '}))>1
% end
%acc='GSM4100721';
% acc='GSM4284323'
url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s', acc);
a = webread(url);
b = strsplit(a, '\n');
c = string(b(contains(b, acc)))';
c = c(contains(c, 'ftp'));

if ~(length(c) >= 2)
    disp(url)
    error('Unknown error.');
end

c1 = c(contains(c, ["jpg", "png"], 'IgnoreCase', true));
if isempty(c1), error('JPG/PNG file not found.'); end
f2 = i_setupfile(c1);
if isempty(f2), error('JPG/PNG file name not processed.'); end
img = imread(f2);

txtnotfound = false;
c1 = c(contains(c, 'txt'));
if isempty(c1)
    c1 = c(contains(c, 'csv'));
    if isempty(c1)
        c1 = c(contains(c, 'tsv'));
        if isempty(c1)
            txtnotfound = true;
            % error('TXT/CSV/TSV file not found.');
        end
    end
end
if ~txtnotfound
    if length(c1) > 1
        warning('More than one TXT/CSV/TSV matrix file found. First one is used.');
        c1 = c1(1);
    end
    f1 = i_setupfile(c1);
    if isempty(f1), error('TXT/CSV/TSV file name not processed.'); end
    [X, g, c] = sc_readtsvfile(f1);
end

y1 = all(contains(g, 'x'));
y2 = all(contains(c, 'x'));
if y1 && ~y2
    X = X.';
    xy = i_strsplit2double(g);
    g = c;
elseif ~y1 && y2
    xy = i_strsplit2double(c);
else
    error('xxxx');
end

xy = xy - mean(xy);
r = [range(xy(:, 1)), range(xy(:, 2))];
a = size(img, 1:2);
xy = xy .* (0.65 * a ./ r);
xy = xy + a / 2;


sce = SingleCellExperiment(X, g, xy);
metainfo = sprintf("Source: %s", acc);
sce = sce.appendmetainfo(metainfo);
end

function f = i_setupfile(c)
try
    tmpd = tempdir;
    [x] = regexp(c(1), '<a href="ftp://(.*)">(ftp', 'match');
    x = string(textscan(x, '<a href="ftp://%s'));
    x = append("https://", extractBefore(x, strlength(x)-5));
    if ~(ismcc || isdeployed)
        x = urldecode(x);
    else
        x = pkg.urldecoding(x);
    end
    fprintf('Downloading %s\n', x)
    files = gunzip(x, tmpd);
    f = files{1};
catch ME
    disp(ME.message)
    f = [];
end
end

% function f=i_setupfile2(c)
%     try
%         tmpd=tempname;
%         [x]=regexp(c(1),'<a href="ftp://(.*)">(ftp','match');
%         x=string(textscan(x,'<a href="ftp://%s'));
%         x=append("https://", extractBefore(x,strlength(x)-5));
%         if ~(ismcc || isdeployed)
%             x=urldecode(x);
%         else
%             x=pkg.urldecoding(x);
%         end
%         fprintf('Downloading %s\n',x)
%         f=websave(tmpd,x);
%     catch
%         f=[];
%     end
% end


function b = i_strsplit2double(a)
b = zeros(length(a), 2);
for k = 1:length(a)
    b(k, :) = str2double(strsplit(a(k), 'x'));
end
end
