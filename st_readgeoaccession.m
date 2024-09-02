function [ste] = st_readgeoaccession(acc)

% scalef=1.0;
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

c1 = c(contains(c, ["jpg", "png", "tif"], 'IgnoreCase', true));
if isempty(c1), error('JPG/PNG/TIF file not found.'); end
f2 = i_setupfile(c1);
if isempty(f2), error('JPG/PNG/TIF file name not processed.'); end
img = imread(f2);
%imshow(img)

posnotfound = false;

% c1 = c(contains(c, 'csv') & contains(c, 'positions'));
% if isempty(c1)
%     c1 = c(contains(c, 'txt') & contains(c, 'positions'));
%     if isempty(c1)
%         c1 = c(contains(c, 'tsv') & contains(c, 'positions'));
%         if isempty(c1)
%             c1 = c(contains(c, 'tsv') & contains(c, 'positions'));
%             if isempty(c1)
%                 posnotfound = true;
%                 % error('TXT/CSV/TSV file not found.');
%             end
%         end
%     end
% end

c1 = c(contains(c, {'csv', 'txt', 'tsv'}) & contains(c, ...
       {'positions', 'coordinates'}));
if isempty(c1), posnotfound = true; end

if ~posnotfound % found position file
    if length(c1) > 1, c1 = c1(1); end
    f1 = i_setupfile(c1);
    if isempty(f1), error('TISSUE_POSITION CSV/TXT/TSV file name not processed.'); end
    f1
    T = readtable(f1);

else
    error('TISSUE_POSITION CSV/TXT/TSV file name not processed.');
    % assert(all(ismember(sce.c_cell_id,string(T.Var1))))
    % [~,idx]=ismember(sce.c_cell_id,string(T.Var1));
    % t=T(idx,:);
    %s=[t.Var3,t.Var4];
end

jsonnotfound = false;
c1 = c(contains(c, 'json'));
if isempty(c1)
    jsonnotfound = true;
end
if ~jsonnotfound
    if length(c1) > 1, c1 = c1(1); end
    f1 = i_setupfile(c1);
    if ~isempty(f1)
        txt = fileread(f1);
        scalef = jsondecode(txt);
        %scalef=value.tissue_hires_scalef;
        %scalef
    end
end

scedone = false;
try
    sce = sc_readgeoaccession(acc);
    scedone = true;
catch
    scedone = false;
end

if ~scedone
    mtxnotfound = false;
    c1 = c(contains(c, 'h5') & contains(c, 'matrix'));
    if isempty(c1)
        mtxnotfound = true;
        % error('TXT/CSV/TSV file not found.');
    end
    if ~mtxnotfound
        if length(c1) > 1, c1 = c1(1); end

        f1 = i_setupfile(c1);
        if isempty(f1), error('MATRIX H5 file name not processed.'); end

        [X, genelist, barcodes] = sc_readhdf5file(f1);
        sce = SingleCellExperiment(X, genelist);

        sce.c_cell_id = barcodes;
        metainfo = sprintf("Source: %s", acc);
        sce = sce.appendmetainfo(metainfo);
    end
end

[~, idx] = ismember(sce.c_cell_id, string(T.Var1));
t = T(idx, :);
xy = [t.Var3, t.Var4];

xy = xy - mean(xy);
r = [range(xy(:, 1)), range(xy(:, 2))];
a = size(img, 1:2);
xy = xy .* (0.65 * a ./ r);
xy = xy + a / 2;

% sce.s=s;
%    ste.NumCells
%    size(s,1)
%    assert(ste.NumCells==size(s,1))
%    ste.s=s;
% end

ste = SpatialTranscriptomicsExperiment(sce, xy, img);
ste.tissue_positions_list = T;
ste.scalefactors_json = scalef;
metainfo = sprintf("Source: %s", acc);
ste = ste.appendmetainfo(metainfo);

end

function [f] = i_setupfile(c)
f = [];
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
    if strcmpi(extractAfter(x, length(x)-3), '.gz')
        files = gunzip(x, tmpd);
        f = files{1};
    else
        f = websave(tempname, x);
    end
catch ME
    disp(ME.message)
    rethrow(ME)
end
end
