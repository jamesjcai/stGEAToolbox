function [idx] = st_cluster_spo(ste, varargin)

p = inputParser;
defaultType = 'sc3';
validTypes = {'sc3'};
checkType = @(x) any(validatestring(x, validTypes));
checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);
addRequired(p, 'ste', @(x) isa(x, 'SpatialTranscriptomicsExperiment'));
addOptional(p, 'method', defaultType, checkType);
% addOptional(p,'plotit',false,@islogical);
% addOptional(p,'sposz',40,checkK);
addOptional(p, 'k', 5, checkK);
parse(p, ste, varargin{:})

method = p.Results.method;
% plotit=p.Results.plotit;
% sposz=p.Results.sposz;
numClass = p.Results.k;

switch method
    case 'sc3'
        sce = ste.sce;
        sce = sce.clustercells(numclass, 'sc3', true);
        [idx] = grp2idx(sce.c_cluster_id);
    otherwise
        error('Invalid option.');
end

end
