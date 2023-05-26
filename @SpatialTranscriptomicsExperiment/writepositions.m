function obj = writepositions(obj,filename)
if nargin<2 || isempty(filename)
    filter = {'*.csv'};
    [filename,filepath] = uiputfile(filter);
    if ischar(filename)
        filename=fullfile(filepath,filename);
    end
end

id1=string(obj.tissue_positions_list.Var1);
id2=obj.sce.c_cell_id;
if ~all(ismember(id2,id1)), error('Unknown barcodes.'); end
[~,idx]=ismember(id2,id1);

obj.tissue_positions_list=obj.tissue_positions_list(idx,:);
%obj.tissue_positions_list.Var3=round(obj.xy(:,1));
%obj.tissue_positions_list.Var4=round(obj.xy(:,2));

writetable(obj.tissue_positions_list, ...
    filename,'Delimiter',',','filetype','text', ...
    'WriteVariableNames',false);
end