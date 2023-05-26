T=readtable('clusterid.txt');
assert(all(ismember(ste.sce.c_cell_id,string(T.ID))))
[~,idx]=ismember(string(T.ID),ste.sce.c_cell_id);
annot_type=string(T.annot_type(idx));
fine_annot_type=string(T.fine_annot_type(idx));
ste.sce.struct_cell_clusterings.annot_type=annot_type;
