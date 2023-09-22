function st_scatter_ste(ste, varargin)

p = inputParser;
checkType = @(x) isa(x, 'SpatialTranscriptomicsExperiment');
addRequired(p, 'ste', checkType);
addOptional(p, 'c', sum(ste.sce.X))
addOptional(p, 'dsize', 15)

parse(p, ste, varargin{:})
c = p.Results.c;
dsize = p.Results.dsize;

import mlreportgen.ppt.*;
if ismcc || isdeployed
    makePPTCompilable();
end
mfolder = fileparts(mfilename('fullpath'));
sce = ste.sce;
xy = ste.xy;
img = ste.img;

%c=sum(sce.X);
%c=log(c(:)+1);

FigureHandle = figure('Name', 'STGEATool - Spatial Transcriptomic Gene Expression Analysis Tool', ...
    'position', round(1.25*[0, 0, 560, 420]), ...
    'visible', 'off');
movegui(FigureHandle, 'center');
set(findall(FigureHandle, 'ToolTipString', 'Link/Unlink Plot'), 'Visible', 'Off')
set(findall(FigureHandle, 'ToolTipString', 'Edit Plot'), 'Visible', 'Off')
set(findall(FigureHandle, 'ToolTipString', 'Open Property Inspector'), 'Visible', 'Off')

a = findall(FigureHandle, 'tag', 'figMenuWindow');
delete(a);
a = findall(FigureHandle, 'tag', 'figMenuDesktop');
delete(a);
a = findall(FigureHandle, 'tag', 'figMenuUpdateFileNew');
delete(a);
a = findall(FigureHandle, 'tag', 'figMenuOpen');
a.MenuSelectedFcn = 'stgeatool';
a = findall(FigureHandle, 'tag', 'figMenuFileSaveAs');
delete(a);
a = findall(FigureHandle, 'tag', 'figMenuFileSave');
a.MenuSelectedFcn = @callback_SAVESTE;
a = findall(FigureHandle, 'tag', 'figMenuGenerateCode');
delete(a);

hAx = axes('Parent', FigureHandle);

m_ext = uimenu(FigureHandle, 'Text', 'E&xternal', 'Accelerator', 'x');
i_addmenu(m_ext, 0, @gui.i_setrenv, 'Check R Environment');
i_addmenu(m_ext, 0, @gui.i_setpyenv, 'Check Python Environment');
i_addmenu(m_ext, 1, @callback_BAYESSPACE, 'Run BayesSpace [PMID:34083791]...');


m_exp = uimenu(FigureHandle, 'Text', 'Ex&perimental', 'Accelerator', 'p');
i_addmenu(m_exp, 0, @callback_SHOWRIDGESPOTS, 'Show ridge cells...');
i_addmenu(m_exp, 0, @callback_SHOWRIDGESPOTS2, 'Show ridge cells 2...');
%i_addmenu(m_exp,1,@callback_TENIFOLDXCT,'scTenifoldXct...');
i_addmenu(m_exp, 1, {@callback_RIDGECCI, 1}, 'Detect cell-cell interactions using RidgeCCI comparison...');
i_addmenu(m_exp, 0, {@callback_RIDGECCI, 2}, 'Detect cell-cell interactions using RidgeCCI permutation...');
i_addmenu(m_exp, 1, @gui.callback_ViewMetaData, 'View Metadata...');
i_addmenu(m_exp, 1, {@gui.i_savemainfig, 3}, 'Save Figure to PowerPoint File...');
i_addmenu(m_exp, 0, {@gui.i_savemainfig, 2}, 'Save Figure as Graphic File...');
i_addmenu(m_exp, 0, {@gui.i_savemainfig, 1}, 'Save Figure as SVG File...');
i_addmenu(m_exp, 1, {@(~, ~) web('https://github.com/jamesjcai/stGEAToolbox')}, 'Visit stGEAToolbox GitHub Site...');
i_addmenu(m_exp, 0, @st.gui.callback_CheckUpdates, 'Check for Updates...');

    currentrotview = 2;
    % dsize=15;


    h2 = plotimg;
    slidshown = true;
    hold on

    h1 = plotspo;
    %dt = datacursormode;
    %dt.UpdateFcn = {@i_myupdatefcnx1};

    %hold on
    box on
    grid on

    dt = datacursormode;
    dt.UpdateFcn = {@i_myupdatefcnx};

    % i_seth1cdata(true,'Library Size');
    % i_seth1cdata(true,'nFeatures, log');

    view(currentrotview);

    DftoolbarHandle = findall(FigureHandle, 'tag', 'FigureToolBar');
    UitoolbarHandle = uitoolbar('Parent', FigureHandle);
    MitoolbarHandle = uitoolbar('Parent', FigureHandle);
    % MitoolbarHandle2 = uitoolbar('Parent', FigureHandle);
    set(UitoolbarHandle, 'Tag', 'FigureToolBar');
    set(MitoolbarHandle, 'Tag', 'AlignToolBar', 'Visible', 'off');
    % set(MitoolbarHandle2, 'Tag', 'AlignToolBar2', 'Visible', 'on');

    i_addbutton(MitoolbarHandle, 'off', {@callback_ZOOM, 1, 'a'}, 'zoomin.gif', 'Enlarge XY');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ZOOM, -1, 'a'}, 'zoomout.gif', 'Reduce XY');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ZOOM, 1, 'x'}, 'zoominxx.gif', 'Enlarge X');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ZOOM, -1, 'x'}, 'zoomoutxx.gif', 'Reduce X ');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ZOOM, 1, 'y'}, 'zoominyy.gif', 'Enlarge Y');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ZOOM, -1, 'y'}, 'zoomoutyy.gif', 'Reduce Y');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ROTATION, 1}, 'rotation2.gif', 'Rotation theta degree');
    i_addbutton(MitoolbarHandle, 'off', {@callback_ROTATION, 0}, 'rotation.gif', 'Rotation 90 degree');
    i_addbutton(MitoolbarHandle, 'off', {@callback_FLIP, 2}, 'flipy.gif', 'Flip along X');
    i_addbutton(MitoolbarHandle, 'off', {@callback_FLIP, 1}, 'flipx.gif', 'Flip along Y');

    ptx = uitoggletool(MitoolbarHandle, 'Separator', 'on');
    ptx.CData = i_get_ptImage('hand.gif');
    ptx.Tooltip = 'Move spot plot';
    ptx.ClickedCallback = @callback_MOVEIT;

    i_addbutton_toggle(1, 0, {@togglebtfun, @callback_TURNMIONOFF, "icon-mat-unfold-more-10.gif", ...
        "icon-mat-unfold-less-10.gif", false}, "Turn on/off user onboarding toolbar");

    %i_addbutton(UitoolbarHandle,'off',@callback_TURNMIONOFF,'reficon.gif','Manual alignment tool on/off');

    i_addbutton(UitoolbarHandle, 'on', @callback_DOTSIZE, 'noun_font_size_591141.gif', 'Increase dot size');
    i_addbutton(UitoolbarHandle, 'off', @callback_DOTALPHA, 'plotpicker-rose.gif', 'Increase dot transparency');
    i_addbutton(UitoolbarHandle, 'off', @callback_RANDCOLOR, "plotpicker-compass.gif", "Pick new color map")
    %i_addbutton(UitoolbarHandle,'off',@gui.callback_PickColorMap,"plotpicker-compass.gif","Pick new color map")
    i_addbutton(UitoolbarHandle, 'off', @callback_2D3D, 'dattut12.gif', '2D/3D');
    i_addbutton(UitoolbarHandle, 'on', @callback_HIDEIMG, 'plotpicker-geobubble.gif', 'Toggle slide image');
    i_addbutton(UitoolbarHandle, 'off', @callback_HIDESPO, 'plotpicker-geobubble4.gif', 'Toggle spots');
    i_addbutton(UitoolbarHandle, 'on', @callback_NOTHING, "IMG00107.GIF", " ");
    %i_addbutton(UitoolbarHandle,'off',@callback_NOTHING,"IMG00107.GIF"," ");
    %i_addbutton(UitoolbarHandle,'off',@callback_NOTHING,"IMG00107.GIF"," ");

    i_addbutton(UitoolbarHandle, 'off', @callback_GENEEXPR, 'list.gif', 'Show Gene Expression');
    i_addbutton(UitoolbarHandle, 'off', @callback_SHOWSTATE, "list2.gif", "Show spot state");
    i_addbutton(UitoolbarHandle, 'off', @callback_PICKGENES, "plotpicker-effects.gif", "Filter genes")

    i_addbutton(UitoolbarHandle, 'on', @callback_CLUSTSNE, "plotpicker-dendrogram.gif", "Clustering using KMEANS with tSNE embeddings of spots")
    i_addbutton(UitoolbarHandle, 'off', @callback_CLUSSPO, "plotpicker-gscatter.gif", "Clustering using SC3 with expression matrix ")
    i_addbutton(UitoolbarHandle, 'off', @callback_CLUSIMG, "clusimg.gif", "Clustering using image segmentation")

    i_addbutton(UitoolbarHandle, 'on', @callback_CTYPECORES, "cellscore.gif", "Calculate Cell Scores from Cell Type Markers")
    i_addbutton(UitoolbarHandle, 'off', @callback_GENELASSO, "plotpicker-kagi.gif", "Marker genes of brushed spots")
    i_addbutton(UitoolbarHandle, 'off', @callback_MARKHEATMAP, "plotpicker-plotmatrix.gif", "Marker gene heatmap")
    % i_addbutton(UitoolbarHandle,'off',{@gui.callback_MarkerGeneHeatmap,sce},"plotpicker-plotmatrix.gif","Marker gene heatmap")
    i_addbutton(UitoolbarHandle, 'on', @callback_CLUSPOP, "plotpicker-geoscatter.gif", "Show cell clusters/groups individually")
    i_addbutton(UitoolbarHandle, 'on', @callback_SAVESTE, "export.gif", "Export & save data");
    i_addbutton(UitoolbarHandle, 'off', @gui.callback_CloseAllOthers, "icon-fa-cut-10.gif", "Close All Other Figures")

    i_addbutton(DftoolbarHandle, 'off', @callback_NOTHING, "IMG00107.GIF", " ");
    i_addbutton(DftoolbarHandle, 'on', @callback_CELLSCORES, "cellscore2.gif", "Calculate Cell Scores from List of Feature Genes")
    i_addbutton(DftoolbarHandle, 'on', @callback_MULTIVIEW, 'plotpicker-arxtimeseries.gif', 'Multiview');
    i_addbutton(DftoolbarHandle, 'on', @callback_SVGANOVA, 'plotpicker-andrewsplot.gif', 'Identify Spatially Variable Genes (SVGs) using ANOVA');
    i_addbutton(DftoolbarHandle, 'off', @callback_DISTGRAD, 'plotpicker-andrewsplot.gif', 'Identify Spatially Variable Genes (SVGs) using Distance Gradient');

    %i_addbutton(DftoolbarHandle,'on',@callback_SHOWRIDGESPOTS,'HDF_object01.gif','Show ridge cells');
    %i_addbutton(DftoolbarHandle,'off',@callback_SHOWRIDGESPOTS2,'HDF_object01.gif','Show ridge cells 2');
    %i_addbutton(DftoolbarHandle,'on',@callback_TENIFOLDXCT,'HDF_object02.gif','scTenifoldXct');
    %i_addbutton(DftoolbarHandle,'off',{@callback_RIDGECCI,1},'HDF_object02.gif','Detect cell-cell interactions using RidgeCCI comparison');
    %i_addbutton(DftoolbarHandle,'off',{@callback_RIDGECCI,2},'HDF_object03.gif','Detect cell-cell interactions using RidgeCCI permutation');

    i_addbutton(DftoolbarHandle, 'on', @callback_SCGEATOOL, 'greenarrowicon.gif', 'SCGEATOOL...');
    i_addbutton(DftoolbarHandle, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save figure to PowerPoint...');


    title(hAx, ste.title);
    set(FigureHandle, 'visible', 'on');
    gui_startpoint = [];
    gui_currenthandle = [];
    guidata(FigureHandle, ste);


    % ------------------------
    % Callback Functions
    % ------------------------


    function i_savefig(~, ~, tag)
        if tag == 1
            filter = {'*.svg'};
            [filename, filepath] = uiputfile(filter);
            if ischar(filename)
                saveas(FigureHandle, [filepath, filename], 'svg');
            end
        elseif tag == 2
            % axx=gca;
            filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
            [filename, filepath] = uiputfile(filter);
            if ischar(filename)
                exportgraphics(FigureHandle, [filepath, filename]);
            end
        elseif tag == 3
            gui.i_export2pptx({FigureHandle}, {'STGEATOOL'});
        end
    end

    function callback_NOTHING(~, ~)
    end

    function callback_TURNMIONOFF(~, ~)
        %         for k=1:length(MitoolbarHandle.Children)
        %             set(MitoolbarHandle.Children(k),'Enable','off');
        %         end
        switch get(MitoolbarHandle, 'Visible')
            case 'on'
                set(MitoolbarHandle, 'Visible', 'off');
                %set(MitoolbarHandle2,'Visible','on');
            case 'off'
                set(MitoolbarHandle, 'Visible', 'on');
                %set(MitoolbarHandle2,'Visible','off');
        end
    end

    function callback_CLUSPOP(~, ~)
        [thisc, ttxt] = gui.i_select1clusterings(sce);
        if isempty(thisc), return; end
        [c] = grp2idx(thisc);
        i_seth1cdata(false, ttxt);
        cmv = 1:max(c);
        idxx = cmv;
        [cmx] = countmember(cmv, c);
        answer = questdlg('Sort by size of cell groups?');
        if strcmpi(answer, 'Yes')
            [~, idxx] = sort(cmx, 'descend');
        end
        totaln = max(c);
        numfig = ceil(totaln/9);
        for nf = 1:numfig
            f = figure('visible', 'off');
            for k = 1:9
                kk = (nf - 1) * 9 + k;
                if kk <= totaln
                    subplot(3, 3, k);
                    gui.i_gscatter3(xy, c, 3, cmv(idxx(kk)));
                    title(sprintf('Cluster #%d\n%d spots (%.2f%%)', ...
                        idxx(kk), cmx(idxx(kk)), ...
                        100*cmx(idxx(kk))/length(c)));
                end
                % colormap(para.oldColorMap);
                % colormap(lines(length(unique(c))));
                colormap(colormap(hAx));
            end
            P = get(f, 'Position');
            set(f, 'Position', [P(1) - 20 * nf, P(2) - 20 * nf, P(3), P(4)]);
            set(f, 'visible', 'on');
            drawnow;
        end
    end

    function callback_CTYPECORES(src, ~)
        [cs, ttxt] = gui.callback_CellTypeMarkerScores(src, [], sce);
        if ~isempty(cs) && ~isempty(ttxt)
            c = cs;
            i_seth1cdata(true, ttxt);
        end
    end

    function callback_GENELASSO(~, ~)
        ptsSelected = logical(h1.BrushData.');
        if ~any(ptsSelected)
            helpdlg("No spots are selected.");
            return;
        end
        [numfig] = gui.i_inputnumg;
        if isempty(numfig), return; end
        fw = gui.gui_waitbar;
        y = double(ptsSelected);
        sce.c = 1 + ptsSelected;
        X = sce.X';
        try
            if issparse(X), X = full(X); end
            [B] = lasso(X, y, 'DFmax', numfig*3, 'MaxIter', 1e3);
        catch ME
            gui.gui_waitbar(fw, true);
            errordlg(ME.message);
            return;
        end
        [~, ix] = min(abs(sum(B > 0)-numfig));
        b = B(:, ix);
        idx = b > 0;
        gui.gui_waitbar(fw);
        if ~any(idx)
            errordlg('No marker gene found');
            return;
        else
            glist = sce.g(idx);
            [~, jx] = sort(b(idx), 'descend');
            glist = glist(jx);
            for k = 1:length(glist)
                i_cascadeexpr(sce, glist(k), xy, k);
            end
        end
    end

    function callback_CELLSCORES(~, ~)

        %{
        actiontype=questdlg('Select a predefined score or define a new score?',...
            'Select/Define Score','Select Predefined Score',...
            'Define New Score','Select Predefined Score');
        switch actiontype
            case 'Select Predefined Score'
                [~,T]=pkg.e_cellscores(sce.X,sce.g,0);

                listitems=T.ScoreType;
                [indx,tf] = listdlg('PromptString',...
                    {'Select Class','',''},...
                    'SelectionMode','single','ListString',listitems,'ListSize',[220,300]);
                if ~tf==1, return; end
                %fw=gui.gui_waitbar;
                [cs]=pkg.e_cellscores(sce.X,sce.g,indx);
                ttxt=T.ScoreType(indx);
                %gui.gui_waitbar(fw);
            case 'Define New Score'
                ttxt='Customized Score';
                % Pd1=pdcd1 tim3=HAVCR2, tcf1=HNF1A  https://www.nature.com/articles/s41577-019-0221-9
                % posgcandidates=["PDCD1","HNF1A","HAVCR2","KLRG1","CD44","LY6C","CTLA","ICOS","LAG3"];
                %posgcandidates=sce.g(randi(length(sce.g),10,1));
                [posg]=gui.i_selectngenes(sce.g);
                if isempty(posg)
                    helpdlg('No feature genes selected.','')
                    return;
                end
                fw=gui.gui_waitbar;
                a=sprintf('%s+',posg);
                a=a(1:min([length(a),50]));
                ttxt=sprintf('%s\n%s',ttxt,a);
                cs=sc_cellscore(sce.X,sce.g,posg);
                gui.gui_waitbar(fw);
            otherwise
                return;
        end
        %}

        try
            [cs, ttxt] = gui.i_computecellscore(sce);
        catch ME
            errordlg(ME.message);
            return;
        end
        if ~isempty(cs) && ~isempty(ttxt)
            c = cs;
            i_seth1cdata(true, ttxt);
        end
    end

    function callback_CLUSIMG(~, ~)
        answer = questdlg('Cluster spots using k-means algorithm with image feature?');
        if ~strcmp(answer, 'Yes'), return; end
        %answer1=questdlg('Select method:','','imsegkmeans', ...
        %    'resnet18kmeans','resnet18labeling','resnet18kmeans');
        %if ~ismember(answer1,{'imsegkmeans','resnet18kmeans','resnet18labeling'}), return; end
        cxdone = false;
        if any(contains(fieldnames(sce.struct_cell_clusterings), 'resnet18kmeans'))
            if ~isempty(sce.struct_cell_clusterings.resnet18kmeans)
                usingexist = questdlg('Using existing clustering result?');
                if isempty(usingexist) || strcmp(usingexist, 'Cancel'), return; end
                if strcmp(usingexist, 'Yes')
                    c = sce.struct_cell_clusterings.resnet18kmeans;
                    cxdone = true;
                end
            end
        end
        if ~cxdone
            numclass = gui.i_inputnumk;
            if isempty(numclass), return; end
            plotit = questdlg('Show intermediate image output?', '', 'Yes', 'No', 'No');
            if strcmp(plotit, 'Yes')
                plotit = true;
            elseif strcmp(plotit, 'No')
                plotit = false;
            else
                return;
            end
            fw = gui.gui_waitbar;
            try
                [idx] = cluster_resnet18(img, xy, numclass, 40, plotit);
                %[idx]=st_cluster_img(ste,'method','resnet18kmeans', ...
                %        'plotit',answer2,'k',numclass);
            catch ME
                errordlg(ME.message);
                gui.gui_waitbar(fw, true);
                return;
            end
            gui.gui_waitbar(fw);
            [c] = grp2idx(idx);
            sce.struct_cell_clusterings(:).('resnet18kmeans') = c;
        end
        i_seth1cdata(false, 'resnet18kmeans cluster id');
    end

    function callback_CLUSTSNE(~, ~)
        answer = questdlg('Cluster spots using k-means with tSNE embedding?');
        if ~strcmp(answer, 'Yes'), return; end
        cxdone = false;
        if any(contains(fieldnames(sce.struct_cell_clusterings), 'kmeans'))
            if ~isempty(sce.struct_cell_clusterings.('kmeans'))
                usingexist = questdlg('Using existing clustering result?');
                if isempty(usingexist) || strcmp(usingexist, 'Cancel'), return; end
                if strcmp(usingexist, 'Yes')
                    c = sce.struct_cell_clusterings.('kmeans');
                    cxdone = true;
                end
            end
        end
        if ~cxdone
            k = gui.i_inputnumk;
            if isempty(k), return; end
            fw = gui.gui_waitbar;
            try
                sce = sce.embedcells('tsne', true, true, 3);
                sce = sce.clustercells(k, 'kmeans', true);
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            gui.gui_waitbar(fw);
            [c] = grp2idx(sce.c_cluster_id);
        end
        % unique(sce.c_cluster_id)
        i_seth1cdata(false, 'tSNE cluter id');
    end

    function callback_MARKHEATMAP(src, events)
        gui.callback_MarkerGeneHeatmap(src, events, sce);
    end

    function callback_CLUSSPO(~, ~)
        answer = questdlg('Cluster spots using SC3 algorithm with expression matrix?');
        if ~strcmp(answer, 'Yes'), return; end
        cxdone = false;
        if any(contains(fieldnames(sce.struct_cell_clusterings), 'sc3'))
            if ~isempty(sce.struct_cell_clusterings.sc3)
                usingexist = questdlg('Using existing clustering result?');
                if isempty(usingexist) || strcmp(usingexist, 'Cancel'), return; end
                if strcmp(usingexist, 'Yes')
                    c = sce.struct_cell_clusterings.sc3;
                    cxdone = true;
                end
            end
        end
        if ~cxdone
            numclass = gui.i_inputnumk;
            if isempty(numclass), return; end
            fw = gui.gui_waitbar;
            try
                sce = sce.clustercells(numclass, 'sc3', true);
            catch ME
                errordlg(ME.message);
                gui.gui_waitbar(fw, true);
                return;
            end
            gui.gui_waitbar(fw);
            [c] = grp2idx(sce.c_cluster_id);
        end
        i_seth1cdata(false, 'sc3 cluter id');
    end

    function Brushed2NewCluster(~, ~)
        answer = questdlg('Make a new cluster out of brushed cells?');
        if ~strcmp(answer, 'Yes'), return; end
        ptsSelected = logical(h1.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return;
        end
        sce.c_cluster_id(ptsSelected) = max(sce.c_cluster_id) + 1;
    end

    function callback_SCGEATOOL(~, ~)
        scgeatool(sce)

        %assignin('base','sce',sce);
        %eval('scgeatool(sce)');

        % f=sc_scatter_sce(sce, ...
        %     'callinghandle',FigureHandle);
        % f.Position=[0.9*f.Position([1 2]) 0.6*f.Position([3 4])];
        % waitfor(f);
        % if ~isempty(guidata(FigureHandle))
        % sce = guidata(FigureHandle);
        % title(hAx,sce.title);
        % end
    end

    function callback_SVGANOVA(~, ~)
        [thisc, ttxt] = gui.i_select1clusterings(sce);
        if isempty(thisc), return; end
        c = grp2idx(thisc);
        i_seth1cdata(false, ttxt);
        fw = gui.gui_waitbar;
        [T] = st_anova(ste, thisc);
        gui.gui_waitbar(fw);
        waitfor(gui.i_exporttable(T, true, 'T'));
        answer = questdlg('Show top 10 genes'' expression?');
        if strcmp(answer, 'Yes')
            for k = 1:10
                i_cascadeexpr(sce, T.gene(k), xy, k);
            end
        end
    end

    function callback_DISTGRAD(~, ~)
        [thisc, ttxt] = gui.i_select1clusterings(sce);
        if isempty(thisc), return; end
        c = grp2idx(thisc);
        i_seth1cdata(false, ttxt);
        idx = listdlg('SelectionMode', 'single', ...
            'ListString', string(1:max(thisc)), ...
            'PromptString', {'Select a target cluster.'});
        if isempty(idx), return; end
        fw = gui.gui_waitbar;
        try
            [T] = st_distgrad(ste, thisc, idx);
        catch ME
            disp(ME.message);
            gui.gui_waitbar(fw, true);
            return;
        end
        gui.gui_waitbar(fw);
        waitfor(gui.i_exporttable(T, true, 'T'));
        answer = questdlg('Show top 10 genes'' expression?');
        if strcmp(answer, 'Yes')
            for k = 1:10
                i_cascadeexpr(sce, T.gene(k), xy, k);
            end
        end
    end

    function [x1, x2] = i_pickridgegrps
        x1 = [];
        x2 = [];
        [thisc, ttxt] = gui.i_select1clusterings(sce);
        if isempty(thisc), return; end

        [c] = grp2idx(thisc);
        assert(isequal(c, thisc));

        i_seth1cdata(false, ttxt);
        [i1, i2] = i_select2grps(c);
        if i1 == 0 || i2 == 0, return; end
        a1 = sprintf('Group %d -> Group %d', i1, i2);
        a2 = sprintf('Group %d -> Group %d', i2, i1);
        answer = questdlg('Select direction (ligand->receptor)', '', a1, a2, a1);
        switch answer
            case a1
                x1 = i1;
                x2 = i2;
            case a2
                x1 = i2;
                x2 = i1;
            otherwise
                return;
        end
    end

    function callback_SHOWRIDGESPOTS(~, ~)
        [x1, x2] = i_pickridgegrps;
        if isempty(x1) || isempty(x2), return; end

        K = gui.i_inputnumk(10);
        if isempty(K), return; end

        [A] = sc_knngraph(ste.xy, K);
        newc1 = zeros(ste.sce.NumCells, 1);
        newc2 = zeros(ste.sce.NumCells, 1);
        fw = gui.gui_waitbar;

        nypairs = 0;
        nnpairs = 0;

        for ix = 1:ste.sce.NumCells
            for jx = 1:ste.sce.NumCells
                if c(ix) == x1 && c(jx) == x2
                    if A(ix, jx) > 0
                        nypairs = nypairs + 1;
                        newc1(ix) = 1;
                        newc1(jx) = 2;
                    else
                        nnpairs = nnpairs + 1;
                        newc2(ix) = 1;
                        newc2(jx) = 2;
                    end
                end
            end
        end
        gui.gui_waitbar(fw);

        %c=newc1;
        %i_seth1cdata(false,'Ridge Spots (1-source, 2-target)');
        f = figure('Visible', 'off');
        subplot(1, 2, 1);
        scatter(ste.xy(:, 1), ste.xy(:, 2), 15, newc1, 'filled');
        axis ij;
        colormap(lines(length(unique(newc1))));
        box on
        subplot(1, 2, 2);
        scatter(ste.xy(:, 1), ste.xy(:, 2), 15, newc2, 'filled');
        axis ij;
        colormap(lines(length(unique(newc2))));
        box on
        f.Position(3) = f.Position(3) * 2;
        movegui(f, 'center');
        set(f, 'visible', 'on');
        helpdlg(sprintf('There are %d ridge spot pairs and %d non-ridge spot pairs. To obtain more ridge spot pairs, increase K of the K-NN.', ...
            nypairs, nnpairs), sprintf('K=%d', K));
    end

    function callback_SHOWRIDGESPOTS2(~, ~)
        [x1, x2] = i_pickridgegrps;
        if isempty(x1) || isempty(x2), return; end

        K = gui.i_inputnumk(10);
        if isempty(K), return; end

        [A] = logical(full(sc_knngraph(ste.xy, K)));
        newc = zeros(ste.sce.NumCells, 1);

        nypairs = 0;
        nnpairs = 0;
        MASK0 = logical(c == x1*(c == x2)'); % 1-2
        MASK1 = MASK0 & A; % neib = 1-2 & knnlink
        MASK2 = MASK0 & (~MASK1); % nonneib = 1-2 & not neib


        MASK = logical(c == x1*(c ~= x1 & c ~= x2)'); % | logical(c==x2*(c~=x1&c~=x2)');
        MASK3 = MASK & A; % neib but not 1-2
        MASK4 = MASK0 & (~MASK1) & (~MASK3); % nonneib = 1-2 & not neib

        isequal(MASK1, MASK3)

        % any(MASK2)
        % cx1=zeros(size(c));
        [i1, i2] = find(MASK1 & MASK3);
        newc(i1) = 1;
        newc(i2) = 2;
        c = newc; % xxx
        i_seth1cdata(false, 'Ridge Spots (1-source, 2-target)');
        %         helpdlg(sprintf('There are %d ridge spot pairs and %d non-ridge spot pairs. To obtain more ridge spot pairs, increase K of the K-NN.',...
        %                 nypairs, nnpairs),sprintf('K=%d',K));
    end

    % function callback_TENIFOLDXCT(~,~)
    %     [x1,x2]=i_pickridgegrps;
    %     if isempty(x1) || isempty(x2), return; end
    %     fw=gui.gui_waitbar;
    %     sce.c_cell_type_tx=string(c);
    %     [T]=run.py_scTenifoldXct(sce,string(x1),string(x2),false);
    %         gui.gui_waitbar(fw);
    %     if ~isempty(T)
    %         gui.i_exporttable(T);
    %     else
    %         helpdlg('No ligand-receptor pairs are identified.','');
    %     end
    % end

    function callback_RIDGECCI(~, ~, methodid)
        if nargin < 3, methodid = 1; end

        [x1, x2] = i_pickridgegrps;
        if isempty(x1) || isempty(x2), return; end

        K = gui.i_inputnumk(5);
        if isempty(K), return; end

        try
            [T, cx] = st_ridgecci(ste, c, x1, x2, K, methodid);
        catch ME
            errordlg(ME.message);
            return;
        end
        waitfor(gui.i_exporttable(T, true, 'T'));
        % i_showccires(T,cx);
        st_showccires2(T, cx, ste);
    end

    %     function callback_RIDGECCIPERM(~,~)
    %         [thisc,ttxt]=gui.i_select1clusterings(sce);
    %         if isempty(thisc), return; end
    %
    %         [c]=grp2idx(thisc);
    %         assert(isequal(c,thisc));
    %
    %         i_seth1cdata(false,ttxt);
    %         [i1,i2]=i_select2grps(c);
    %         if i1==0 || i2==0, return; end
    %         a1=sprintf('Group %d -> Group %d',i1,i2);
    %         a2=sprintf('Group %d -> Group %d',i2,i1);
    %         answer=questdlg('Select direction (ligand->receptor)','',a1,a2,a1);
    %         switch answer
    %             case a1
    %                 x1=i1; x2=i2;
    %             case a2
    %                 x1=i2; x2=i1;
    %             otherwise
    %                 return;
    %         end
    %         %fw=gui.gui_waitbar;
    %         try
    %             [T,RES]=st_ridgecci_perm(ste,c,x1,x2);
    %         catch ME
    %             %gui.gui_waitbar(fw,true);
    %             errordlg(ME.message);
    %             return;
    %         end
    %         %gui.gui_waitbar(fw);
    %         waitfor(gui.i_exporttable(T,true,'T'));
    %         answer=questdlg('Show top 5 gene pairs'' expression?');
    %         if strcmp(answer,'Yes')
    %             for k=1:2:10
    %                 j=(k+1)/2;
    %                 i_cascadeexpr(sce,T.lgene(j),xy,k);
    %                 title(sprintf('Ligand %d: %s',j,T.lgene(j)));
    %                 i_cascadeexpr(sce,T.rgene(j),xy,k+1);
    %                 title(sprintf('Receptor %d: %s',j,T.rgene(j)));
    %                 res=RES{j};
    %                 figure;
    %                 gui.i_boxplot([res{1};res{2}],[zeros(size(res{1})); ones(size(res{2}))]);
    %             end
    %         end
    %     end

    %     function i_showccires(T,cx)
    %         answer=questdlg('Show top 5 gene pairs'' expression?');
    %         if strcmp(answer,'Yes')
    %             OUTppt=[tempname,'.pptx'];
    %             ppt = Presentation(OUTppt);
    %             open(ppt);
    %
    %             images={};
    %
    %             img0=[tempname,'.png'];
    %             saveas(FigureHandle,img0);
    %             slide4 = add(ppt,'Title and Content');
    %             replace(slide4,'Content',Picture(img0));
    %             images = [images {img0}];
    %
    %             c=cx{1};
    %             i_seth1cdata(false,'Ridge Spots (1-source, 2-target)');
    %             img1=[tempname,'.png'];
    %             saveas(FigureHandle,img1);
    %  %           slide4 = add(ppt,'Title and Content');
    %  %           replace(slide4,'Content',Picture(img1));
    %             images = [images {img1}];
    %
    %             c=cx{2};
    %             i_seth1cdata(false,'Ridge Spots (1-source, 2-target)');
    %             img2=[tempname,'.png'];
    %             saveas(FigureHandle,img2);
    % %            slide4 = add(ppt,'Title and Content');
    % %            replace(slide4,'Content',Picture(img1));
    %             images = [images {img2}];
    %
    %             slide4 = add(ppt,'Comparison');
    %             replace(slide4,'Left Content',Picture(img1));
    %             replace(slide4,'Right Content',Picture(img2));
    %
    % %             f=figure('Visible','off');
    % %             scatter(xy(:,1),xy(:,2),10,cx,'filled');
    % %             axis ij;
    % %             colormap(lines(length(unique(cx))));
    % %             img1=[tempname,'.png'];
    % %             saveas(f,img1);
    % %             slide4 = add(ppt,'Title and Content');
    % %             replace(slide4,'Content',Picture(img1));
    % %             images = [images {img1}];
    %
    %
    %             for k=1:2:10
    %                 j=(k+1)/2;
    %                 img1=[tempname,'.png'];
    %                 img2=[tempname,'.png'];
    %                 images = [images {img1}];
    %                 images = [images {img2}];
    %                 f1=i_cascadeexpr(sce,T.lgene(j),xy,k);
    %                 title(sprintf('Ligand %d: %s',j,T.lgene(j)));
    %                 saveas(f1,img1);
    %                 f2=i_cascadeexpr(sce,T.rgene(j),xy,k+1);
    %                 title(sprintf('Receptor %d: %s',j,T.rgene(j)));
    %                 saveas(f2,img2);
    %                 slide4 = add(ppt,'Comparison');
    %                 replace(slide4,'Title',sprintf('Ligand %d - Receptor %d',j,j));
    %                 replace(slide4,'Left Text',sprintf('Ligand %d: %s',j,T.lgene(j)));
    %                 replace(slide4,'Right Text',sprintf('Receptor %d: %s',j,T.rgene(j)));
    %                 replace(slide4,'Left Content',Picture(img1));
    %                 replace(slide4,'Right Content',Picture(img2));
    %             end
    %             fw=gui.gui_waitbar;
    %             close(ppt);
    %             len = length(images);
    %             for i = 1:len
    %                 delete(images{i});
    %             end
    %             gui.gui_waitbar(fw);
    %             rptview(ppt);
    %         end
    %     end

    function callback_SHOWSTATE(~, ~)
        [thisc, ttxt] = gui.i_select1clusterings(sce);

        if isempty(thisc), return; end
        [c, cL] = grp2idx(thisc);
        dtp = findobj(h1, 'Type', 'datatip');
        delete(dtp);
        row = dataTipTextRow('', cL(c));
        h1.DataTipTemplate.DataTipRows = row;

        %             for i = 1:max(c)
        %                 idx = find(c == i);
        %                 siv = xy(idx, :);
        %                 si = mean(siv, 1);
        %                 [k] = dsearchn(siv, si);
        %                 datatip(h1, 'DataIndex', idx(k));
        %             end


        i_seth1cdata(false, ttxt);
    end

    function callback_PICKGENES(~, ~)
        answer = inputdlg('Expressed in less than % of cells (e.g., 0.05=5%) or # of cells (e.g., 10).', ...
            'Remove Genes', [1, 55], {'0.05'});
        if isempty(answer), return; end
        if iscell(answer)
            a = str2double(answer{1});
            if a > 0 && a < intmax
                n0 = sce.NumGenes;
                fw = gui.gui_waitbar;
                sce = sce.selectkeepgenes(1, a);
                gui.gui_waitbar(fw);
                n1 = sce.NumGenes;
                helpdlg(sprintf('%d genes removed.', n0-n1), '');
                title(hAx, sce.title);
            end
        end
    end

    function callback_DOTSIZE(~, ~)
        if dsize > 50
            dsize = 5;
        else
            dsize = dsize + 5;
        end
        h1.SizeData = dsize;
    end

    function callback_DOTALPHA(~, ~)
        a = h1.MarkerEdgeColor;
        h1.MarkerEdgeColor = h1.MarkerFaceColor;
        h1.MarkerFaceColor = a;
        % a=h1.MarkerFaceAlpha;
        % a=0.9*a;
        % if a<0.3, a=1.0; end
        % h1.MarkerFaceAlpha=a;
    end

    function callback_RANDCOLOR(~, ~)
        n = length(unique(c));
        colormap(i_getcolormaps(n));
    end

    function callback_MULTIVIEW(~, ~)
        [thisc, ttxt] = gui.i_select1clusterings(sce);
        if isempty(thisc), return; end
        [thiss, clable] = gui.i_select1embedding(sce);
        if isempty(thiss), return; end
        f0 = figure('Visible', 'off');
        subplot(1, 2, 1);
        hh1 = scatter(xy(:, 1), xy(:, 2), 10, thisc, 'filled');
        box on;
        grid on;
        title(ttxt);
        subplot(1, 2, 2);
        hh2 = gui.i_gscatter3(sce.s, thisc, 1, 1);
        box on;
        title(clable);
        f0.Position(3) = f0.Position(3) * 2;
        movegui(f0, 'center');
        % st_multigroupings(xy,sce);
        hBr = brush(f0);
        hBr.ActionPostCallback = {@onBrushAction, hh1, hh2};
        set(f0, 'Visible', true);
    end

    function onBrushAction(~, event, hh1, hh2)
        if isequal(event.Axes.Children, hh1)
            hh2.BrushData = hh1.BrushData;
        elseif isequal(event.Axes.Children, hh2)
            hh1.BrushData = hh2.BrushData;
        end
    end

    function h = plotspo
        h = scatter3(hAx, xy(:, 1), ...
            xy(:, 2), ones(size(xy(:, 1))), ...
            dsize, c, 'filled');
    end

    function h = plotimg
        methodid = 2;
        switch methodid
            case 1
                xImage = [0, 1; ...
                    0, 1] * size(img, 1);
                yImage = [1, 1; ...
                    0, 0] * size(img, 2);

                xImage = [0, 1; ...
                    0, 1] * size(img, 1);
                yImage = [0, 0; ...
                    1, 1] * size(img, 2);

                zImage = [-1, -1; -1, -1] + 1; % The z data for the image corners
                h = surf(hAx, xImage, yImage, zImage, ... % Plot the surface
                'CData', img, 'FaceColor', 'texturemap');
                %axis tight
                %axis image
            case 2
                h = imagesc(img);
                %axis ij
                axis tight
            case 3
                h = imshow(img);
        end
        slidshown = true;
    end

    function callback_ZOOM(~, ~, inout, dirtag)
        if nargin < 4, dirtag = 'a'; end
        % prompt = {'Enter matrix size:'};
        % dlgtitle = 'Input';
        % dims = [1 35];
        % definput = {'20','hsv'};
        % answer = inputdlg(prompt,dlgtitle,dims,definput)
        f = 0.01;
        a = xy - mean(xy);
        switch dirtag
            case 'a'
                if inout > 0
                    a = (1 + f) * a;
                else
                    a = (1 - f) * a;
                end
            case 'x'
                if inout > 0
                    a(:, 1) = (1 + f) * a(:, 1);
                else
                    a(:, 1) = (1 - f) * a(:, 1);
                end
            case 'y'
                if inout > 0
                    a(:, 2) = (1 + f) * a(:, 2);
                else
                    a(:, 2) = (1 - f) * a(:, 2);
                end
            otherwise
                error('Invalid option DIRTAG.');
        end
        xy = a + mean(xy);
        delete(h1);
        set(ptx, 'State', 'off');
        h1 = plotspo;
    end

    function callback_FLIP(~, ~, tag)
        delete(h1);
        set(ptx, 'State', 'off');
        xy(:, tag) = -(xy(:, tag)) + 2 * mean(xy(:, tag));
        h1 = plotspo;
    end

    function callback_MOVEIT(~, ~)
        if isempty(get(h1, 'ButtonDownFcn'))
            set(h1, 'ButtonDownFcn', @startmovit);
            brush('off');
            rotate3d('off');
            pan('off');
            zoom('off');
            datacursormode('off');
        else
            set(h1, 'ButtonDownFcn', '');
        end
    end

    function startmovit(src, ~)
        set(gcf, 'PointerShapeCData', nan(16, 16));
        set(gcf, 'Pointer', 'custom');
        gui_currenthandle = src;
        thisfig = gcbf();
        set(thisfig, 'WindowButtonMotionFcn', @movit);
        set(thisfig, 'WindowButtonUpFcn', @stopmovit);
        gui_startpoint = get(gca, 'CurrentPoint');
        set(gui_currenthandle, 'UserData', ...
            {get(gui_currenthandle, 'XData'), get(gui_currenthandle, 'YData')});
    end

    function movit(~, ~)
        try
            if isequal(gui_startpoint, []), return; end
        catch
        end
        pos = get(gca, 'CurrentPoint') - gui_startpoint;
        XYData = get(gui_currenthandle, 'UserData');
        set(gui_currenthandle, 'XData', XYData{1}+pos(1, 1));
        set(gui_currenthandle, 'YData', XYData{2}+pos(1, 2));
        drawnow;
        xy = [get(gui_currenthandle, 'XData')', get(gui_currenthandle, 'YData')'];
    end

    function stopmovit(~, ~)
        thisfig = gcbf();
        set(gcf, 'Pointer', 'arrow');
        set(thisfig, 'WindowButtonUpFcn', '');
        set(thisfig, 'WindowButtonMotionFcn', '');
        drawnow;
        xy = [get(h1, 'XData')', get(h1, 'YData')'];
        % disp('xy updated.');
    end

    function callback_2D3D(~, ~)
        currentrotview = currentrotview + ...
            (currentrotview == 2) - (currentrotview == 3);
        view(currentrotview);
    end

    function callback_ROTATION(~, ~, inputtheta)
        if nargin < 3, inputtheta = false; end
        if inputtheta
            theta = [];
            dgr = inputdlg('Theta (0-360):', 'Input Angle', [1, 40], {'-90'});
            if ~isempty(dgr)
                dgr = dgr{1};
                dgr = str2double(dgr);
                if abs(dgr) > 0 && abs(dgr) < 360
                    theta = pi * dgr / 180;
                end
            else
                return;
            end
            if isempty(theta)
                errordlg('Invalid degree value.');
                return;
            end
        else
            theta = -pi / 2;
        end
        a = xy - mean(xy);
        xy = ([cos(theta), -sin(theta); sin(theta), cos(theta)] * a')' + mean(xy);
        delete(h1);
        set(ptx, 'State', 'off');
        h1 = plotspo;
    end

    function callback_HIDESPO(~, ~)
        set(h1, 'visible', ~get(h1, 'visible'));
    end

    function callback_HIDEIMG(~, ~)
        if slidshown
            a = get(h2, 'CData');
            if ~isempty(a)
                aa = a(1, 1, :);
                aa(:, :, 1) = 240;
                aa(:, :, 2) = 240;
                aa(:, :, 3) = 240;
                b = repmat(aa, [size(a, 1), size(a, 2), 1]);
                set(h2, 'CData', b);
            end
        else
            set(h2, 'CData', img);
        end
        slidshown = ~slidshown;
    end

    function callback_GENEEXPR(~, ~)
        %answer = questdlg('Show expression of single or mulitple genes?',...
        %    'Single/Multiple Genes','Single','Multiple','Cancel','Single');
        % answer='Multiple';
        % switch answer
        %     case 'Single'
        %         [gsorted]=gui.i_sortgenenames(sce);
        %         if isempty(gsorted), return; end
        %         [indx,tf] = listdlg('PromptString',{'Select a gene:','',''},...
        %             'SelectionMode','single','ListString',gsorted);
        %         if tf==1
        %             i_cascadeexpr(sce,gsorted(indx),xy);
        %         end
        %     case 'Multiple'
        [glist] = gui.i_selectngenes(sce);
        if isempty(glist)
            helpdlg('No gene selected.', '');
            return;
        end
        for k = 1:length(glist)
            i_cascadeexpr(sce, glist(k), xy, k);
        end
        %         otherwise
        %             return;
        % end
    end

    function [f] = i_cascadeexpr(sce, g, xy, k)
        if nargin < 4, k = 1; end
        f = figure('Visible', 'off');
        tb = uitoolbar('Parent', f);
        P = get(f, 'Position');
        set(f, 'Position', [P(1) - 120, P(2), round(0.75*P(3)), round(0.75*P(4))]);
        ix = find(sce.g == g);
        cx = sce.X(ix(1), :);
        h = scatter(xy(:, 1), xy(:, 2), dsize*0.75, cx, 'filled');
        axis ij
        gui.i_setautumncolor(cx, 'parula');
        title(g)
        cb = colorbar;
        cb.Label.String = 'Expression Level';
        pkg.i_addbutton2fig(tb, 'off', {@i_RescaleExpr, h, cb}, 'IMG00067.GIF', 'Scale expression level using log2-transformation');
        pkg.i_addbutton2fig(tb, 'off', {@i_ResetExpr, h, cx, cb}, 'IMG00074.GIF', 'Reset expression level to UMI');
        pkg.i_addbutton2fig(tb, 'on', {@i_genecards, g}, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
        pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
        P = get(f, 'Position');
        set(f, 'Position', [P(1) - 20 * k, P(2) - 20 * k, P(3), P(4)]);
        set(f, 'visible', 'on');
        drawnow;
        box on
        grid on
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcnx2, cx};
    end

    function i_RescaleExpr(~, ~, h, cb)
        set(h, 'CData', log2(1+get(h, 'CData')));
        cb.Label.String = 'log2(UMI+1)';
    end

    function i_ResetExpr(~, ~, h, v, cb)
        set(h, 'CData', v);
        cb.Label.String = 'Expression Level';
    end

    function i_seth1cdata(iscontinuous, ttxt)
        cx = c;
        if nargin < 1
            if length(unique(cx)) <= 20
                iscontinuous = true;
            else
                iscontinuous = false;
            end
        end
        if nargin < 2, ttxt = []; end
        set(h1, 'CData', cx);
        colormap default
        if iscontinuous
            a = colormap('parula');
            a(1, :) = [.8, .8, .8];
        else
            a = lines(length(unique(cx)));
        end
        colormap(hAx, a);
        if ~isempty(ttxt)
            cb = colorbar;
            cb.Label.String = strrep(ttxt, '_', '\_');
        else
            colorbar off
        end
    end

    function i_addbutton(Handle, sepTag, callbackFnc, imgFil, tooltipTxt)
        if nargin < 4, sepTag = 'off'; end
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc = str2func(callbackFnc);
        end
        pt = uipushtool(Handle, 'Separator', sepTag);
        pt.CData = i_get_ptImage(imgFil);
        pt.Tooltip = tooltipTxt;
        %if iscell(callbackFnc)
        %    pt.ClickedCallback = {callbackFnc{1},callbackFnc{2}};
        %else
        %    pt.ClickedCallback = callbackFnc;
        %end
        pt.ClickedCallback = callbackFnc;
    end


    function callback_Brush4MarkersLASSO(~, ~)
        ptsSelected = logical(h1.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return;
        end
        [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce);
        if ~letdoit, return; end
        [numfig] = gui.i_inputnumg;
        if isempty(numfig), return; end
        fw = gui.gui_waitbar;
        y = double(ptsSelected);
        sce.c = 1 + ptsSelected;
        X = sce.X';
        try
            if issparse(X), X = full(X); end
            [B] = lasso(X, y, 'DFmax', numfig*3, 'MaxIter', 1e3);
        catch ME
            gui.gui_waitbar(fw, true);
            errordlg(ME.message);
            return;
        end

        [~, ix] = min(abs(sum(B > 0)-numfig));
        b = B(:, ix);
        idx = b > 0;
        gui.gui_waitbar(fw);

        if ~any(idx)
            errordlg('No marker gene found')
            return;
        else
            markerlist = sce.g(idx);
            [~, jx] = sort(b(idx), 'descend');
            markerlist = markerlist(jx);
            axesh = FigureHandle.findobj('type', 'Axes');
            [ax, bx] = view(axesh);
            for kk = 1:length(markerlist)
                gui.i_cascadefig(sce, markerlist(end-(kk - 1)), ax, bx, kk);
            end
        end
        fprintf('%d marker genes: ', length(markerlist));
        fprintf('%s ', markerlist)
        fprintf('\n')
    end

    function [OKPressed] = callback_SAVESTE(~, ~)
        ste.xy = xy;
        if ~any(strcmp(sce.list_cell_attributes, 'ste.xy'))
            sce.list_cell_attributes{end+1} = 'ste.xy';
            sce.list_cell_attributes{end+1} = xy;
        end
        ste.sce = sce;
        OKPressed = false;
        answer = questdlg('Export & save data to:', '', ...
            'Workspace', 'MAT file', 'Workspace');
        switch answer
            case 'Workspace'
                labels = {'Save STE to variable named:', ...
                    'Save STE.SCE to variable named:', 'Save STE.XY to variable named:'};
                vars = {'ste', 'sce', 'xy'};
                values = {ste, ste.sce, ste.xy};
                [~, OKPressed] = export2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', logical([1, 0, 0]));
            case 'MAT file'
                [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as');
                if isequal(file, 0) || isequal(path, 0)
                    return;
                else
                    filename = fullfile(path, file);
                    fw = gui.gui_waitbar;
                    save(filename, 'ste', '-v7.3');
                    gui.gui_waitbar(fw);
                    OKPressed = true;
                end
            otherwise
                return;
        end
    end

    function callback_BAYESSPACE(~, ~)
        [glist] = gui.i_selectngenes(sce);
        if isempty(glist)
            helpdlg('No feature genes selected.', '');
            return;
        end
        fw = gui.gui_waitbar;
        try
            [T, X, ste1] = st.run.BayesSpace(ste, glist);
        catch ME
            gui.gui_waitbar(fw, true);
            errordlg(ME.message);
            return;
        end
        gui.gui_waitbar(fw);
        try
            if ~isempty(T)
                %figure;
                %scatter(T.row,T.col,15,T.spatial_cluster,'filled');
                %axis ij;
                if ~isempty(ste1)
                    if ~(ismcc || isdeployed)
                        labels = {'Save STE to variable named:', ...
                            'Save BayesSpace output to variable named:'};
                        vars = {'ste_enhanced', 'T_bayesspace'};
                        values = {ste1, T};
                        [~, ~] = export2wsdlg(labels, vars, values, ...
                            'Save Data to Workspace', logical([1, 1]));
                    end
                    answer = questdlg('Show STE with enhancement of spatial resolution?', '');
                    if strcmp(answer, 'Yes')
                        st_scatter_ste(ste1, 'c', T.spatial_cluster, 'dsize', 10);
                    end
                end
                if ~isempty(X)
                    answer = questdlg('Show feature genes in spatially-enhanced STE?', '');
                    if strcmp(answer, 'Yes')
                        sce1 = SingleCellExperiment(X, glist);
                        % sce1.struct_cell_clusterings.bayesspace=T.spatial_cluster;
                        for k = 1:length(glist)
                            i_cascadeexpr(sce1, glist(k), [T.row, T.col], k);
                        end
                    end
                end
            end
        catch ME
            errordlg(ME.message);
        end
    end

    % --------------------------

    function [txt] = i_myupdatefcnx1(~, event_obj)
        % pos = event_obj.Position;
        txt = '';
        idx = event_obj.DataIndex;
        if idx <= sce.NumCells
            [cx, cL] = grp2idx(c);
            txt = cL(cx(idx));
        end
    end

    function [txt] = i_myupdatefcnx2(~, event_obj, cx)
        idx = event_obj.DataIndex;
        txt = cx(idx);
    end

    function i_genecards(~, ~, g)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', g));
    end

    function [co] = i_getcolormaps(n)
        mycmap = pkg.i_mycolormap(n);
        co = {parula(n), lines(n), summer(n), ...
            jet(n), copper(n), winter(n), hsv(n), ...
            mycmap, ...
            gui.linspecer(min([n, 12]), 'qualitative'), ...
            gui.linspecer(n, 'sequential'), ...
            gui.linspecer(n, 'red'), gui.linspecer(n, 'gray'), gui.linspecer(n, 'green')};
        indx = randi(length(co));
        co = co{indx};
    end

    function [i1, i2] = i_select2grps(thisc)
        i1 = 0;
        i2 = 0;
        if numel(unique(thisc)) == 1
            warndlg("Cannot compare with an unique group");
            return;
        end
        [~, cLi] = grp2idx(thisc);
        listitems = string(cLi);
        n = length(listitems);
        if n < 2
            errordlg('Need at least two groups.');
            return;
        end
        [indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
            'SelectionMode', 'multiple', ...
            'ListString', listitems, ...
            'InitialValue', [n - 1, n]);
        if tfx == 1
            if numel(indxx) ~= 2
                errordlg('Please select 2 groups');
                return;
            end
            i1 = indxx(1);
            i2 = indxx(2);
        else
            return;
        end
    end

    function [txt] = i_myupdatefcnx(~, event_obj)
        % pos = event_obj.Position;
        idx = event_obj.DataIndex;
        if idx <= ste.NumSpots
            txt = string(c(idx));
        else
            txt = '';
        end
    end

    function i_addmenu(menuHdl, sepTag, callbackFnc, tooltipTxt)
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc = str2func(callbackFnc);
        end
        uimenu(menuHdl, 'Text', tooltipTxt, ...
            'Separator', sepTag, ...
            'Callback', callbackFnc);
    end

    function i_addbutton_toggle(toolbarHdl, sepTag, callbackFnc, tooltipTxt)
        imgFil = callbackFnc{3};
        %if ischar(callbackFnc{1}) || isstring(callbackFnc{1})
        %    callbackFnc=str2func(callbackFnc{1});
        %end
        if sepTag == 1
            septag = 'on';
        else
            septag = 'off';
        end
        if toolbarHdl == 0
            barhandle = MitoolbarHandle;
        elseif toolbarHdl == 1
            barhandle = UitoolbarHandle;
        end
        pt = uitoggletool(barhandle, 'Separator', septag);
        pt.CData = i_get_ptImage(imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
    end

    function togglebtfun(src, ~, func, imgFil1, imgFil2, actiondelay)
        if nargin < 6, actiondelay = true; end
        if src.State == "off"
            imgFil = imgFil1;
        elseif src.State == "on"
            imgFil = imgFil2;
        end
        src.CData = i_get_ptImage(imgFil);
        if actiondelay
            if src.State == "off"
                func(src);
            else
                uiwait(helpdlg('To execute the function, click the button again or locate and click the same button in the toolbar above. Hover over the button to view a description of its function.', ''));
                end
            else
                func(src);
            end
        end

            function [ptImage] = i_get_ptImage(imgFil)
                try
                    [imgx, mapx] = imread(fullfile(mfolder, 'resources', imgFil));
                    ptImage = ind2rgb(imgx, mapx);
                catch
                    ptImage = rand(16, 16, 3);
                end
        end

        end
