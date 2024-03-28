function stgeatool(sce, img, xy)

    %cdgea_st;
    olddir = pwd();
    cdgea;
    cd(olddir);
    
    if ~(ismcc || isdeployed)
        if ~exist('grp2idx.m', 'file')
            errordlg('Statistics and Machine Learning Toolbox is required.');
            return;
        end
        if ~exist('scgeatool.m', 'file')
            errordlg('scGEAToolbox is required. Visit https://github.com/jamesjcai/scGEAToolbox to learn how to install.');
            return;
        end
    end
    
    if usejava('jvm') && ~feature('ShowFigureWindows')
        error('MATLAB is in a text mode. This function requires a GUI-mode.');
        end
        promotesave = false;
        if nargin < 1
            list = {'STE Data File (*.mat)...', ...
                '10x Visium ''outs'' Folder...', ...
                '----------------------------------', ...
                'GEO Accession Number...', ...
                '----------------------------------', ...
                'Load STE Variable from Workspace...', ...
                'Load Example Data...'};
    
            [indx, tf] = listdlg('ListString', list, ...
                'SelectionMode', 'single', ...
                'PromptString', {'Select an input data type:'}, ...
                'ListSize', [230, 200], ...
                'Name', 'STGEATOOL', 'InitialValue', length(list));
            if tf ~= 1, return; end
            ButtonName = list{indx};
    
            switch ButtonName
                case 'STE Data File (*.mat)...'
                    promotesave = false;
                    [fname, pathname] = uigetfile( ...
                        {'*.mat', 'STE Data Files (*.mat)'; ...
                        '*.*', 'All Files (*.*)'}, ...
                        'Pick a STE Data File');
                    % if ~(fname), return; end
                    if isequal(fname, 0), return; end
                    stefile = fullfile(pathname, fname);
                    try
                        fw = gui.gui_waitbar;
                        load(stefile, 'ste');
                    catch ME
                        gui.gui_waitbar(fw, true);
                        errordlg(ME.message);
                        return;
                    end
                    gui.gui_waitbar(fw);
    
                case '10x Visium ''outs'' Folder...'
                    selpath = uigetdir;
                    if selpath == 0, return; end
                    try
                        fw = gui.gui_waitbar;
    
                        [ste] = st_read10xdir(selpath);
                        gui.gui_waitbar(fw);
                    catch ME
                        gui.gui_waitbar(fw, true);
                        errordlg(ME.message);
                        return;
                    end
                    %if ~ftdone, errordlg('Input Error'); return; end
                    %metainfo=sprintf("Source: %s",selpath);
                    %sce=sce.appendmetainfo(metainfo);
                    % ste=SpatialTranscriptomicsExperiment(sce,img,xy);
                case 'GEO Accession Number...'
                    acc = inputdlg({'Input number (e.g., GSM5764426, GSM5764424 or GSM5213483):'}, ...
                        'GEO Accession', [1, 50], {'GSM4565826'});
                    if isempty(acc), return; end
                    acc = deblank(acc{1});
                    if isempty(acc) || ~strlength(acc) > 4, return; end
                    if strlength(acc) > 4 && ~isempty(regexp(acc, 'G.+', 'once'))
                        accv = unique(strsplit(acc, {',', ';', ' '}), 'stable');
                        if length(accv) > 1
                            %                         dmanswer=questdlg('Download and merge data sets?',...
                            %                             '','Yes','Cancel','Yes');
                            %                         if ~strcmp(dmanswer,'Yes'), return; end
                            %                         [sce]=pkg.pipeline_multisamplesmerge(accv);
                        else
                            try
                                fw = gui.gui_waitbar;
                                % [sce]=sc_readgeoaccession(acc);
                                [ste] = st_readgeoaccession(acc);
                                gui.gui_waitbar(fw);
                            catch ME
                                gui.gui_waitbar(fw);
                                errordlg(ME.message);
                                return;
                            end
                        end
                    end
                    % ste=SpatialTranscriptomicsExperiment(sce,img,xy);
                case 'Load STE Variable from Workspace...'
                    a = evalin('base', 'whos');
                    b = struct2cell(a);
                    valididx = ismember(b(4, :), 'SpatialTranscriptomicsExperiment');
                    if isempty(valididx)
                        helpdlg('No STE in the Workspace.', '');
                        return;
                    end
                    a = a(valididx);
                    [indx, tf] = listdlg('PromptString', {'Select STE variable:'}, ...
                        'liststring', b(1, valididx), 'SelectionMode', 'single');
                    if tf == 1
                        ste = evalin('base', a(indx).name);
                    else
                        return;
                    end
                case 'Load Example Data...'
                    answerstruced = questdlg('Load processed or raw data?', ...
                        '', 'Processed', 'Raw', 'Cancel', 'Processed');
                    if ~(strcmp(answerstruced, 'Processed') || strcmp(answerstruced, 'Raw'))
                        return;
                    end
                    switch answerstruced
                        case 'Processed'
                            f = 'new_test.mat';
                        case 'Raw'
                            f = 'testSte_unaligned.mat';
                    end
                    pw1 = fileparts(mfilename('fullpath'));
                    fprintf('Loading STE Data File example_data/%s...', f);
                    tic;
                    file1 = fullfile(pw1, 'example_data', f);
                    if ~exist(file1, "file")
                        errordlg("Example data file does not exist.");
                        return;
                    end
                    load(file1, 'ste');
                    %if strcmp(answerstruced,'Raw')
                    %    sce=SingleCellExperiment(sce.X,sce.g);
                    %end
                    fprintf('Done.\n');
                    toc;
                otherwise
                    return;
            end
        else
            if isa(sce, 'SpatialTranscriptomicsExperiment')
                st_scatter_ste(sce);
                return;
            end
            if nargin < 1 || isempty(sce) || ~isa(sce, 'SingleCellExperiment')
                error('Invalid inputs.');
            end
            if nargin < 3 || isempty(xy)
                xy = sce.s(:, 1:2);
            end
            if nargin < 2 || isempty(img)
                img = zeros(128, 128);
                img(32:96, 32:96) = rand(96-32+1, 96-32+1);
            end
            try
                ste = SpatialTranscriptomicsExperiment(sce, xy, img);
            catch
                stgeatool;
            end
        end
        if isempty(ste), return; end
        try
            st_scatter_ste(ste);
        catch ME
            disp(ME.identifier);
            errordlg(ME.message);
        end
    end
    