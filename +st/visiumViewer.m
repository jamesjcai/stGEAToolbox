function visiumViewer
% visiumViewer
% UI app to load 10x Genomics Space Ranger (v4.x) outputs and overlay spots on slide image.
% - Choose an outs/ folder produced by spaceranger count.
% - Loads spatial images (hires/lowres/aligned) + tissue_positions CSV + scalefactors JSON.
% - Projects spot centers to chosen image resolution and overlays as circles or filled dots.
%
% Tested on R2021b+. Requires no toolboxes (Image Processing Toolbox optional).
%
% UI controls:
%   • "Load outs/..." button – select the Space Ranger outputs directory.
%   • Image dropdown – choose which image to display (hires/lowres/aligned).
%   • "In-tissue only" – hide off-tissue spots.
%   • "Draw circles (outline)" – draw true-diameter circles (slower); off = filled dots (faster).
%   • Size slider – scale marker/circle size.
%   • Alpha slider – transparency of spots.

    % --- State ---
    S = struct();
    S.outsDir = '';
    S.img = [];
    S.imgName = '';
    S.scalef = struct('tissue_hires_scalef', [], 'tissue_lowres_scalef', [], ...
                      'spot_diameter_fullres', [], 'fiducial_diameter_fullres', []);
    S.pos = table();   % parsed positions with columns: barcode, in_tissue, x_fullres, y_fullres
    S.axesImg = [];
    S.axesOverlay = [];
    S.imHandle = [];
    S.spotHandle = [];
    S.imageChoices = {};
    S.drawCircles = false;

    % --- UI ---
    fig = uifigure('Name','Visium / Space Ranger v4 Viewer','Position',[100 100 1100 750]);
    gl  = uigridlayout(fig,[6 6]);
    gl.RowHeight = {40, 40, 40, 30, '1x', 30};
    gl.ColumnWidth = {160, 160, 160, 160, 180, '1x'};

    btnLoad   = uibutton(gl,'Text','Load outs/...','ButtonPushedFcn',@onLoad);
    btnLoad.Layout.Row = 1; btnLoad.Layout.Column = [1 2];

    ddImage   = uidropdown(gl,'Items',{'(none)'},'Value','(none)','ValueChangedFcn',@onImageChange);
    ddImage.Layout.Row = 1; ddImage.Layout.Column = [3 4];
    ddImage.Tooltip = 'Choose which spatial image to display';

    cbInTissue = uicheckbox(gl,'Text','In-tissue only','Value',true,'ValueChangedFcn',@onRedraw);
    cbInTissue.Layout.Row = 1; cbInTissue.Layout.Column = 5;

    cbCircles  = uicheckbox(gl,'Text','Draw circles (outline)','Value',false,'ValueChangedFcn',@onCircles);
    cbCircles.Layout.Row = 1; cbCircles.Layout.Column = 6;


    lbl = uilabel(gl,'Text','Size','HorizontalAlignment','right');
    lbl.Layout.Row = 2;
    lbl.Layout.Column = 1;

    % uilabel(gl,'Text','Size','HorizontalAlignment','right').Layout.Row = 2; ans.Layout.Column = 1; %#ok<*NASGU>
    sSize = uislider(gl,'Limits',[0.25 4],'Value',1,'ValueChangedFcn',@onRedraw,'MajorTicks',[]);
    sSize.Layout.Row = 2; sSize.Layout.Column = [2 4];

    lbl2 = uilabel(gl,'Text','Alpha','HorizontalAlignment','right');
    lbl2.Layout.Row = 2;
    lbl2.Layout.Column = 5;
    
    % uilabel(gl,'Text','Alpha','HorizontalAlignment','right').Layout.Row = 2; ans.Layout.Column = 5;
    sAlpha = uislider(gl,'Limits',[0.05 1],'Value',0.7,'ValueChangedFcn',@onRedraw,'MajorTicks',[]);
    sAlpha.Layout.Row = 2; sAlpha.Layout.Column = [5 6];

    msg = uitextarea(gl,'Editable','off','Value',{'Load an outs/ folder to begin...'});
    msg.Layout.Row = 3; msg.Layout.Column = [1 6];

    ax = uiaxes(gl);  % we'll use one axes and two layers (image + overlay)
    ax.Layout.Row = [4 5]; ax.Layout.Column = [1 6];
    axis(ax,'ij'); axis(ax,'image'); ax.Visible = 'off';
    hold(ax,'on');
    S.axesImg = ax;
    S.axesOverlay = ax;

    % Footer
    footer = uilabel(gl,'Text','10x Visium/Visium HD • Space Ranger v4 viewer (MATLAB)','HorizontalAlignment','center');
    footer.Layout.Row = 6; footer.Layout.Column = [1 6];

    % Nested callbacks capture S by reference via guidata-ish pattern using appdata
    setappdata(fig,'S',S);

    % ----------- Callbacks -----------
    function onLoad(~,~)
        S = getappdata(fig,'S');
        d = uigetdir(pwd,'Select Space Ranger outs/ directory');
        if isequal(d,0), return; end
        S.outsDir = d;
        try
            [S.img, S.imageChoices, S.imgName] = loadBestImage(S.outsDir);
            S.scalef = readScalefactors(fullfile(S.outsDir,'spatial','scalefactors_json.json'));
            S.pos    = readPositions(fullfile(S.outsDir,'spatial'));
            ddImage.Items = S.imageChoices;
            if ~isempty(S.imgName) && any(strcmp(S.imgName, S.imageChoices))
                ddImage.Value = S.imgName;
            else
                ddImage.Value = S.imageChoices{1};
            end
            ax.Visible = 'on';
            cla(ax); S.imHandle = imshow(S.img,'Parent',ax); hold(ax,'on');
            title(ax, sprintf('%s', ddImage.Value), 'Interpreter','none');
            setappdata(fig,'S',S);
            drawSpots();
            setMsg({'Loaded outs/:', S.outsDir, ...
                    sprintf('Positions: %d barcodes', height(S.pos)), ...
                    sprintf('Image: %s', ddImage.Value)});
        catch ME
            setMsgErr(ME);
        end
    end

    function onImageChange(~,~)
        S = getappdata(fig,'S');
        try
            if isempty(S.outsDir), return; end
            [S.img, ~, S.imgName] = loadSpecificImage(S.outsDir, ddImage.Value);
            cla(ax); S.imHandle = imshow(S.img,'Parent',ax); hold(ax,'on');
            title(ax, sprintf('%s', S.imgName), 'Interpreter','none');
            setappdata(fig,'S',S);
            drawSpots();
        catch ME
            setMsgErr(ME);
        end
    end

    function onCircles(src,~)
        S = getappdata(fig,'S');
        S.drawCircles = logical(src.Value);
        setappdata(fig,'S',S);
        drawSpots();
    end

    function onRedraw(~,~)
        drawSpots();
    end

    % ----------- Drawing -----------
    function drawSpots
        S = getappdata(fig,'S');
        if isempty(S.img) || isempty(S.pos), return; end
        if isgraphics(S.spotHandle), delete(S.spotHandle); end

        % pick scale factor based on selected image
        scale = pickScaleFactor(S.scalef, ddImage.Value);

        % map fullres coords to current image
        x = S.pos.x_fullres * scale;
        y = S.pos.y_fullres * scale;

        % filter in-tissue
        if cbInTissue.Value && any(strcmpi(S.pos.Properties.VariableNames,'in_tissue'))
            mask = S.pos.in_tissue > 0;
            x = x(mask); y = y(mask);
        end

        % choose rendering mode
        alpha = sAlpha.Value;
        sizeScale = sSize.Value;

        % estimate spot diameter in current image pixels
        d_full = S.scalef.spot_diameter_fullres;
        if isempty(d_full) || isnan(d_full), d_full = 100; end % fallback
        d_img  = d_full * scale * sizeScale;
        r_img  = max(d_img/2, 1);

        if S.drawCircles
            % Draw outline circles (slower for many spots)
            % Try to batch with line objects
            theta = linspace(0,2*pi,64);
            X = x(:)' + r_img.*cos(theta);
            Y = y(:)' + r_img.*sin(theta);
            % For variable radii (if later needed) use expansion; here single r_img scalar
            S.spotHandle = plot(S.axesOverlay, X, Y, '-', 'LineWidth', 0.75, 'Color', [1 0 0 alpha]);
        else
            % Fast filled scatter; map diameter to MarkerSize (points^2). Approximate:
            % MarkerSize in scatter is area in points^2; relate pixels to points using heuristic.
            % We'll choose a visually reasonable mapping:
            ms = max( (d_img.^2) / 30, 8 ); % tweakable
            S.spotHandle = scatter(S.axesOverlay, x, y, ms, 'filled', ...
                                   'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
        end

        drawnow limitrate;
    end

    % ----------- Helpers -----------
    function setMsg(lines)
        msg.Value = cellstr(lines);
    end

    function setMsgErr(ME)
        msg.Value = {['Error: ' ME.message], ' ', 'Stack (top):', ...
                     sprintf('%s (line %d)', ME.stack(1).name, ME.stack(1).line)};
    end

    function [img, choices, chosen] = loadBestImage(outsDir)
        spatialDir = fullfile(outsDir,'spatial');
        candidates = { ...
            'tissue_hires_image.png', ...
            'tissue_lowres_image.png', ...
            'aligned_tissue_image.jpg', ...
            'cytassist_image.tif', ...
            'cytassist_image.tiff' ...
        };
        existing = candidates(cellfun(@(f) exist(fullfile(spatialDir,f),'file')>0, candidates));
        if isempty(existing)
            error('No spatial images found in %s', spatialDir);
        end
        % Prefer hires if available
        pick = existing{1};
        img = imread(fullfile(spatialDir, pick));
        choices = existing;
        chosen = pick;
    end

    function [img, choices, chosen] = loadSpecificImage(outsDir, name)
        [~, choices, ~] = loadBestImage(outsDir); %#ok<ASGLU>
        spatialDir = fullfile(outsDir,'spatial');
        p = fullfile(spatialDir, name);
        if exist(p,'file')==0
            [img, choices, chosen] = loadBestImage(outsDir);
            return;
        end
        img = imread(p);
        chosen = name;
        % Rebuild choices list in case on first call
        all = { ...
            'tissue_hires_image.png', ...
            'tissue_lowres_image.png', ...
            'aligned_tissue_image.jpg', ...
            'cytassist_image.tif', ...
            'cytassist_image.tiff' ...
        };
        choices = all(cellfun(@(f) exist(fullfile(spatialDir,f),'file')>0, all));
        if isempty(choices), choices = {name}; end
    end

    function scalef = readScalefactors(jsonPath)
        scalef = struct('tissue_hires_scalef',[], 'tissue_lowres_scalef',[], ...
                        'spot_diameter_fullres',[], 'fiducial_diameter_fullres',[]);
        if exist(jsonPath,'file')==0, return; end
        J = jsondecode(fileread(jsonPath));
        flds = fieldnames(scalef);
        for i=1:numel(flds)
            if isfield(J, flds{i})
                scalef.(flds{i}) = J.(flds{i});
            end
        end
    end

    function scale = pickScaleFactor(scalef, imageName)
        % Positions are in FULLRES pixel space.
        % To project to current image, multiply by the appropriate scale factor.
        if contains(lower(imageName),'lowres') && ~isempty(scalef.tissue_lowres_scalef)
            scale = scalef.tissue_lowres_scalef;
        elseif contains(lower(imageName),'hires') && ~isempty(scalef.tissue_hires_scalef)
            scale = scalef.tissue_hires_scalef;
        else
            % aligned/cytassist images don't have explicit scale; fallback to hires if present
            if ~isempty(scalef.tissue_hires_scalef), scale = scalef.tissue_hires_scalef;
            elseif ~isempty(scalef.tissue_lowres_scalef), scale = scalef.tissue_lowres_scalef;
            else, scale = 1; % last resort
            end
        end
    end

    function T = readPositions(spatialDir)
        % Space Ranger v4 may write either tissue_positions.csv (w/ header)
        % or tissue_positions_list.csv (older, no header).
        cand = {'tissue_positions.csv','tissue_positions_list.csv'};
        file = '';
        for i=1:numel(cand)
            p = fullfile(spatialDir, cand{i});
            if exist(p,'file'), file = p; break; end
        end
        if isempty(file)
            error('No tissue_positions CSV found in %s', spatialDir);
        end

        % Try reading with header first
        opts = detectImportOptions(file, 'Delimiter',',');
        T0 = readtable(file, opts);

        if any(strcmpi(T0.Properties.VariableNames,'pxl_col_in_fullres'))
            % Headered format
            T = table();
            T.barcode = T0.barcode;
            T.in_tissue = T0.in_tissue;
            T.x_fullres = double(T0.pxl_col_in_fullres);
            T.y_fullres = double(T0.pxl_row_in_fullres);
        else
            % Possibly headerless legacy: 6 columns
            % 1=barcode 2=in_tissue 3=array_row 4=array_col 5=pxl_col_in_fullres 6=pxl_row_in_fullres
            if width(T0) < 6
                % Force no header read
                Traw = readtable(file,'Delimiter',',','ReadVariableNames',false,'TextType','string');
                if width(Traw) < 6
                    error('Unrecognized tissue_positions format.');
                end
                T = table();
                T.barcode   = Traw.Var1;
                T.in_tissue = double(Traw.Var2);
                T.x_fullres = double(Traw.Var5);
                T.y_fullres = double(Traw.Var6);
            else
                % Header exists but names differ (be defensive)
                vn = T0.Properties.VariableNames;
                col_x = find(contains(lower(vn),'pxl_col'));
                col_y = find(contains(lower(vn),'pxl_row'));
                col_bar = find(contains(lower(vn),'barcode'));
                col_tis = find(contains(lower(vn),'in_tissue'));
                if isempty(col_x) || isempty(col_y)
                    error('Could not identify pixel columns in %s', file);
                end
                T = table();
                T.barcode = T0.(vn{col_bar(1)});
                T.in_tissue = double(T0.(vn{col_tis(1)}));
                T.x_fullres = double(T0.(vn{col_x(1)}));
                T.y_fullres = double(T0.(vn{col_y(1)}));
            end
        end
    end
end
