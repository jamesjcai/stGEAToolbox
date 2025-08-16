function spaceRangerViewer()
    % spaceRangerViewer: GUI to view 10x Genomics Space Ranger v4 outputs
    %
    % Shows histology image with aligned spot overlays.

    % Create UI figure
    fig = uifigure('Name','Space Ranger Viewer','Position',[100 100 900 700]);
    ax = uiaxes(fig,'Position',[50 80 800 600]);
    btn = uibutton(fig,'push','Text','Load Space Ranger Folder',...
        'Position',[50 20 200 40],'ButtonPushedFcn',@(src,evt)onLoad(ax));

end

function onLoad(ax)
    folder = uigetdir(pwd,'Select Space Ranger output folder');
    if folder == 0
        return;
    end
    
    spatialDir = fullfile(folder,'spatial');
    
    %--- Load image ---
    imgFile = fullfile(spatialDir,'tissue_hires_image.png');
    if ~isfile(imgFile)
        uialert(ax.Parent,'tissue_hires_image.png not found','Error');
        return;
    end
    img = imread(imgFile);

    %--- Load scalefactors ---
    scaleFile = fullfile(spatialDir,'scalefactors_json.json');
    if ~isfile(scaleFile)
        uialert(ax.Parent,'scalefactors_json.json not found','Error');
        return;
    end
    txt = fileread(scaleFile);
    S = jsondecode(txt);
    scale = S.tissue_hires_scalef;  % scale for hires image

    %--- Load spot positions (parquet or csv) ---
    posFile = fullfile(spatialDir,'tissue_positions.parquet');
    if isfile(posFile)
        T = parquetread(posFile);
    else
        posFile = fullfile(spatialDir,'tissue_positions_list.csv');
        if ~isfile(posFile)
            uialert(ax.Parent,'No tissue_positions file found','Error');
            return;
        end
        opts = detectImportOptions(posFile);
        T = readtable(posFile,opts);
    end

    % Normalize variable names (parquet vs csv differ)
    if any(strcmpi(T.Properties.VariableNames,'pxl_col_in_fullres'))
        x = T.pxl_col_in_fullres * scale;
        y = T.pxl_row_in_fullres * scale;
    elseif width(T) >= 3
        % CSV format: barcode, in_tissue, row, col, pxl_col_in_fullres, pxl_row_in_fullres
        x = T{:,5} * scale;
        y = T{:,6} * scale;
    else
        error('Unsupported tissue_positions format.');
    end
    
    %--- Plot ---
    imshow(img,'Parent',ax); hold(ax,'on');
    scatter(ax,x,y,20,'r','filled','MarkerFaceAlpha',0.4);
    hold(ax,'off');
    title(ax,'Histology with Spots');
end
