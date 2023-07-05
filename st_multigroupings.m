function [f0]=st_multigroupings(s_n,sce)

answer=questdlg('Cluster algorithm?','','K-means','SC3','K-means');
switch answer
    case 'K-means'
        k=gui.i_inputnumk(5);
        if isempty(k), return; end
        fw=gui.gui_waitbar;
        try
        sce2=sce.embedcells('tSNE',true);
        sce2.c=sc_cluster_s(sce2.s,k);
        gui.gui_waitbar(fw);
        catch ME
            gui.gui_waitbar(fw,true);
            errordlg(ME.message);
            return;
        end
    case 'SC3'
        k=gui.i_inputnumk(5);
        if isempty(k), return; end
        fw=gui.gui_waitbar;
        try
        sce2=sce.embedcells('tSNE',true);
        sce2.c=sc_cluster_x(sce2.X,k);
        gui.gui_waitbar(fw);
                catch ME
            gui.gui_waitbar(fw,true);
            errordlg(ME.message);
            return;
        end
    otherwise
        return;
end


f0=figure('Visible',false);

ax1=subplot(1,2,1);
%h1=scatter(s_n(:,1),s_n(:,2),10,sum(sce.X),'filled');
h1=scatter(s_n(:,1),s_n(:,2),10,sce2.c,'filled');
dt = datacursormode(f0);
dt.UpdateFcn = {@i_myupdatefcnx12};
box on
view(2);

ax2=subplot(1,2,2);
% h2=scatter(sce.s(:,1),sce.s(:,2),10,sce.c,'filled');
% h2=gui.i_gscatter3(sce.s,sce.c,1,1);
%h2=scatter3(sce2.s(:,1),sce2.s(:,2),sce2.s(:,3),10,sce.c,'filled');
h2=gui.i_gscatter3(sce2.s,sce2.c,1,1);
box on
view(3);

dt = datacursormode(f0);
dt.UpdateFcn = {@i_myupdatefcnx12};
sgtitle(sce.title);


% kc1=numel(unique(c1));
% colormap(ax1,pkg.i_mycolorlines(kc1));
% kc2=numel(unique(c2));
% colormap(ax2,pkg.i_mycolorlines(kc2));

%evalin('base','h=findobj(gcf,''type'',''axes'');');
%evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
%rotate3d(f0,'on');
f0.Position(3)=f0.Position(3)*2;

hBr=brush(f0);
hBr.ActionPostCallback = {@onBrushAction,h1,h2};


% tb = uitoolbar(f0);
set(f0,'Visible',true);


function [txt] = i_myupdatefcnx12(~, ~)
    txt='';
end

function onBrushAction(~,event,h1,h2)
    if isequal(event.Axes.Children,h1)
        h2.BrushData=h1.BrushData;
    elseif isequal(event.Axes.Children,h2)
        h1.BrushData=h2.BrushData;
    end
end
end


