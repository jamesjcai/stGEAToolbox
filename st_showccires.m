function st_showccires(T,~,cx,ste)
import mlreportgen.ppt.*
%  [T,RES,cx]=st_ridgecci(ste,c,x1,x2,K,methodid);
% RES{1}=nv,vv}

FigureHandle=figure;
    answer=questdlg('Show top 5 gene pairs'' expression?');
    if strcmp(answer,'Yes')
        OUTppt=[tempname,'.pptx'];
        ppt = Presentation(OUTppt);
        open(ppt);

        images={};

        img0=[tempname,'.png'];
        saveas(FigureHandle,img0);
        slide4 = add(ppt,'Title and Content');
        replace(slide4,'Content',Picture(img0));
        images = [images {img0}];

        c=cx{1};
        % i_seth1cdata(false,'Ridge Spots (1-source, 2-target)');
        %FigureHandle=figure;
        scatter(ste.xy(:,1),ste.xy(:,2),5,c,'filled'); 
        axis ij; 
        colormap(lines(length(unique(c))));

        img1=[tempname,'.png'];
        saveas(FigureHandle,img1);
%           slide4 = add(ppt,'Title and Content');
%           replace(slide4,'Content',Picture(img1));
        images = [images {img1}];
        
        c=cx{2};
        %i_seth1cdata(false,'Ridge Spots (1-source, 2-target)');
        scatter(ste.xy(:,1),ste.xy(:,2),5,c,'filled'); 
        axis ij; 
        colormap(lines(length(unique(c))));
        
        img2=[tempname,'.png'];
        saveas(FigureHandle,img2);
%            slide4 = add(ppt,'Title and Content');
%            replace(slide4,'Content',Picture(img1));
        images = [images {img2}];

        slide4 = add(ppt,'Comparison');
        replace(slide4,'Left Content',Picture(img1));
        replace(slide4,'Right Content',Picture(img2));
        
%             f=figure('Visible','off');
%             scatter(xy(:,1),xy(:,2),10,cx,'filled');
%             axis ij;
%             colormap(lines(length(unique(cx))));
%             img1=[tempname,'.png'];
%             saveas(f,img1);
%             slide4 = add(ppt,'Title and Content');
%             replace(slide4,'Content',Picture(img1));
%             images = [images {img1}];

        for k=1:2:10
            j=(k+1)/2;
            img1=[tempname,'.png'];
            img2=[tempname,'.png'];
            images = [images {img1}];
            images = [images {img2}];
            f1=i_cascadeexpr(ste.sce,T.lgene(j),ste.xy,k);
            title(sprintf('Ligand %d: %s',j,T.lgene(j)));
            saveas(f1,img1);
            f2=i_cascadeexpr(ste.sce,T.rgene(j),ste.xy,k+1);
            title(sprintf('Receptor %d: %s',j,T.rgene(j)));
            saveas(f2,img2);
            slide4 = add(ppt,'Comparison');
            replace(slide4,'Title',sprintf('Ligand %d - Receptor %d',j,j));
            replace(slide4,'Left Text',sprintf('Ligand %d: %s',j,T.lgene(j)));
            replace(slide4,'Right Text',sprintf('Receptor %d: %s',j,T.rgene(j)));
            replace(slide4,'Left Content',Picture(img1));
            replace(slide4,'Right Content',Picture(img2));
        end
        fw=gui.gui_waitbar;
        close(ppt);
        len = length(images);
        for i = 1:len
            delete(images{i});
        end
        gui.gui_waitbar(fw);
        rptview(ppt);            
    end
end


    function [f]=i_cascadeexpr(sce,g,xy,k)
        if nargin<4, k=1; end
        f=figure('Visible','off');
        tb = uitoolbar('Parent', f);
        P = get(f,'Position');
        set(f,'Position',[P(1)-120 P(2) round(0.75*P(3)) round(0.75*P(4))]);
        ix=find(sce.g==g);
        cx=sce.X(ix(1),:);
        dsize=15;
        h=scatter(xy(:,1),xy(:,2),dsize*0.75,cx,'filled');
        axis ij
        gui.i_setautumncolor(cx,'parula');
        title(g)
        cb=colorbar;
        cb.Label.String = 'UMI';
        pkg.i_addbutton2fig(tb,'off',{@i_RescaleExpr,h,cb},'IMG00067.GIF','Scale expression level using log2-transformation');
        pkg.i_addbutton2fig(tb,'off',{@i_ResetExpr,h,cx,cb},'IMG00074.GIF','Reset expression level to UMI');
        pkg.i_addbutton2fig(tb,'on',{@i_genecards,g},'fvtool_fdalinkbutton.gif','GeneCards...');
        P = get(f,'Position');
        set(f,'Position',[P(1)-20*k P(2)-20*k P(3) P(4)]);
        set(f,'visible','on');
        drawnow;
        box on
        grid on
%        dt = datacursormode;
%        dt.UpdateFcn = {@i_myupdatefcnx2,cx};
    end
