function [idx]=st_cluster_img(ste,varargin)

p = inputParser;
defaultType = 'resnet18kmeans';
validTypes = {'imsegkmeans','resnet18kmeans','resnet18labeling'};
checkType = @(x) any(validatestring(x,validTypes));
checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);
addRequired(p,'ste',@(x) isa(x,'SpatialTranscriptomicsExperiment'));
addOptional(p,'method',defaultType,checkType);
addOptional(p,'plotit',false,@islogical);
addOptional(p,'sposz',40,checkK);
addOptional(p,'k',5,checkK);
parse(p,ste,varargin{:})

method=p.Results.method;
plotit=p.Results.plotit;
sposz=p.Results.sposz;
numClass=p.Results.k;

% if nargin<4, methodid=2; end
% if nargin<3 || isempty(sposz), sposz=40; end
% if nargin<2 || isempty(k), k=5; end

if contains(method,'resnet18')
try
    net = resnet18();
catch ME
    % disp(ME.message)
    answer=questdlg('This function requires the Deep Learning Toolbox Model for ResNet-18 Network support package for the pretrained weights. Install this support package?');
    if strcmp(answer,'Yes')
        matlab.addons.supportpackage.internal.explorer.showSupportPackages('RESNET18','tripwire');
    end
    rethrow(ME);
end
end
switch method
    case 'imsegkmeans'
        L = imsegkmeans(ste.img,numClass);
        l=sub2ind(size(ste.img,[1 2]), ...
            round(ste.xy(:,1)),round(ste.xy(:,2)));
        idx=double(L(l));
        if plotit
            figure;
            imshow(labeloverlay(ste.img,L));
            a=lines(length(unique(idx)));
            colormap(a);
            title(sprintf('k=%d',length(unique(idx))));
            colorbar;
        end
    case 'resnet18kmeans'
        % ref: https://www.mathworks.com/matlabcentral/fileexchange/75342-image-clustering-and-dimension-reduction-using-cnn
        % net = resnet18();
        % net =resnet18('Weights','none');
        [imarray]=st_sampleimg(ste.img,ste.xy,sposz);
        fea=zeros(length(imarray),1000);
        for kk=1:length(imarray)            
            fea(kk,:)=squeeze(activations(net,imresize(imarray{kk},[224 224]),'fc1000'));
        end
        idx=kmeans(fea,numClass,"Start","plus");

        if plotit
            figure;
%             [~,s]=pca(fea);
%             subplot(1,2,1)
%             gscatter(s(:,1),s(:,2),idx);
%             subplot(1,2,2)
%             % perform t-sne for the dimension reduction
            T=tsne(fea);
            scatter(T(:,1),T(:,2),[],idx);
            % gscatter(T(:,1),T(:,2),idx);
            dt = datacursormode;
            dt.UpdateFcn = {@i_myupdatefcnx};
        end
        return;
    case 'resnet18labeling'
        % net = resnet18();
        [imarray]=st_sampleimg(ste.img,ste.xy,sposz);
        idx=strings(length(imarray),1);
        for kk=1:length(imarray)
            idx(kk)=string(classify(net,imresize(imarray{kk},[224,224])));
        end
        if plotit
            a=randperm(length(imarray));
            figure;
            for kk=1:25
                subplot(5,5,kk)
                imshow(imarray{a(kk)});
                title(idx(a(kk)));
            end
        end
end
idx=grp2idx(idx);

    function [txt]=i_myupdatefcnx(~, event_obj)
        persistent myupdatefcn3fig
        if isempty(myupdatefcn3fig) || ~isvalid(myupdatefcn3fig)
            myupdatefcn3fig=figure('Position',[400 400 100 100]);
        end
        if isvalid(myupdatefcn3fig) && isa(myupdatefcn3fig,'matlab.ui.Figure')
           figure(myupdatefcn3fig);
        end
        imagesc(imarray{event_obj.DataIndex});
        axis tight
        axis equal
        txt = {num2str(event_obj.DataIndex)};
    end
end



% if methodid
%     figure;
%     B = labeloverlay(ste.img,L);
%     imshow(B);
% end
% idx=size(ste.xy,1);
% A=ste.img(:,:,1);
% for kk=1:length(idx)    
%     x1=round(ste.xy(kk,:)-sposz/2);
%     x2=round(ste.xy(kk,:)+sposz/2);
%     row_start=x1(1);
%     row_end=x2(1);
%     col_start=x1(2);
%     col_end=x2(2);
%     l = stack (A,row_start,row_end,col_start,col_end);
%     idx(kk)=mode(L(l));

% function result = stack (A,row_start,row_end,col_start,col_end)
% % A = [4     4     4     4     4     4     4
% %      4     1     1     1     1     3     0
% %      4     1     3     3     1     3     0
% %      4     1     3     3     1     3     0
% %      4     1     1     1     1     3     0
% %      4     4     4     4     4     4     4];
% %  row_start=3; col_start=4;
% %  row_end=5; col_end=6;
% height=(size(A,1));
% result=(row_start:row_end)+(height)*((col_start:col_end)'-1);
% result=transpose(result); result=result(:);
% end