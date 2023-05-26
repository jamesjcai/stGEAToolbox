function [idx]=cluster_resnet18(img,xy,numclass,sposz,plotit)

%REF: https://doi.org/10.1101/2022.04.14.488259

if nargin<4, plotit=false; end
if nargin<3, numclass=4; end
if nargin<4, sposz=40; end

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


% ref: https://www.mathworks.com/matlabcentral/fileexchange/75342-image-clustering-and-dimension-reduction-using-cnn
% net = resnet18();
% net =resnet18('Weights','none');
[imarray]=st_sampleimg(img,xy,sposz);
fea=zeros(length(imarray),1000);
parfor kk=1:length(imarray)            
    fea(kk,:)=squeeze(activations(net,imresize(imarray{kk},[224 224]),'fc1000'));
end
idx=kmeans(fea,numclass,"Start","plus");

if plotit
    figure;
    s=tsne(fea);
    scatter(s(:,1),s(:,2),[],idx);
    dt = datacursormode;
    dt.UpdateFcn = {@i_myupdatefcnx};
    fnm=tempname;
    save(fnm,'s','idx','imarray');
    fprintf('f=''%s.mat'';\ni_showimgsamples(f);\n',fnm);
end

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
