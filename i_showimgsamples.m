function i_showimgsamples(fnm)

if nargin < 1
    mfolder = fileparts(mfilename('fullpath'));
    fnm = fullfile(mfolder, 'example_data', 'sampled_img_resnet18.mat');
end
load(fnm, 's', 'idx', 'imarray');

figure;
scatter(s(:, 1), s(:, 2), [], idx);
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

    function [txt] = i_myupdatefcnx(~, event_obj)
        persistent myupdatefcn3fig
        if isempty(myupdatefcn3fig) || ~isvalid(myupdatefcn3fig)
            myupdatefcn3fig = figure('Position', [400, 400, 100, 100]);
        end
        if isvalid(myupdatefcn3fig) && isa(myupdatefcn3fig, 'matlab.ui.Figure')
            figure(myupdatefcn3fig);
        end
        imagesc(imarray{event_obj.DataIndex});
        axis tight
        axis equal
        txt = {num2str(event_obj.DataIndex)};
end

end