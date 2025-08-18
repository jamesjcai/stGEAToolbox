function removeMenu(figHandle, menuLabel)
% removeMenu Remove a top-level or nested menu from a figure by label
%
%   removeMenu(figHandle, menuLabel)
%
%   figHandle : handle to a figure
%   menuLabel : char or string, label of the menu item to remove
%
%   Example:
%       f = figure;
%       removeMenu(f,'Window');
%       removeMenu(f,'Desktop');
%       removeMenu(f,'New Figure');

    if nargin < 2
        error('Usage: removeMenu(figHandle, menuLabel)');
    end

    % find all menus in this figure
    m = findall(figHandle,'Type','uimenu');
    if isempty(m)
        warning('No menus found in figure.');
        return;
    end

    % match label (case-insensitive, exact match)
    labels = get(m, 'Label');
    if ischar(labels), labels = {labels}; end
    idx = find(strcmpi(labels, menuLabel));

    if isempty(idx)
        warning('Menu "%s" not found.', menuLabel);
    else
        delete(m(idx));
    end
end
