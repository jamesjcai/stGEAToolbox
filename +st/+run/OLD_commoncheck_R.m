function [ok, msg] = OLD_commoncheck_R(rscriptdir)

ok = false;
msg = [];

if isempty(pkg.FindRpath)
    msg = ('Rscript.exe is not found');
    return;
end

folder = fileparts(mfilename('fullpath'));
wrkpth = fullfile(folder, rscriptdir);
cd(wrkpth);
fprintf('CURRENTWDIR = "%s"\n', wrkpth);
[~, cmdout] = pkg.RunRcode('require.R');
if strfind(cmdout, 'there is no package') > 0
    msg = cmdout;
    return;
end
ok = true;
end
