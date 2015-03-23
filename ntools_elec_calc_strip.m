function [elec, data] = ntools_elec_calc_strip(ini_cell,subjectpath,sph)

% find the strip electrodes on the outer brain surface using nearest points

if isempty(ini_cell)
    disp('no strip electrodes detected');
    elec = []; data = [];
    return;
end
fprintf('Calculating the strip electrodes....'); tic;
strip = cell2mat(ini_cell(:,2:4));

if strcmp(sph,'both')
    surf_lh = fs_read_surf([subjectpath '/surf/lh.pial-outer-smoothed']);
    if ~isfield(surf_lh,'coords')
        surf_lh.coords = surf_lh.vertices;
    end
    surf_rh = fs_read_surf([subjectpath '/surf/rh.pial-outer-smoothed']);
    if ~isfield(surf_rh,'coords')
        surf_rh.coords = surf_rh.vertices;
    end    
    surf = [surf_lh.coords;surf_rh.coords];
else
    surf_h = fs_read_surf([subjectpath '/surf/',sph,'.pial-outer-smoothed']);
    if ~isfield(surf_h,'coords')
        surf_h.coords = surf_h.vertices;
    end    
    surf = surf_h.coords;
end
k = dsearchn(surf,strip);
data = surf(k,:);

% data = ICP_finite(surf,strip,struct('Optimizer','fminsearch'));
ini_cell(:,1) = regexprep(ini_cell(:,1),'(?<!\d)(\d)(?!\d)','0$1');
elec = [ini_cell(:,1), num2cell(data), repmat({'S'},[length(strip),1])];

fprintf('Done. (%f seconds) \n\n', toc);