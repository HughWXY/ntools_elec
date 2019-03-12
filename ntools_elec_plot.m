function ntools_elec_plot(varargin)

% a stand-alone program that shows ieeg electrodes on the pial surface and
% save the images into textfile folder/images/. Default saving name is
% PatientID(NYxx)_space_elecLabel_viewpoint_hemisphere.png
% 
% space could be T1 or MNI
% elecLabel could be grid, strip, depth, or grid and strip
% viewpoint could be left, right, top, below, back or front
% hemisphere could be lh, rh, or both
% 
% required input:
% elec_text: text file with xyz electrode coordinates
% pial_mat: matlab structure with pial surface
%
% optional input:
% plot: part to plot 
%       'G': Grid only 
%       'S': Strip only 
%       'D': Depth Only 
%       'GS': Both Grid and Strip  
%       'brain': brain image only
% showlabel: to show the electrodes' labels (1) or not (0)
% saveimg: to save the images (1) or not (0),if set to (10), save in silence
% aparc: plot the image with aparc.annot, all electrodes and labels will be white
%
% Usage: run ntools_elec_plot in command window
% the gui is prompted out for file selection
% or: ntools_elec_plot(fname, {'lh.mat', 'rh.mat'})
%     ntools_elec_plot(fname, l(r)h.mat);
%     ntools_elec_plot(fname, l(r)h.mat,'showlabel',1,'saveimg',10);
%     ntools_elec_plot(fname,lh.mat,'plot','G','aparc','lh.aparc.annot');
%
% also see: http://ieeg.pbworks.com/Viewing-Electrode-Locations
%
% written by  Hugh Wang, Xiuyuan.Wang@nyumc.org, May 13, 2009
%
% modified on May 14, 2009 by Hugh
% make judgement on the input file type and not sensitive to the order of 
% input variable. add the option to show the electrodes' labels or not.
%
% modified on July 22nd, 2009 by Hugh
% for subjects who has electrodes on both hemisphere, loading the both
% pial.mat file will generate the image with whole brain and save the
% images from 6 views (left, right, top, below, back & front)
%
% modified on Aug 8th, 2009 by Hugh
% show only lh(rh)'s electrodes if choose lh(rh)_surf.mat
%
% modified on Jan 28th, 2010 by Hugh
% default saving name
%
% modified on Jan 13th, 2014 by Hugh
% more optional inputs 'plot' and 'aparc'
% 
% modified on Feb 11th, 2019 by Hugh
% support for experimental grid


%% Get the elec info
if nargin==0
    [FileName,PathName] = uigetfile('*.txt','Select the electrodes text file',pwd'); 
    [surfname, surfpath] = uigetfile('*.mat','Select the patient brain surf',PathName,'MultiSelect','on');
    surf = strcat(surfpath,surfname);      
elseif nargin>=2
    aa = strfind(varargin{1},'/');
    if ~isempty(aa)
        FileName = varargin{1}(aa(end)+1:end);
        PathName = varargin{1}(1:aa(end));
    else
        FileName = varargin{1};
        PathName = [pwd,'/'];
    end
    surf = varargin{2}; 

    try labelshow = varargin{find(strcmp('showlabel',varargin))+1}; catch err; end
    try genimg = varargin{find(strcmp('saveimg',varargin))+1}; catch err; end
    try plt = varargin{find(strcmp('plot',varargin))+1}; catch err; end
    try aparc = varargin{find(strcmp('aparc',varargin))+1}; catch err; end

end

if exist(fullfile(PathName, FileName),'file')
    fid = fopen(fullfile(PathName, FileName));
    elec_all = textscan(fid,'%s %f %f %f %s');
    fclose(fid);
    
    elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];
else
%     elec_cell = [];
    disp('No electrode was found. Please check you input text file.')
    return
end

%% Get the filename info
b = strfind(FileName,'_');
Pname = FileName(1:b(1)-1);

if contains(upper(FileName),'T1')
    space = '_T1_';
elseif contains(upper(FileName),'MNI')
    space = '_MNI_';
else
    space = '_';
end

if length(surf)==2
    sph = 'both';
else
    sph = regexpi(surf,'[r,l]h','match');
    sph = char(sph{:});
end

%% Separate Grid, Strip and Depth electrodes

if isempty(char(elec_all{5}(:)))
    g = strncmpi('G',elec_cell(:,1),1);
    d = strncmpi('D',elec_cell(:,1),1);
else
    g = strncmpi('G',elec_all{5},1);
    d = strncmpi('D',elec_all{5},1);
    eg = strncmpi('EG',elec_all{5},2);
end

elec_grid = elec_cell(g,:);
elec_expGrid = elec_cell(eg,:);
elec_depth = elec_cell(d,:);
elec_cell(logical(g+d+eg),:) = [];


%% Plot the elecs
% what to plot
if ~exist('plt','var')
    plt = menu('What part do you want to plot?','Grid only', 'Strip only','Depth Only','Both Grid and Strip','Brain only');
end
% show label
if ~exist('labelshow','var')
    labelshow = menu('Do you want to show the label?','Yes','No');
end

% plot with aparc
if ~exist('aparc','var')
    aparc = menu('Do you want to plot with aseg parcellations? (note: all electrodes and labels will be white)','Yes','No');
    if aparc==1 
%         aparc = uigetdir(getenv('SUBJECTS_DIR'),'Select subject FS recon folder');
        [aparc,aparc_path]= uigetfile('*.annot','Select the aparc annotation file',PathName,'MultiSelect','on');
        aparc = fullfile(aparc_path,aparc);
    else
        aparc = [];
    end
elseif isempty(aparc) || isnumeric(aparc)
    aparc = [];
end

% save image
if ~exist('genimg','var')
    genimg = menu('Do you want to save the images?','Yes', 'No');
end

% load pial surfs and corresponding aparc annotation file
if strcmp(sph,'both')
    surf_brain.sph1 = load(surf{1});
    surf_brain.sph2 = load(surf{2});

    if ~isempty(aparc)
        hemi1 = regexpi(surf{1},'[r,l]h','match'); hemi1 = char(hemi1{:});
%         aparc_annot.hemi1 = char(aparc(strncmpi(hemi1,aparc,2)));
        idx = strfind(aparc,hemi1);
        idx(cellfun(@isempty,idx)) = {0};
        aparc_annot.hemi1 = char(aparc(logical(cell2mat(idx))));
        
        hemi2 = regexpi(surf{2},'[r,l]h','match'); hemi2 = char(hemi2{:});
%         aparc_annot.hemi2 = char(aparc(strncmpi(hemi2,aparc,2)));
        idx = strfind(aparc,hemi2);
        idx(cellfun(@isempty,idx)) = {0};
        aparc_annot.hemi2 = char(aparc(logical(cell2mat(idx))));
    else
        aparc_annot = [];
    end
    
else 
    surf_brain = load(surf);
%     if ~isempty(aparc), aparc_annot = strcat(aparc,'/label/',sph,'.aparc.annot'); else aparc_annot = []; end
    aparc_annot = aparc;
end

% main plot
if (isequal(plt,1) || strcmpi(plt,'G')) && ~isempty(elec_grid)
    showpart = 'G';
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),'r',labelshow,aparc_annot);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% experimental Grid
    if ~isempty(elec_expGrid)
        hold on;
        expGrid_radius = 0.5;
        
        for i=1:size(elec_expGrid,1)
            if isempty(aparc_annot)
                expGrid_color = 'm';
            else % re-define strip color if plot with aparc
                expGrid_color = 'w';
            end            
            
            plotSpheres(elec_expGrid{i,2},elec_expGrid{i,3},elec_expGrid{i,4},expGrid_radius,expGrid_color);
            if labelshow==1
                [xx, yy, zz] = adjust_elec_label([elec_expGrid{i,2:4}],expGrid_radius); % default radius = 2
                text('Position',[xx yy zz],'String',elec_expGrid(i,1),'Color','w','VerticalAlignment','top');
            end    
        end
        hold off;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (isequal(plt,2) || strcmpi(plt,'S')) && ~isempty(elec_cell)
    showpart = 'S';
    nyu_plot(surf_brain,sph,cell2mat(elec_cell(:,2:4)),char(elec_cell(:,1)),'b',labelshow,aparc_annot);
elseif (isequal(plt,3) || strcmpi(plt,'D')) && ~isempty(elec_depth)
    showpart = 'D';
    nyu_plot(surf_brain,sph,cell2mat(elec_depth(:,2:4)),char(elec_depth(:,1)),'g',labelshow,[],1.5,0.3);
elseif (isequal(plt,4) || strcmpi(plt,'GS')) && ~isempty(elec_grid) && ~isempty(elec_cell)
    showpart = 'GS';
    elec = cell2mat(elec_cell(:,2:4));
    elec_name = char(elec_cell(:,1));
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),'r',labelshow,aparc_annot); hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% experimental Grid
    if ~isempty(elec_expGrid)
        hold on;
        expGrid_radius = 0.5;
        
        for i=1:size(elec_expGrid,1)
            if isempty(aparc_annot)
                expGrid_color = 'm';
            else % re-define strip color if plot with aparc
                expGrid_color = 'w';
            end
            
            plotSpheres(elec_expGrid{i,2},elec_expGrid{i,3},elec_expGrid{i,4},expGrid_radius,expGrid_color);
            if labelshow==1
                [xx, yy, zz] = adjust_elec_label([elec_expGrid{i,2:4}],expGrid_radius); % default radius = 2
                text('Position',[xx yy zz],'String',elec_expGrid(i,1),'Color','w','VerticalAlignment','top');
            end    
        end
        hold off;        
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:size(elec,1)
        if isempty(aparc_annot)
            stripcolor = 'b';
        else % re-define strip color if plot with aparc
            stripcolor = 'w';
        end
        plotSpheres(elec(i,1),elec(i,2),elec(i,3),2,stripcolor);
        if labelshow==1
            [xx, yy, zz] = adjust_elec_label(elec(i,:)); % default radius = 2
            text('Position',[xx yy zz],'String',elec_name(i,:),'Color','w','VerticalAlignment','top');
        end
    end
    hold off;  
elseif isequal(plt,5) || strcmpi(plt,'brain')
    showpart = 'Brain';
    nyu_plot(surf_brain,[],[],[]);
else
    disp('sorry, the electrodes you choose to plot are not on the surface you loaded');
    return;
end

%% save images

if genimg==1 || genimg==10
    if ~exist([PathName 'images/'],'dir')
        mkdir([PathName 'images/']);
    end
    
    if labelshow==1
        label = '_label';
    else
        label = [];
    end
    
    if ~isempty(aparc)
        % split aparc file name and get the middle one in output image file
        % name
        aparc_fnames = regexp(aparc,'\.','split');
        aparc = ['_',aparc_fnames{3}];
    end
    
    format = 'png';
    outimgfname = [PathName,'images/',Pname,space,showpart,'_XaxisviewsX_',sph,label,aparc];
    
    % save the images in silence
    if genimg==10, set(gcf,'Visible','off'); end
    
    if strcmp(sph,'lh')
        view(270, 0);
        saveas(gcf,regexprep(outimgfname,'XaxisviewsX','lateral'),format);
        view(90,0);
        saveas(gcf,regexprep(outimgfname,'XaxisviewsX','mesial'),format);
        
    elseif strcmp(sph,'rh')
        view(270, 0);
        saveas(gcf,regexprep(outimgfname,'XaxisviewsX','mesial'),format);
        view(90,0);
        saveas(gcf,regexprep(outimgfname,'XaxisviewsX','lateral'),format);
        
    elseif strcmp(sph,'both')
        view(270, 0);
        saveas(gcf,regexprep(outimgfname,'XaxisviewsX','left'),format);
        view(90,0);
        saveas(gcf,regexprep(outimgfname,'XaxisviewsX','right'),format);
    end
    view(0,0);
    saveas(gcf,regexprep(outimgfname,'XaxisviewsX','posterior'),format);

    view(180,0);
    saveas(gcf,regexprep(outimgfname,'XaxisviewsX','frontal'),format);

    view(90,90);
    saveas(gcf,regexprep(outimgfname,'XaxisviewsX','dorsal'),format);

    view(90,-90);
    set(light,'Position',[1 0 -1]);
    saveas(gcf,regexprep(outimgfname,'XaxisviewsX','ventral'),format);
else 
    return;
end

end

%% subfunctions 
%% nyu_plot
function nyu_plot(surf_brain,sph,elec,elecname,color,label,aparc_annot,radius,alpha)

if ~exist('color','var')
    color = 'w';
end
if ~exist('label','var')
    label = 0;
end
if ~exist('alpha','var')
    alpha = 1;
end
if ~exist('radius','var')
    radius = 2;
end
if ~exist('aparc_annot','var')
    aparc_annot = [];
end

figure;

col = [.7 .7 .7];

if strcmp(sph,'both')
    sub_sph1.vert = surf_brain.sph1.coords;
    sub_sph1.tri = surf_brain.sph1.faces;

    sub_sph2.vert = surf_brain.sph2.coords;
    sub_sph2.tri = surf_brain.sph2.faces;
    
    if isempty(aparc_annot) || ~exist(aparc_annot.hemi1,'file') || ~exist(aparc_annot.hemi2,'file')
        col1=repmat(col(:)', [size(sub_sph1.vert, 1) 1]);
        col2=repmat(col(:)', [size(sub_sph2.vert, 1) 1]);
    else
        [~,albl1,actbl1]=fs_read_annotation(aparc_annot.hemi1);
        [~,aa] = ismember(albl1,actbl1.table(:,5));
        aa(aa==0) = 1;
        col1 = actbl1.table(aa,1:3)./255;
        
        [~,albl2,actbl2]=fs_read_annotation(aparc_annot.hemi2);
        [~,bb] = ismember(albl2,actbl2.table(:,5));
        bb(bb==0) = 1;
        col2 = actbl2.table(bb,1:3)./255;
                
        % re-define the electrode color if plot with aparc
        color = 'w';
    end    
    
    trisurf(sub_sph1.tri, sub_sph1.vert(:, 1), sub_sph1.vert(:, 2),sub_sph1.vert(:, 3),...
        'FaceVertexCData', col1,'FaceColor', 'interp','FaceAlpha',alpha);
    hold on;
    trisurf(sub_sph2.tri, sub_sph2.vert(:, 1), sub_sph2.vert(:, 2), sub_sph2.vert(:, 3),...
        'FaceVertexCData', col2,'FaceColor', 'interp','FaceAlpha',alpha);
else    
    if isfield(surf_brain,'coords')==0
        sub.vert = surf_brain.surf_brain.coords;
        sub.tri = surf_brain.surf_brain.faces;
    else
        sub.vert = surf_brain.coords;
        sub.tri = surf_brain.faces;
    end
    
    if isempty(aparc_annot) || ~exist(aparc_annot,'file')
        col=repmat(col(:)', [size(sub.vert, 1) 1]);
    else
        [~,albl,actbl]=fs_read_annotation(aparc_annot);
        [~,cc] = ismember(albl,actbl.table(:,5));
        cc(cc==0) = 1;
        col = actbl.table(cc,1:3)./255;
        
        % re-define the electrode color if plot with aparc
        color = 'w';
    end    
        
    trisurf(sub.tri, sub.vert(:, 1), sub.vert(:, 2), sub.vert(:, 3),...
        'FaceVertexCData', col,'FaceColor', 'interp','FaceAlpha',alpha);
end

shading interp;
lighting gouraud;
material dull;
light;
axis off
hold on;
for i=1:size(elec,1)
    plotSpheres(elec(i,1),elec(i,2),elec(i,3),radius,color);
    if label==1
        [x, y, z] = adjust_elec_label(elec(i,:),radius);
        text('Position',[x y z],'String',elecname(i,:),'Color','w','VerticalAlignment','top');
    end
end
set(light,'Position',[-1 0 1]); 
    if strcmp(sph,'lh')
        view(270, 0);      
    elseif strcmp(sph,'rh')
        view(90,0);        
    elseif strcmp(sph,'both')
        view(90,90);
    end
set(gcf, 'color','black','InvertHardCopy', 'off');
axis tight;
axis equal;
end

%% adjust_elec_label
function [x, y, z] = adjust_elec_label(elec,radius)

if ~exist('radius','var')
    radius = 2;
end

if elec(1)>0
    x = elec(1)+radius;
else
    x = elec(1)-radius;
end

if elec(3)>0
    z = elec(3)+radius;
else
    z = elec(3)-radius;
end

y = elec(2);

end

%% plotSpheres
function [shand]=plotSpheres(spheresX, spheresY, spheresZ, spheresRadius,varargin)

if nargin>4
    col=varargin{:};
end

spheresRadius = ones(length(spheresX),1).*spheresRadius;
% set up unit sphere information
numSphereFaces = 25;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);

% set up basic plot
sphereCount = length(spheresRadius);

% for each given sphere, shift the scaled unit sphere by the
% location of the sphere and plot
for i=1:sphereCount
sphereX = spheresX(i) + unitSphereX*spheresRadius(i);
sphereY = spheresY(i) + unitSphereY*spheresRadius(i);
sphereZ = spheresZ(i) + unitSphereZ*spheresRadius(i);
shand=surface(sphereX, sphereY, sphereZ,'FaceColor',col,'EdgeColor','none','AmbientStrength',0.7);
end

end
