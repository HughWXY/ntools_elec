% function ntools_elec(varargin)
clear; close;

% check for spm, fsl and freesurfer
if isempty(which('spm')), error('Please install SPM and set up properly.\n'); end
if isempty(getenv('FSLDIR')),error('Please install FSL and set up properly.\n'); end
if isempty(getenv('FREESURFER_HOME')), error('Please instll Freesurfer and set up properly\n'); end

%% get the subject info
Ssdir=getenv('SUBJECTS_DIR');
disp('Select the subject FreeSurfer reconstruction folder');
Subject_path = uigetdir(Ssdir,'Select the subject FreeSurfer reconstruction folder');

if Subject_path==0
    disp('No Freesurfer folder was selected');
    return
end

ii = strfind(Subject_path,'/');
Ssdir = Subject_path(1:ii(end));
Sdir = Subject_path(ii(end)+1:end);

%% -------NYU settings: subject name
Sub = regexp(Sdir,'NY+[0-9a-zA-Z]+_','match');
if ~isempty(Sub)
    ss2 = Sdir(length(char(Sub))+1:end);
else
    ss2 = Sdir;
end

%% read the inital text file
disp('Select the initial text file');
[FileName,PathName] = uigetfile({'*.xlsx';'*.txt';'*.xls'},'Select the initial text file',pwd);
jj = strfind(PathName,'/');
Sname = PathName(jj(end-1)+1:jj(end)-1);

%% -------NYU settings: subject name
Sname = [Sname '_' ss2];
disp(Sname)

%% read the removed elecs text file
disp('Select the removed elec text file');
[removed_elec_file, removed_elec_path] = uigetfile({'*.xlsx';'*.txt';'*.xls'},...
    'Select the removed elec text file',PathName);

if isnumeric(removed_elec_file) || isempty(removed_elec_file)
    removed_elec = [];
else
    [~,~,ext] = fileparts(removed_elec_file);
    
    if strcmpi(ext,'.txt')
        fid = fopen(fullfile(removed_elec_path,removed_elec_file));
        removed_elec = textscan(fid,'%s','CommentStyle','%');
        fclose(fid);
        
    elseif strcmpi(ext,'.xls') || strcmpi(ext,'.xlsx')
        [~,~,removed_elec] = xlsread(fullfile(removed_elec_path,removed_elec_file));
    end
    
    % upper-case and zero-padding
    removed_elec = upper(removed_elec);
    removed_elec = regexprep(removed_elec,'(?<!\d)(\d)(?!\d)','0$1');
end
%%
disp('Select the T1 pre-operation image');
[preop_img_file,preop_img_path] = uigetfile({'*.nii.gz';'*.nii';'*.img'},...
    'Select the T1 pre-operation image',PathName);
    
% move other files to backup folder
backup_dir = [PathName,'backup_',datestr(now,29)];
if ~exist(backup_dir,'dir'), mkdir(backup_dir); end
all_files = dir(PathName);
all_files_isdir = [all_files(:).isdir];
all_files_name = {all_files(:).name};
all_files_name = all_files_name(all_files_isdir==0);
idx = ~cellfun(@isempty,strfind(all_files_name,'.log'))+...
    ~cellfun(@isempty,strfind(all_files_name,'pial_surf.mat'))+...
    ~cellfun(@isempty,strfind(all_files_name,'.ppt'))+...
    ~cellfun(@isempty,strfind(all_files_name,'.pptx'))+...
    ~cellfun(@isempty,strfind(all_files_name,FileName))+...
    ~cellfun(@isempty,strfind(all_files_name,preop_img_file))+...
    ~cellfun(@isempty,strfind(all_files_name,'aparc'))+... % keep aparc.annot file
    ~cellfun(@isempty,strfind(all_files_name,'removed'))+... % keep removed elec file
    ~cellfun(@isempty,strfind(all_files_name,'missing')); % keep the missing_coor.txt
all_files_name = all_files_name(~logical(idx));
for i=1:length(all_files_name), movefile([PathName, all_files_name{i}],backup_dir,'f'); end
    
% start diary 
diary_file = [PathName,'localization_process_',datestr(now,29),'.log'];
diary on;
diary(diary_file)
fprintf('\n================================================================\n');
fprintf('Starting localization process for %s at %s\n',Sname,datestr(now,31));
fprintf('Freesurfer Recon dir: %s\n',Subject_path);
fprintf('Initial location text file: %s\n',fullfile(PathName,FileName));

% read in preop T1
hdr = ntools_elec_load_nifti([preop_img_path preop_img_file],1);
if ~isequal(hdr.pixdim(2),hdr.pixdim(3),hdr.pixdim(4))
    warning('T1 voxel mm dimensions not equal. Will affect the accuracy of distance calculation');
end
scale = mean(hdr.pixdim(2:4));

%% read in initial text/xls file
[~,~,ext] = fileparts(FileName);

if strcmpi(ext,'.txt') 
    fid = fopen([PathName, FileName]);
    ini_elec_all = textscan(fid,'%s %f %f %f %s','CommentStyle','%'); 
    ini_elec_all = [ini_elec_all{1},num2cell(ini_elec_all{2}),...
                    num2cell(ini_elec_all{3}),num2cell(ini_elec_all{4}),ini_elec_all{5}];
    fclose(fid);
    
elseif strcmpi(ext,'.xls') || strcmpi(ext,'.xlsx')
    [~,~,ini_elec_all] = xlsread(fullfile(PathName,FileName));
end

%% transform voxel coordinates to ras coordinates if necessary
coor = cell2mat(ini_elec_all(:,2:4));
if all(floor(coor(:))==coor(:)) % if all coordinates are integer
    voxcoor = coor;
    [~,msg] = unix(sprintf('mri_info --vox2ras %s',fullfile(preop_img_path,preop_img_file)));
    vox2ras = str2num(msg);
    RAScoor = vox2ras*[voxcoor,ones(size(voxcoor,1),1)]';
    RAScoor = RAScoor(1:3,:)';
else
    RAScoor = coor;
end



%%
[status,msg] = unix(sprintf('mri_info --tkr2scanner %s',fullfile(preop_img_path,preop_img_file)));
if ~status
    transform = str2num(msg);
    scanner2tkr = -transform(1:3,4)'; % tkrRAS = scannerRAS + scanner2tkr
end

%% transform scanner RAS to tkr surface RAS
tkrRAS = round(RAScoor + repmat(scanner2tkr,size(ini_elec_all,1),1));
tkrRAS = num2cell(tkrRAS);
ini_elec_all_tkrRAS = ini_elec_all;
ini_elec_all_tkrRAS(:,2:4) = tkrRAS;

%% split elecs by type
if 5 == size(ini_elec_all_tkrRAS,2)
    g = strncmpi('G',ini_elec_all_tkrRAS(:,5),1);
    d = strncmpi('D',ini_elec_all_tkrRAS(:,5),1);
else
    g = strncmpi('G',ini_elec_all_tkrRAS(:,1),1);
    d = strncmpi('D',ini_elec_all_tkrRAS(:,1),1);
end

ini_grid_tkrRAS = ini_elec_all_tkrRAS(g,:);
ini_depth_tkrRAS = ini_elec_all_tkrRAS(d,:);
ini_elec_all_tkrRAS(or(g,d),:) = [];
ini_strip_tkrRAS = ini_elec_all_tkrRAS;


%% hemi
hemi = menu('Electrodes on which hemisphere?','left hemi','right hemi','both hemi');
if hemi==1
    sph_s = 'lh'; fprintf('\nElectrodes are on Left hemisphere\n');
elseif hemi==2
    sph_s = 'rh'; fprintf('\nElectrodes are on Right hemisphere\n');
else
    sph_s = 'both'; fprintf('\nElectrodes are on Both hemispheres\n');
end

%% Calculate the electrodes locations
% outer-brain surface check and create
ntools_elec_outer_brain(Subject_path)

% calculate grids
[elec_grid_tkrRAS,grid_stats] = ntools_elec_calc_grid(ini_grid_tkrRAS,Subject_path,scale,[]);
% remove the electrodes that were cut off during operation
if ~isempty(removed_elec)
    elec_grid_tkrRAS(ismember(elec_grid_tkrRAS(:,1),removed_elec),:) = []; 
end

%%
% calculate depth elecs
elec_depth_tkrRAS = ntools_elec_calc_depth(ini_depth_tkrRAS);

% calculate strip elecs
elec_strip_tkrRAS = ntools_elec_calc_strip(ini_strip_tkrRAS,Subject_path,sph_s);

% save the electrodes locations into a text file
fname_t1 = [PathName Sname '_coor_T1_' datestr(now,29) '.txt'];
ntools_elec_savetxt(fname_t1,elec_grid_tkrRAS,elec_strip_tkrRAS,elec_depth_tkrRAS);
fid = fopen(fname_t1);
text = textscan(fid,'%s %f %f %f %s');
fclose(fid);
name = text{1}; x = text{2}; y = text{3}; z = text{4}; label = text{5};

% save all into binary nifti image
fname_bin = [PathName,Sname,'_elec_bin_T1_' datestr(now,29), '.nii.gz'];
elec_all_tkrRAS = [x y z];
elec_all_scannerRAS = elec_all_tkrRAS - repmat(scanner2tkr,size(elec_all_tkrRAS,1),1);
elec_vox = ntools_elec_savebin(elec_all_scannerRAS,hdr,fname_bin);

%% transform into mni space
elec_mni = ntools_elec_dartel_warp(fname_bin,fullfile(preop_img_path,preop_img_file));
fname_mni = [PathName Sname '_coor_MNI_' datestr(now,29) '.txt'];
ntools_elec_savetxt(fname_mni,[name num2cell(elec_mni) label]);


%% Save the surf.mat and plot

if strcmp(sph_s,'both')
    surf_brain_lh = fs_load_subj(Sdir,'lh','pial',0,Ssdir);
    if ~isfield(surf_brain_lh,'coords')
        surf_brain_lh.coords = surf_brain_lh.vertices;
    end
    surf_brain_rh = fs_load_subj(Sdir,'rh','pial',0,Ssdir);
    if ~isfield(surf_brain_rh,'coords')
        surf_brain_rh.coords = surf_brain_rh.vertices;
    end    
    
    lh_mat = [PathName,Sdir,'_lh_pial_surf.mat'];
    rh_mat = [PathName,Sdir,'_rh_pial_surf.mat'];

    save(lh_mat,'-struct','surf_brain_lh','coords','faces');
    save(rh_mat,'-struct','surf_brain_rh','coords','faces');

    % split the text file into lh/rh hemi, then save with
    % NYU_ntools_elec_autoplot
    
    aparc = {[Subject_path,'/label/lh.aparc.split_STG_MTG.annot'];...
             [Subject_path,'/label/rh.aparc.split_STG_MTG.annot']}; 
    aparc2 = {[Subject_path,'/label/lh.aparc.a2009s.annot'];...
              [Subject_path,'/label/rh.aparc.a2009s.annot']};
    if ~exist(aparc{1},'file') || ~exist(aparc{2},'file')
        fs_mris_divide_parcellation(Sdir)
    end
        
    copyfile(aparc{1},PathName,'f');
    copyfile(aparc{2},PathName,'f');
    copyfile(aparc2{1},PathName,'f');
    copyfile(aparc2{2},PathName,'f');    
    
    clear surf_brain_lh surf_brain_rh 
else
    surf_brain = fs_load_subj(Sdir,sph_s,'pial',0,Ssdir);
    if ~isfield(surf_brain,'coords')
        surf_brain.coords = surf_brain.vertices;
    end    
    surf_mat = [PathName,Sdir,'_',sph_s,'_pial_surf.mat'];

    save(surf_mat,'-struct','surf_brain','coords','faces');
    
    % NYU settings: auto save images with/without aparc
    aparc = [Subject_path,'/label/',sph_s,'.aparc.split_STG_MTG.annot']; 
    if ~exist(aparc,'file'), fs_mris_divide_parcellation(Sdir); end
    aparc2 = [Subject_path,'/label/',sph_s,'.aparc.a2009s.annot'];
    copyfile(aparc,PathName,'f'); 
    copyfile(aparc2,PathName,'f');
%     NYU_ntools_elec_autoplot(fname_t1,surf_mat,aparc);
%     NYU_ntools_elec_autoplot(fname_t1,surf_mat,aparc2);
    
    clear surf_brain
end

% close diary
fprintf('\nElectrodes Localization finished for %s',Sname);
fprintf('\n================================================================\n');
diary off

% somehow an empty "diary" file was left after diary off
dfile = fullfile(PathName,'diary');
if exist(dfile,'file'), delete(dfile); end

%% -------NYU settings: save to a mat file
% matfile = ['/home/halgdev/projects/nyuproj/loc/NY_struct/',Sname,'_',datestr(now,'mmddyy')];
% save(matfile);

% end

