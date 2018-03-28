% function ntools_elec_autocoreg
clear; clc
% using flirt dof 6 to coregister the postop MRI to preop MRI (usually FS
% recon T1)

coreg = '/home/wangx11/Synology/Loc/'; % your localization folder
subj = getenv('SUBJECTS_DIR');

t1 = menu('Select the pre-operation image','Pick up my own T1','Select Freesurfer T1.mgz');
if t1==1
    [preop, preoppath] = uigetfile('*.*', 'Select the pre-operation image',coreg);
else
    [preop, preoppath] = uigetfile('*.*', 'Select the pre-operation image',subj);
end
preop = fullfile(preoppath,preop);

[aseg,asegpath] = uigetfile('*.mgz','Select subject aseg file',preoppath);
if ~isnumeric(aseg), aseg = fullfile(asegpath,aseg); else aseg = []; end;

[elecT1, elecpath] = uigetfile('*.*', 'Select the post-operation T1 image with electrodes',coreg);
if ~isnumeric(elecT1),elecT1 = fullfile(elecpath,elecT1); else disp('No electrode image selected'); return; end

[elecT2, elecpathT2] = uigetfile('*.*', 'Select the post-operation T2 image with electrodes',elecpath);
if ~isnumeric(elecT2),elecT2 = fullfile(elecpathT2,elecT2); else elecT2 = []; end

% convert format and orientation
% preop
[~,name_preop,ext] = fileparts(preop);
[~,name_preop] = fileparts(name_preop);
preop_nii = fullfile(elecpath,[name_preop,ext]);
preop_nii = regexprep(preop_nii,ext,'.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',preop,preop_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

% aseg
if ~isempty(aseg)
    [~,name_aseg,ext] = fileparts(aseg);
    [~,name_aseg] = fileparts(name_aseg);
    aseg_nii = fullfile(elecpath,[name_aseg,ext]);
    aseg_nii = regexprep(aseg_nii,ext,'.nii.gz');
    cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',aseg,aseg_nii);
    [status, msg] = unix(cmd_convert);
    if status, disp(msg); return; end
end

% elecT1
elecT1_nii = regexprep(elecT1,'.img','.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',elecT1,elecT1_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

if ~isempty(elecT2)
    elecT2_nii = regexprep(elecT2,'.img','.nii.gz');
    cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',elecT2,elecT2_nii);
    [status, msg] = unix(cmd_convert);
    if status, disp(msg); return; end 
end

% coregister
[~,name_elecT1] = fileparts(elecT1_nii);
[~,name_elecT1] = fileparts(name_elecT1);
elec_preopT1 = fullfile(elecpath,[name_elecT1,'_preop.nii.gz']);
% elec_preop_brainT1 =  [elecpath,name_elecT1,'_preop_brain.nii.gz'];
% elec_preop_omat = [elecpath,name_elecT1,'_preop.mat'];
unix(sprintf('flirt -in %s -ref %s -out %s -dof 6 -interp trilinear',elecT1_nii,preop_nii,elec_preopT1),'-echo');

preop_brain = fullfile(elecpath,[name_preop,'_brain.nii.gz']);
preop_brain_mask = fullfile(elecpath,[name_preop,'_brain_mask.nii.gz']);
unix(sprintf('bet %s %s -f 0.5 -g 0 -m',preop_nii,preop_brain),'-echo');

if ~isempty(elecT2)
    [~,name_elecT2] = fileparts(elecT2_nii);
    [~,name_elecT2] = fileparts(name_elecT2);
    elec_preopT2 = fullfile(elecpath,[name_elecT2,'_preop.nii.gz']);
    unix(sprintf('flirt -in %s -ref %s -out %s -dof 6 -interp trilinear',elecT2_nii,preop_nii,elec_preopT2),'-echo');
end

if ~isempty(aseg)
    % remove cerebellum
    cerebellum = [elecpath,'cerebellum.nii.gz'];
    unix(sprintf('mri_binarize --i %s --match 6 7 8 45 46 47 170 171 172 173 174 175--o %s',aseg_nii,cerebellum));
    % dialte cerebellum
    unix(sprintf('fslmaths %s -dilM %s',cerebellum,cerebellum));
    % inverse mask
    cerebellum_inv = [elecpath,'cerebellum_inv.nii.gz'];
    unix(sprintf('fslmaths %s -sub 1 -abs %s',cerebellum,cerebellum_inv));
%     elec_preop_cortex = regexprep(elec_preop_brainT1,'brain','cortex');
    elec_preop_cortexT1 = fullfile(elecpath,[name_elecT1,'_preop_cortex.nii.gz']);
    unix(sprintf('fslmaths %s -mas %s -mas %s %s',elec_preopT1,preop_brain_mask,cerebellum_inv,elec_preop_cortexT1));
    
    if ~isempty(elecT2)
        elec_preop_cortexT2 = fullfile(elecpath,[name_elecT2,'_preop_cortex.nii.gz']);
        unix(sprintf('fslmaths %s -mas %s -mas %s %s',elec_preopT2,preop_brain_mask,cerebellum_inv,elec_preop_cortexT2));
    end
    
    % check results
    if isempty(elecT2)
        unix(sprintf('fsleyes %s %s',elec_preopT1, elec_preop_cortexT1));
    else
        unix(sprintf('fsleyes %s %s %s',elec_preopT1, elec_preop_cortexT1, elec_preop_cortexT2));
    end
else
    % check results
    elec_preop_brainT1 =  fullfile(elecpath,[name_elecT1,'_preop_brain.nii.gz']);
    unix(sprintf('fslmaths %s -mas %s %s',elec_preopT1,preop_brain_mask,elec_preop_brainT1));
    unix(sprintf('fsleyes %s %s',elec_preopT1, elec_preop_brainT1));
end

