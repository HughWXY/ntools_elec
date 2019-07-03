function elec_mni = ntools_elec_dartel_warp(elec_vox,preop_t1)

fprintf('start DARTEL warping process......\n')

tic
[preop_t1_path,preop_t1_img,ext] = fileparts(preop_t1);
dartel_dir = fullfile(preop_t1_path,'dartel');

if ~exist(dartel_dir,'dir')
    mkdir(dartel_dir);
end

vox_img = gunzip(elec_vox,dartel_dir);

switch ext
    case '.gz'
        preop_t1 = gunzip(preop_t1,dartel_dir);
    case '.nii'
        copyfile(preop_t1,dartel_dir,'f');
        preop_t1 = cellstr(fullfile(dartel_dir,[preop_t1_img ext]));
    otherwise
        elec_mni = [];
        disp('Please convert the Preop_T1 image into .nii or .nii.gz format first.')
        return
end

flowfield = fullfile(dartel_dir,['u_rc1',preop_t1_img,ext]);

%% dartel warping new seg

nte_path = fileparts(which('ntools_elec'));
spm_path = fileparts(which('spm'));

if ~exist(flowfield,'file')
    %% create deformation flow field and deform the image
    nrun = 1; % enter the number of runs here
    jobfile = {fullfile(nte_path,'dartel','dartel_warp_spm12_job.m')};
    jobs = repmat(jobfile, 1, nrun);

    inputs = cell(15, nrun);
    for crun = 1:nrun

        inputs{1, crun} = preop_t1; % Segment: Volumes - cfg_files
        inputs{2, crun} = {[spm_path,'/tpm/TPM.nii,1']}; % Segment: Tissue probability map - cfg_files
        inputs{3, crun} = {[spm_path,'/tpm/TPM.nii,2']}; % Segment: Tissue probability map - cfg_files
        inputs{4, crun} = {[spm_path,'/tpm/TPM.nii,3']}; % Segment: Tissue probability map - cfg_files
        inputs{5, crun} = {[spm_path,'/tpm/TPM.nii,4']}; % Segment: Tissue probability map - cfg_files
        inputs{6, crun} = {[spm_path,'/tpm/TPM.nii,5']}; % Segment: Tissue probability map - cfg_files
        inputs{7, crun} = {[spm_path,'/tpm/TPM.nii,6']}; % Segment: Tissue probability map - cfg_files
        inputs{8, crun} = {[nte_path '/dartel/Template_1.nii']}; % Run Dartel (existing Templates): Template - cfg_files
        inputs{9, crun} = {[nte_path '/dartel/Template_2.nii']}; % Run Dartel (existing Templates): Template - cfg_files
        inputs{10, crun} = {[nte_path '/dartel/Template_3.nii']}; % Run Dartel (existing Templates): Template - cfg_files
        inputs{11, crun} = {[nte_path '/dartel/Template_4.nii']}; % Run Dartel (existing Templates): Template - cfg_files
        inputs{12, crun} = {[nte_path '/dartel/Template_5.nii']}; % Run Dartel (existing Templates): Template - cfg_files
        inputs{13, crun} = {[nte_path '/dartel/Template_6.nii']}; % Run Dartel (existing Templates): Template - cfg_files
        inputs{14, crun} = {[spm_path '/canonical/single_subj_T1.nii']}; % Deformations: Image to base Id on - cfg_files
        inputs{15, crun} = vox_img; % Deformations: Apply to - cfg_files

    end

    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});

else
    
    %% deform the image with flow field
    % Deformations: Flow field - cfg_files
    % Deformations: Image to base Id on - cfg_files
    % Deformations: Apply to - cfg_files
    nrun = 1; % enter the number of runs here
    jobfile = {fullfile(nte_path,'dartel','dartel_deform_spm12_job.m')};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(3, nrun);
    for crun = 1:nrun
        inputs{1, crun} = {flowfield}; % Deformations: Flow field - cfg_files
        inputs{2, crun} = {[spm_path '/canonical/single_subj_T1.nii']}; % Deformations: Image to base Id on - cfg_files
        inputs{3, crun} = vox_img; % Deformations: Apply to - cfg_files
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});

end


toc

%% get elecs' locations

[pathstr, name, ext] = fileparts(char(vox_img)); 
warped_vox_img = fullfile(pathstr, ['w',name, ext]);
hdr_ch2 = ntools_elec_load_nifti(warped_vox_img);

s = max(unique(hdr_ch2.vol));
elec_ch2_vox = zeros(s,3);

for ii=1:s
    ind = find(hdr_ch2.vol==ii);
    if ~isempty(ind)
        [a, b, c] = ind2sub(hdr_ch2.dim(2:4)',ind);
        a_avg = mean(a);
        b_avg = mean(b);
        c_avg = mean(c);
        elec_ch2_vox(ii,:) = [a_avg b_avg c_avg];
    else
        elec_ch2_vox(ii,:) = [abs(hdr_ch2.quatern_x) abs(hdr_ch2.quatern_y) abs(hdr_ch2.quatern_z)];
        disp(['Electrode #' num2str(ii) ' not found. Probably it shares the same position with anther electrode. Reset to [0 0 0].']);
    end
end

elec_ch2_vox = [elec_ch2_vox ones(s,1)]; 
elec_ch2_ras = hdr_ch2.vox2ras*elec_ch2_vox';
elec_mni = elec_ch2_ras';
elec_mni = elec_mni(:,1:3);
