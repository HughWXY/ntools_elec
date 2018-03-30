% List of open inputs
% Segment: Volumes - cfg_files
% Segment: Tissue probability map - cfg_files
% Segment: Tissue probability map - cfg_files
% Segment: Tissue probability map - cfg_files
% Segment: Tissue probability map - cfg_files
% Segment: Tissue probability map - cfg_files
% Segment: Tissue probability map - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Deformations: Image to base Id on - cfg_files
% Deformations: Apply to - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/home/wangx11/matlab/ntools_elec/dartel/dartel_warp_spm12_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(15, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Volumes - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability map - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability map - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability map - cfg_files
    inputs{5, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability map - cfg_files
    inputs{6, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability map - cfg_files
    inputs{7, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability map - cfg_files
    inputs{8, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{9, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{10, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{11, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{12, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{13, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{14, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Image to base Id on - cfg_files
    inputs{15, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Apply to - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
