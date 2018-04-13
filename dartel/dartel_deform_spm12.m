% List of open inputs
% Deformations: Flow field - cfg_files
% Deformations: Image to base Id on - cfg_files
% Deformations: Apply to - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/home/wangx11/matlab/ntools_elec/dartel/dartel_deform_spm12_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Flow field - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Image to base Id on - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Apply to - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
