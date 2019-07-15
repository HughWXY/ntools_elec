function elec_mni = ntools_elec_fnirt_warp(elec_vox,preop_t1)
% Spatial normalization to MNI152_1mm with FSL
% by Josh Chen: jingyun.chen@nyulangone.org
% 2019-07-15 beta version

% % for testing only
% elec_vox_old=elec_vox;
% elec_vox=fname_bin;
% preop_t1=fullfile(preop_img_path,preop_img_file);

fprintf('start fnirt warp process......\n')

tic
[preop_t1_path,preop_t1_img,ext] = fileparts(preop_t1);
fnirt_dir = fullfile(preop_t1_path,'fnirt');

if ~exist(fnirt_dir,'dir')
    mkdir(fnirt_dir);
end

vox_img = gunzip(elec_vox,fnirt_dir);

switch ext
    case '.gz'
        preop_t1 = gunzip(preop_t1,fnirt_dir);
    case '.nii'
        copyfile(preop_t1,fnirt_dir,'f');
        preop_t1 = cellstr(fullfile(fnirt_dir,[preop_t1_img ext]));
    otherwise
        elec_mni = [];
        disp('Please convert the Preop_T1 image into .nii or .nii.gz format first.')
        return
end

% extract brain mask
unix(['cd ' fnirt_dir '; first_flirt ' char(preop_t1) ' first_flirt -cort;imrm first_flirt*stage*; '...
      'convert_xfm -omat first_flirt_cort_inv.mat -inverse first_flirt_cort.mat; '...
      'flirt -in ${FSLDIR}/data/standard/MNI152_T1_1mm_brain_mask_dil.nii.gz -ref ' char(preop_t1) ' -out bmask -applyxfm -init first_flirt_cort_inv.mat; '...
      'fslmaths ' char(preop_t1) ' -mas bmask bmasked_T1;']);
  
% fnirt to MNI152
[pathstr, name, ext] = fileparts(char(vox_img)); 
unix(['cd ' fnirt_dir '; flirt -ref ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -in bmasked_T1 -omat affine_transf.mat; '...
      'fnirt --in=' char(preop_t1) ' --aff=affine_transf.mat --cout=nonlinear_transf --config=T1_2_MNI152_2mm; '...
      'applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_1mm --in=' char(preop_t1) ' --warp=nonlinear_transf --out=wT1; '... % warp T1
      'applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_1mm --in=' char(vox_img) ' --warp=nonlinear_transf --out=w' char(name) ' --interp=nn;']);

toc

%% get elecs' locations

warped_vox_img = fullfile(pathstr, ['w',name, ext, '.gz']);
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