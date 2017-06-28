
function [elec_txt, mri_mat] = meta_file_MRI_recon(pt_ID)

switch pt_ID
    case 'NY451'
        elec_txt = '/space/mdeh1/5/halgdev/projects/nyuproj/loc/NY451/NY451_fmri_050614_coor_T1_2014-12-22.txt';
        mri_mat  = '/space/mdeh1/5/halgdev/projects/nyuproj/loc/NY451/NY451_fmri_050614_rh_pial_surf.mat';
    case 'NY439'
        elec_txt = '/space/mdeh1/5/halgdev/projects/nyuproj/loc/NY439/NY439_fmri_020614_coor_T1_LH_2015-06-19.txt';
        mri_mat  = '/space/mdeh1/5/halgdev/projects/nyuproj/loc/NY439/NY439_fmri_020614_lh_pial_surf.mat';
    case 'NY523'
        elec_txt = '/space/mdeh5/1/halgdev/projects/pdelprato/Results/SL/NY523_Meso/MRI_recon_AddMeso/NY523_fmri_062014_coor_T1_LH_2016-02-21.txt';
        mri_mat  = '/space/mdeh5/1/halgdev/projects/pdelprato/Results/SL/NY523_Meso/MRI_recon_AddMeso/NY523_fmri_062014_lh_pial_surf.mat';
end