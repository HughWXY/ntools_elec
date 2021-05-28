# Electrode Localization Toolbox Manual

developed by Hugh Wang (Xiuyuan.Wang at nyumc.org) and Jingyun Chen (Jingyun.Chen at nyulangone.org), Department of Neurology, New York University School of Medicine

**Attention**:
>The Electrode Localization Toolbox was tested and only runs under 64bit Lunix or OSX (Thanks to David Groppe  for compiling the mex files under OSX 10.9.5 / 2015-6-3)
> 
>Preferred software version (or newer): Matlab R2013a, FSL 5.0, and FREESURFER 5.3.0 
>
>Please add SPM (currently only works with SPM12) and ntools_elec folder into Matlab path
>
>Please install FSL and FREESURFER properly

## Autocoregistration

Autocoregistration is a program that coregisters the post-implantation image with the pre-implantation image using fsl commands to coregister, match to standard MNI images, and skull strip. Input can be both .nii and .mgz files. It also using Freesurfer subject aseg file to remove the cerebellum so that it will give a clear view for the electrodes on inferior temporal lobe with 3D rendering in MRIcro.

Run ntools_elec_autocoreg In MATLAB command window. In the pop-ups:
>
>Select Freesurfer T1.mgz
>
>Select $SUBJECTS_DIR/subject/mri/T1.mgz
>
>Select $SUBJECTS_DIR/subject/mri/aseg.mgz (if available)
>
>Select elec_MRI_T1 image
>
>Select elec_MRI_T2 image (if available)

This program takes a few minutes. An elec_preop.nii.gz file should be written, as well as elec_preop.mat, elec_preop_brain.nii.gz, elec_preop_brain_cortex.nii.gz, T1.nii.gz.

## Electrode Localization

Electrode localization is composed of two parts: manual localization of some of the electrodes and automated projection of the rest of the electrodes.

### Manual localization

![manual localization interface](images/elec_manual-fsleyes.png)

Initial electrodes voxel coordinates are given from Fsleyes in above example. Use the crosshair to find the center of the electrode (best approximation) and enter the coordinates found in red box in above picture into the electrodes coordinates text or excel file.

The idea here is to localize the following:

**Depth electrodes: the first and last electrode.**

**Strip electrodes: every electrode.**

**Grid electrodes: any 2 or 3 electrodes in the grid. Usually aim for corner electrodes.**

**Notes**:

1. If the edges of grid are clear on coregistered mri, it is always preferred to select all 4 corner electrodes in either clockwise or anticlockwise order.

2. If you want to use 3 initial points for the grid, make sure these 3 points can shape a rectangular triangle, and the perpendicular foot has to be in the middle of the three.

3. If there are 2 grids labeled with the same letter, or one grid broken up into several parts, treat them as separate grids and label them differently, for example splitting GA1 to GA64 as GA1 to GA32 and GB1 to GB32 (which is actually GA33~64), adjusting the numbering of the rest grids. After a final image is created, edit the output electrodes coordinates text file and manually adjust each electrode to its original label according to the sketch.

An example of initial electrodes coordinates text:

![example electrode initial excel file](images/initial_text.png)

In the 5th column of the initial text file, indicate what type of the electrode it is (G for grid, D for depth and S for strip).

Alternatively, manual localization can also be done with [3D Slicer](https://www.slicer.org/), which allows direct pinpointing electrodes in 3D viewer. See step-by-step instructions [here](Manual_Localization_Instructions.pdf).

### Automated Projection

Run ntools_elec in Matlab command window. In the pop-ups:

1. Select the subject's FREESURFERRECON folder

![Picture1](images/Picture1.png)

2. Select the text or xls(xlsx) file that has initial electrode coordinates

![Picture2](images/Picture2.png)

3. Select the text or xls(xlsx) file that has the removed electrodes. This could be handy in some cases that a few electrodes has been removed from a 8*8 grid. Instead splitting the grid into several parts, just calculate the grid as a whole and exclude the removed electrodes

4. choose the preop T1 image

![Picture3](images/Picture3.png)

5. Select the hemisphere where electrodes lying on ( lh, rh, or both )

6. Input the inter-electrode distance (mm) (Press Enter for default distance 10mm)

7. Select the size of the grid

![Picture6](images/Picture6.png)

   Verify the orientation of the initial grid points ( If you pick up 3 or 4 initial points for the grid, this window won't pop up)

![Picture6.5](images/Picture6.5.png)

8. Calculating the mean and standard deviation of distance between nearest two electrodes. Program will automatically select the best projection with standard deviation less than 1mm and distance mean as close as the number typed in step 6 as possible. If no result satisfies the criteria, manually input the NO. of result you think is good enough or type number '0' to exit and recheck the initial text file.

![Picture7](images/Picture7.png)

9. Viewing the results on the brain surface (see **ntools_elec_plot** below)

### Outputs

lh.aparc.a2009s.annot: Freesurfer subject Destrieux Atlas annotation in surface space

lh.aparc.split_STG_MTG.annot: Freesurfer subject Desikan Atlas with split (rostral, middle, caudal) superior/medial temporal gyri in surface space

[subj]\_coor\_MNI\_[date].txt: electrode coordinates in MNI volume and surface space

[subj]\_coor\_T1\_[date].txt: electrode coordinates in subject SURFACE space (**not fit with T1.nii.gz**)

[subj]\_elec\_bin\_T1\_[date].nii.gz: electrode image in subject preop T1 VOLUME space(**not fit with \_pial\_surf.mat**), with intensity corresponding to the row number of coor_T1 text file, indicating the electrode labels. The coordinates of electrodes from elec_bin is usually different from the coor_T1 text file, but this difference is consistent across electrodes, which is the shift between volume space and surface space.

[subj]\_lh\_pial\_surf.mat: subject triangulated surface mat file

T1.nii.gz: subject preop T1 MPRAGE preprocessed by Freesurfer

## Useful programs for electrode MRI visualization

[MRICRCO](http://www.sph.sc.edu/comd/rorden/mricro.html) (PC, Linux & MAC)

[MRICRON](http://www.sph.sc.edu/comd/rorden/mricron) (PC, Linux & MAC)

[Fsleyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes) (Linux & MAC)

[3DSlicer](http://www.slicer.org/) (PC, Linux, MAC)

## Electrode Visualization

### ntools_elec_plot

![Picture_plot](images/Picture_plot.png)

a stand-alone Matlab script for viewing the electrodes on brain surface. In MATLAB command window start _ntools_elec_plot_

>Select the [subj]\_coor\_T1\_[date].txt file
>
>Select the surface mat file
>
>Choose to plot only the grid, strips, depth or both grid and stips.
>
>Select if (or not) to show the electrode labels

Select if (or not) to plot with segmentation annotations (see below for plotting with annotation)

![Picture_plot_aparc](images/Picture_plot_aparc.png)

Select if (or not) to save the figures into images. If yes, it will create a folder "images" under MATLAB current working directory and save the plots there

### ntools_elec_plotGroup

![Picture_plot_group](images/Picture_plot_group.png)

a stand-alone Matlab script for viewing certain electrodes on brain surface. It requires a different electrode coordinates text file: in the fifth column, instead of G/D/S, there are numbers e.g. 1~8, so you can choose to plot 2,4,5. Numbers can be any positive integers. Example text file shows below:

![Picture_text_group](images/Picture_text_group.png)

### ntools_elec_saveAnatomical

ntools_elec_saveAnatomical paints the electrodes onto subject’s pial surface and output the anatomical regions (in percentage) where each electrode locates. A full FREESURFER reconstruction is required for this feature.

The input variable ‘hemi’ can only be ‘lh’ or ‘rh’ or ‘depth’. ‘lh’ and ‘rh’ will find the anatomical regions percentage for cortex electrodes (grids and strips), and ‘depth’ option will find the percentage for depth electrodes. The electrode coordinates text file from ntools_elec.m is combined left and right hemisphere, therefore a manual split into two separate text files is necessary.

Usage: ntools_elec_saveAnatomical(subjID, hemi, elec_coor_text, elec_voxel_img)

e.g.:
>
>ntools_elec_saveAnatomrical(‘NY001’,’lh’,’NY001_lh_coor.txt’)
>
>ntools_elec_saveAnatomrical(‘NY001’,’rh’,’NY001_rh_coor.txt’)
>
>ntools_elec_saveAnatomrical(‘NY001’,’depth’,’NY001_lh_rh_coor.txt’,’NY001_elec_bin.nii’)

Output:

[subj]\_T1\_lh(rh)\_split\_STG\_MTG\_AnatomicalRegions.txt: anatomical percentage coverage of each electrode, using aparc.split_STG_MTG.annot in default

### other types of grid supported

There are two more types of grid supported in ntools_elec: EG and MG

EG stands for experimental grid, usually in square or rectangular shape but the inter-electrode distance is much smaller than a regular grid as well as the size of the electrode. If the initial grid points are labeled 'EG' in the 5th column, it will be plotted as 1mm diameter sphere in magenta color. A typical EG grid has 16x8 electrodes, with 3mm inter-electrode distance.  

![Picture_text_group](images/NY705_T1_GS_lateral_lh.png)

MG stands for meso grid, which is PMT model 2110-128-021, an 8x8 regular grid with another 64 channels in between regular electrodes. Users only need to select the initial positions for the regular grid and the script will calculate for all the rest ones.

![Picture_text_group](images/NY717_T1_GS_lateral_lh.png)



