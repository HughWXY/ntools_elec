function annotfile = ntools_elec_saveAnnot(cfg)
%
% Generates a FreeSurfer annotation file of electrode locations. 
%
% Inputs:
%   cfg             structure that contains:
%   subject         subject from SUBJECTS_DIR
%   surf            surface from fs_load_subj
%   elec_names      electrode names 
%   elec_coords     electrode RAS coords
%   hemi            hemisphere
%   outdir          where all label files and color lookup table are written
%
% Outputs:
% -FS annotation file in outdir
%   no copy to $SUBJECTS_DIR/subject/label/hemi.iEEG_electrodes.annot
% -FS label files for each electrode (in outdir)
% -Color lookup table (in outdir)
%
% Example:
% cd /space/mdeh1/7/halgdev/analysis/iEEG_NYU/NY68/NY68_CMC/
% cfg = [];
% cfg.subject = 'NY68_110707';
% cfg.Ssdir= '/home/halgdev/projects/nyuproj/subjects';
% cfg.elec_names = elec_names;
% cfg.elec_coords = elec_coords;
% cfg.hemi = 'lh';
% cfg.outdir = './labels/';
% electrodes2annot(cfg)
%
% author: btq
% date: 04/08/09
%
% modified by Hugh, 03/13/15
% - force creating new outdir
% - remove the fsaverage, elec_text, Ssdir options
% - read in electrodes var instead of text file
% 
%

subject = cfg.subject;
surf = cfg.surf;
elec_names = cfg.elec_names;
elec_coords = cfg.elec_coords;
hemi = cfg.hemi;

% force to create new output dir
if exist(cfg.outdir,'dir'), unix(sprintf('rm -rf %s',cfg.outdir)); end
mkdir(cfg.outdir);

clut_file=sprintf('%s/%s_%s_CLUT.txt',cfg.outdir,hemi,subject);


%begin writing the CLUT
fidclut=fopen(clut_file,'w');
if (fidclut == -1)
     error(sprintf('ERROR: could not open %s for writing, check path and permissions',clut_file));
end
fprintf(fidclut,'#$Id %s_%s_CLUT.txt, v 1.00 %s$\n\n',hemi,subject,date);
fprintf(fidclut,'#No. Label Name:                            R   G   B   A\n\n');
%fprintf(fidclut,'0\tunknown\t\t\t0\t0\t0\t0\n');

vertex_list=[];

%loop through each electrode and find vertex at minimum distance to
%electrode RAS.  Include neighboring vertices to make electrode more
%visible (though not consistent in size, can range from 5mm-12mm)

for i=1:size(elec_coords,1)
        b_z=abs(surf.coords(:,3)-elec_coords(i,3));
        b_y=abs(surf.coords(:,2)-elec_coords(i,2));
        b_x=abs(surf.coords(:,1)-elec_coords(i,1));
        d=sqrt((b_x.^2+b_z.^2+b_y.^2)); %edist
        cindex=find(d==min(d));
        vertex_list(end+1)=cindex;
        dvert=cindex-1;%subtract 1 to get back to FS 0-based indexing
        fid=fopen(sprintf('%s/%s.%s.label',cfg.outdir,hemi,elec_names{i}),'w');
        fprintf(fid,'#!ascii label , from subject %s vox2ras=TkReg coords=pial\n',subject);
        fprintf(fid,'%d\n',sum(surf.nbrs(cindex,:)>0)+1);
        fprintf(fid,'%g  %g  %g  %g 0.000000\n',dvert,surf.coords(cindex,1),...
            surf.coords(cindex,2),...
            surf.coords(cindex,3));
        for nvert=1:length(surf.nbrs(cindex,:))
            if surf.nbrs(cindex,nvert) > 0
                fprintf(fid,'%g  %g  %g  %g 0.000000\n',...
                    surf.nbrs(cindex,nvert)-1,...
                    surf.coords(surf.nbrs(cindex,nvert),1),...
                    surf.coords(surf.nbrs(cindex,nvert),2),...
                    surf.coords(surf.nbrs(cindex,nvert),3));
            end
        end
        fclose(fid);
        fprintf(fidclut,'%g\t%s\t\t\t%3.3g\t%3.3g\t%3.3g\t0\n',...
            i-1,elec_names{i},round(255*rand(1)),round(255*rand(1)),round(255*rand(1)));       
end

% if only one electrode on the hemi, repeat it again in the CLUT file,
% otherwise mris_label2annot will complain
if size(elec_coords,1)==1
   fprintf(fidclut,'%g\t%s\t\t\t%3.3g\t%3.3g\t%3.3g\t0\n',...
            i,elec_names{i},round(255*rand(1)),round(255*rand(1)),round(255*rand(1)));
end

fclose(fidclut);
      
%Issue command to convert labels to annotation file
annotfile = [cfg.outdir,'/',hemi,'.iEEG_electrodes.annot'];
cmd=sprintf('mris_label2annot --subject %s --hemi %s --ctab %s --annot-path %s --ldir %s --no-unknown',...
    subject,hemi,clut_file,annotfile,cfg.outdir);
[status, msg] = unix(cmd);
if status, disp(msg); return; end

