function [anatomical_text, EOI] = ntools_elec_saveAnatomical(subj,hemi,elec_text,elec_bin)

% paint the electrodes onto subject's pial surface and output the
% anatomical regions (in percentage) where each electrode locates 
%
% Usage: ntools_elec_saveAnatomical(subj,hemi,elec_bin,elec_text)
% default cortical parcellation: ?h.aparc.annot
%
% Input:
% subj: subject ID in SUBJECTS_DIR
% hemi: hemisphere (lh rh or depth)
% elec_text: subject electrode location file in T1 space
% elec_bin: subject electrode nifti image (optional)
%
% Output:
% cortical_text: anatomical regions of each electrode in percentage
%
% created by Hugh Wang, 3/11/2015, Xiuyuan.Wang@nyumc.org
%
% modified by Hugh Wang, 5/21/2015
% make elec_bin optional to avoid wrong indexing when input text is hemi
% splited
%
% to do:
% add hemi info the electrode label 

%% read in the text file and parse the G/S/D electrodes

fprintf('\n%s\n',subj);

if ~exist('elec_bin','var'), elec_bin = []; end % if no binary image input, just map the cortical electrodes

fid = fopen(elec_text);
elec_all = textscan(fid,'%s %f %f %f %s');
elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];

% Separate Grid, Strip and Depth electrodes

if isempty(char(elec_all{5}(:)))
%     g = strncmpi('G',elec_cell(:,1),1);
    d = strncmpi('D',elec_cell(:,1),1);
else
%     g = strncmpi('G',elec_all{5},1);
    d = strncmpi('D',elec_all{5},1);
end

elec_gs = elec_cell(~d,:);
elec_depth = elec_cell(d,:);


%% process with G/S

PathName = fileparts(elec_text);if isempty(PathName), PathName = '.'; end;
cfg = [];
cfg.outdir = [PathName '/labels'];
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end;
    
hippo_elec = cell(1);
entorhinal_elec = cell(1);
% lateral_frontal_elec = cell(1);
insula_elec = cell(1);
amygdala_elec = cell(1);

if strcmpi(hemi,'depth'), elec_gs = []; end % for both hemi, jump to depth part

if ~isempty(elec_gs)
    % load pial surface
    [surf] = fs_load_subj(subj,hemi,'pial');
    [surf] = fs_calc_triarea(surf);
    surf.coords = surf.vertices;
    
    PathName = fileparts(elec_text); if isempty(PathName), PathName = '.'; end;
    cfg = [];
    cfg.subject = subj;
    cfg.surf= surf;
    cfg.elec_names = elec_gs(:,1);
    cfg.elec_coords = cell2mat(elec_gs(:,2:4));
    cfg.hemi = hemi;
    cfg.outdir = [PathName '/labels'];
    cfg.fsavg = 0;
    
    annotfile = ntools_elec_saveAnnot(cfg);
    
    %% loading ?h.aparc.annot and get the region percentage of each elec
    
    % read annotation
    [~, elec_label, elec_colortable] = fs_read_annotation(annotfile);
    
    [~, label, colortable] = fs_read_annotation([getenv('SUBJECTS_DIR'),subj,'/label/',hemi,'.aparc.split_STG_MTG.annot']);
    
    
    %%   
    for i=1:size(elec_gs,1)
        note = [];
        % if only one electrode available, use the 2nd because it was
        % overwritten in the annotation file
        if size(elec_gs,1)==1
            elec_nbrs = elec_label==elec_colortable.table(2,5);
        else
            elec_nbrs = elec_label==elec_colortable.table(i,5);
        end
        
        % didn't find neighbours in the label, happens when the electrode share
        % the same location with others
        if sum(elec_nbrs)==0
            elec_gs(i,5) = {[]};
            continue;
        end
        
        total_area = surf.vertex_area(elec_nbrs);
        aparc_table_idx = label(elec_nbrs);
        
        % remove 0 index
        aparc_table_idx = aparc_table_idx(aparc_table_idx>0);
        if isempty(aparc_table_idx), continue; end;
        
        aparc_idx_uniq = unique(aparc_table_idx);
        for j=1:length(aparc_idx_uniq)
            aparc_region = colortable.struct_names{colortable.table(:,5)==aparc_idx_uniq(j)};
            aparc_area = sum(total_area(aparc_table_idx==aparc_idx_uniq(j)));
            aparc_area_ratio = aparc_area/sum(total_area)*100;
            
            note = [note, sprintf('%0.2f%% %s ',aparc_area_ratio,aparc_region)];
        end
        
        if sum(strcmp('parahippocampal',strsplit(note,' ')))
            hippo_elec(end+1) = elec_gs(i,1);
        end
        if sum(strcmp('entorhinal',strsplit(note,' ')))
            entorhinal_elec(end+1) = elec_gs(i,1);
        end
%         if sum(strcmp('caudalmiddlefrontal',strsplit(note,' '))) || ...
%                 sum(strcmp('rostralmiddlefrontal',strsplit(note,' '))) || ...
%                 sum(strcmp('parsopercularis',strsplit(note,' '))) || ...
%                 sum(strcmp('parstriangularis',strsplit(note,' '))) || ...
%                 sum(strcmp('parsorbitalis',strsplit(note,' '))) || ...
%                 sum(strcmp('lateralorbitofrontal',strsplit(note,' '))) || ...
%                 sum(strcmp('frontalpole',strsplit(note,' ')))
%             lateral_frontal_elec(end+1) = elec_gs(i,1);
%         end
        if sum(strcmp('insula',strsplit(note,' ')))
            insula_elec(end+1) = elec_gs(i,1);
        end
        
        elec_gs(i,5) = cellstr(note);
        clear elec_nbrs total_area aparc*
    end

end

%% process depth electrodes 
if ~isempty(elec_depth) && ~isempty(elec_bin)
    hdr = ntools_elec_load_nifti(elec_bin);
    depth_row = find(d);
    
    % check the orientation of elec_bin
    [~,msg] = unix(sprintf('mri_info %s',elec_bin));
    orientation = regexp(msg,'Orientation\s+:\s+\w+','match');
    if ~isempty(orientation)
        orientation  = strtrim(orientation{1}(end-3:end));
        if ~strcmpi(orientation,'RAS') && ~strcmpi(orientation,'LAS')
            cont = input(sprintf('The orientation of the elec_bin is %s, this is unusual. Continue? [y/n]: ',...
                orientation),'s');
            if strcmpi(cont,'n'), return; end;
        end
    else
        disp('orientation is not available. please check elec_bin image');
        return;
    end
    
    % convert aseg
    aseg_mgz = [getenv('SUBJECTS_DIR'),subj,'/mri/aparc+aseg.mgz'];
    aseg_nii = [cfg.outdir,'/aparc+aseg.nii.gz'];
    if ~exist(aseg_nii,'file')
        [status,msg] = unix(sprintf('mri_convert --out_orientation %s %s %s',orientation,aseg_mgz,aseg_nii));
        if status, disp(msg); return; end;
    end
    
    aseg = ntools_elec_load_nifti(aseg_nii);

    [seg_idx, seg_name] = xlsread('/home/wangx11/matlab/ntools_elec/aparc_aseg_idx_name.xlsx');

    for k=1:length(depth_row)
        note = [];
        
        seg_num = aseg.vol(hdr.vol==depth_row(k));
        
        % continue if seg_num is empty
        if isempty(seg_num), elec_depth(k,5) = {[]}; continue; end;
        
        unique_seg_num = unique(seg_num);
        for m=1:length(unique_seg_num)
            seg_num_vox = sum(seg_num==unique_seg_num(m));
            note = [note, sprintf('%.2f%% %s ',100*seg_num_vox/27, seg_name{seg_idx==unique_seg_num(m)})];
        end
        % find out hippocampus depth elecs
        if sum(strcmp('Left-Hippocampus',strsplit(note,' '))) || ...
           sum(strcmp('Right-Hippocampus',strsplit(note,' ')))
            hippo_elec(end+1) = elec_depth(k,1);
        end
        % find out amygdala depth elecs
        if sum(strcmp('Left-Amygdala',strsplit(note,' '))) || ...
           sum(strcmp('Right-Amygdala',strsplit(note,' ')))
            amygdala_elec(end+1) = elec_depth(k,1);
        end        
        % find out insula depth elecs
        if sum(strcmp('ctx-lh-insula',strsplit(note,' '))) || ...
           sum(strcmp('ctx-rh-insula',strsplit(note,' ')))
            insula_elec(end+1) = elec_depth(k,1);
        end
        
        elec_depth(k,5) = cellstr(note);
    end 
else
    elec_depth = [];
end

%% save into text file

hippo_elec{1} = length(hippo_elec)-1;
entorhinal_elec{1} = length(entorhinal_elec)-1;
% lateral_frontal_elec{1} = length(lateral_frontal_elec)-1;
insula_elec{1} = length(insula_elec)-1;
amygdala_elec{1} = length(amygdala_elec)-1;

anatomical_text = [PathName,'/',subj,'_T1_',hemi,'_split_STG_MTG_AnatomicalRegions.txt'];
elec_cell = [elec_gs;elec_depth];
ntools_elec_savetxt(anatomical_text,elec_cell);

fid = fopen(anatomical_text,'a');
fprintf(fid,'\n');
fprintf(fid,'%% Total number of electrodes %d\n',size(elec_cell,1));

fprintf(fid,'%% Number of hippocampus electrodes %d\n',hippo_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(hippo_elec), fprintf(fid,'%s ',hippo_elec{ll}); end; fprintf(fid,'\n');

fprintf(fid,'%% Number of entorhinal cortex electrodes %d\n',entorhinal_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(entorhinal_elec), fprintf(fid,'%s ',entorhinal_elec{ll}); end; fprintf(fid,'\n');

% fprintf(fid,'%% Number of lateral frontal electrodes %d\n',lateral_frontal_elec{1});
% fprintf(fid,'%% ');
% for ll = 2:length(lateral_frontal_elec), fprintf(fid,'%s ',lateral_frontal_elec{ll}); end; fprintf(fid,'\n');

fprintf(fid,'%% Number of insula electrodes %d\n',insula_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(insula_elec), fprintf(fid,'%s ',insula_elec{ll}); end; fprintf(fid,'\n');

fprintf(fid,'%% Number of amygdala electrodes %d\n',amygdala_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(amygdala_elec), fprintf(fid,'%s ',amygdala_elec{ll}); end; fprintf(fid,'\n');

fclose(fid);


% EOI = [length(elec_cell),hippo_elec{1},entorhinal_elec{1},lateral_frontal_elec{1}];
EOI = [length(elec_cell),hippo_elec{1},entorhinal_elec{1},insula_elec{1},amygdala_elec{1}];


