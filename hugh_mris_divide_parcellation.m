function hugh_mris_divide_parcellation(varargin)

% split STG and MTG in lh/rh.aparc.annot and save as
% lh/rh.aparc.split_STG_MTG.annot

if isempty(varargin)
    subj = input('subject ID: ','s');
elseif nargin==1
    subj = varargin{1};
else
    disp('only one subject at a time')
    return
end
    

sdir = getenv('SUBJECTS_DIR');
hemi = {'lh','rh'};

for i=1:length(hemi)
    fname = ['label/',hemi{i},'.aparc.split_STG_MTG.annot'];
    outfile = fullfile(sdir,subj,fname);
    
    cmd = sprintf('mris_divide_parcellation %s %s aparc.annot /home/wangx11/matlab/ntools_elec/splittable_STG_MTG.txt %s',...
        subj,hemi{i},outfile);
    [status,msg] = unix(cmd);
    if status, disp(msg); continue; end
    
    [vertices, label, colortable] = read_annotation(outfile);
    
    % change colortable.struct_names
    colortable.struct_names = regexprep(colortable.struct_names,'superiortemporal_div1','cSTG');
    colortable.struct_names = regexprep(colortable.struct_names,'superiortemporal_div2','mSTG');
    colortable.struct_names = regexprep(colortable.struct_names,'superiortemporal_div3','rSTG');
    
    colortable.struct_names = regexprep(colortable.struct_names,'middletemporal_div1','cMTG');
    colortable.struct_names = regexprep(colortable.struct_names,'middletemporal_div2','mMTG');
    colortable.struct_names = regexprep(colortable.struct_names,'middletemporal_div3','rMTG');
    
    % change colortable.table
    idx1 = find(strcmp(colortable.struct_names,'mSTG'));
    mSTG_id = colortable.table(idx1,5);
    colortable.table(idx1,1:3) = [24,200,24];
    colortable.table(idx1,5) = colortable.table(idx1,1:4)*[1;2^8;2^16;2^24];
    
    idx2 = find(strcmp(colortable.struct_names,'mMTG'));
    mMTG_id = colortable.table(idx2,5);
    colortable.table(idx2,1:3) = [110,160,220];
    colortable.table(idx2,5) = colortable.table(idx2,1:4)*[1;2^8;2^16;2^24];  
    
    % replace ID in label
    label(label==mSTG_id) = colortable.table(idx1,5);
    label(label==mMTG_id) = colortable.table(idx2,5);
    
    % replace ID in vertices
    vertices(vertices==mSTG_id) = colortable.table(idx1,5);
    vertices(vertices==mMTG_id) = colortable.table(idx2,5);
    
    % save new annotation
    write_annotation(outfile,vertices,label,colortable);
    
end
