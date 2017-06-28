function NYU_ntools_elec_autoplot_XJ(fname,surfmat,aparc)

% automatically save elec images with/without labels/aparc in background
% only plot one hemisphere at a time.

if exist(fname,'file')
    fid = fopen(fname);
    elec_all = textscan(fid,'%s %f %f %f %s');
%     elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];
    label = elec_all{5};
else
    disp('No electrode was found. Please check you input text file.')
    return
end

unique_label = unique(label);

if ismember('G',unique_label) && ismember('S',unique_label) && ~ismember('M',unique_label)
    plt = 'GS';
elseif ismember('G',unique_label) && ~ismember('S',unique_label)
    plt = 'G';
elseif ~ismember('G',unique_label) && ismember('S',unique_label)    
    plt = 'S';
elseif ismember('G',unique_label) && ismember('S',unique_label) && ismember('M',unique_label)
    plt = 'GSM';
else
    plt = [];
end

if ~isempty(plt)
    ntools_elec_plot_AddMeso(fname,surfmat,'plot',plt,'saveimg',10,'showlabel',1,'aparc',[]);
    ntools_elec_plot_AddMeso(fname,surfmat,'plot',plt,'saveimg',10,'showlabel',0,'aparc',[]);
    ntools_elec_plot_AddMeso(fname,surfmat,'plot',plt,'saveimg',10,'showlabel',0,'aparc',aparc);
    ntools_elec_plot_AddMeso(fname,surfmat,'plot',plt,'saveimg',10,'showlabel',1,'aparc',aparc);
end

% if ismember('D',unique_label)
%     ntools_elec_plot(fname,surfmat,'plot','D','saveimg',10,'showlabel',1,'aparc',[]);
%     ntools_elec_plot(fname,surfmat,'plot','D','saveimg',10,'showlabel',0,'aparc',[]);  
% end
    

    