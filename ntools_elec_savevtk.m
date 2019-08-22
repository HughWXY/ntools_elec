function coorvtk = ntools_elec_savevtk(coortxt)

% convert coodinates text files to vtk format

fid = fopen(coortxt,'r');
elec_all = textscan(fid,'%s %f %f %f %s');
fclose(fid);

% elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];
labels = elec_all{1};
X = elec_all{2};
Y = elec_all{3};
Z = elec_all{4};
types = elec_all{5};


[pname,fname] = fileparts(coortxt);
if isempty(pname), pname = pwd; end

coorvtk = fullfile(pname,[fname,'.vtk']);
fid = fopen(coorvtk,'w+');
% write vtk header
fprintf(fid,'# tvk DataFile Version 2.0\n');
fprintf(fid,'%s\n',fname);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');
fprintf(fid,'POINTS %d float\n',numel(labels));

% write electrode coordinates
for i=1:numel(labels)
    fprintf(fid,'%f %f %f\n',X(i),Y(i),Z(i));
end

% write vertices
fprintf(fid,'VERTICES 1 %d\n',numel(labels));
for i=0:numel(labels)-1
    fprintf(fid,'%d ',i);
end
fprintf(fid,'\n');

fclose(fid);


