function [surf] = fs_calc_triarea(surf)
% fs_calc_triarea - compute triangle and vertex area for a fs surface structure
% 
% [surf] = fs_calc_triarea(surf)
% 
% Input:
% surf is a structure containing:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%   faces:  vertex numbers for each face (3 corners)
%   vertices: x,y,z coordinates for each vertex
%   nbrs:   vertex numbers of neighbors for each vertex
%
% Output:
% surf is a structure containing:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%   faces:  vertex numbers for each face (3 corners)
%   vertices: x,y,z coordinates for each vertex
%   nbrs:   vertex numbers of neighbors for each vertex
%   vertex_area: surface area for each vertex
%   face_area:   surface area for each face
%   avg_dist: average inter-vertex distance
%
% created:        03/29/06 Don Hagler
% last modified:  01/08/10 Don Hagler
%
% see also: fs_read_surf, fs_read_trisurf, fs_find_neighbors
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(surf,'vertices') | ~isfield(surf, 'faces')
  fprintf('%s: ERROR: input surf must contain vertices and faces\n',mfilename);
  return;
end;

% calculate surface area for each triangle face
fprintf('%s: calculating triangle & vertex areas...',mfilename); tic;
surf.face_area = zeros(surf.nfaces,1);
surf.vertex_area = zeros(surf.nverts,1);
surf.avgdist = 0;
for f=1:surf.nfaces
  v = surf.faces(f,:);
  v_vertices = surf.vertices(v,:);
  edges = [v_vertices(2,:)-v_vertices(1,:);...
           v_vertices(3,:)-v_vertices(1,:);...
           v_vertices(3,:)-v_vertices(2,:)];
  surf.face_area(f) = 0.5*norm(cross(edges(1,:),edges(2,:)));
  surf.vertex_area(v) = surf.vertex_area(v) + surf.face_area(f);
  surf.avgdist = surf.avgdist + sum(sqrt(sum(edges.^2,2)));
end;
surf.avgdist = surf.avgdist / (3*surf.nfaces);
surf.vertex_area = surf.vertex_area / 3; % each vertex shares face with 2 other vertices
t=toc; fprintf('done (%0.2f sec)\n',t);

return;

