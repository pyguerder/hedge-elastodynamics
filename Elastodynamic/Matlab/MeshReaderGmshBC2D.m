%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{Nv}, @var{VX}, @var{VY}, @var{K}, @var{EToV}, @var{BCType}] =} MeshReaderGmshBC2D (@var{filename})
%% Read the specified 2D Gmsh mesh file.
%%
%% The parameters read are as follows:
%%
%% @itemize
%% @item @var{Nv}: Number of vertices in the grid.
%% @item @var{VX}: List of x-coordinates for vertices in grid [1*Nv].
%% @item @var{VY}: List of y-coordinates for vertices in grid [1*Nv].
%% @item @var{K}:  Number of elements
%% @item @var{EToV}:   EToV(k,i) is the vertex number in (VX,VY,VZ) for vertex i in element k [K*Nfaces].
%% @item @var{BCType}: BCType(k,i) is the boundary code for face i of element k [K*Nfaces].
%% @end itemize
%%
%% Example:
%%
%% @example
%%   [Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D('square.msh')
%% @end example
%%
%% @end deftypefn

%% Author: Bart Vandewoestyne <Bart.Vandewoestyne@kuleuven-kortrijk.be>

%% References:
%%
%% [1] Jan S. Hesthaven and Tim Warburton, `Nodal Discontinuous Galerkin
%%     Methods - Algorithms, Analysis, and Applications', Springer, 2008,
%%     ISBN 978-0-387-72065-4.
%%
%% [2] http://www.geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format

function [Nv, VX, VY, K, EToV, EToMat, BCType, psources, receivers, materials] = MeshReaderGmshBC2D(filename)

fid = fopen(filename, 'rt');

receivers = [];
psources = [];
materials = cell(0);

%%%
%%% PhysicalNames section.
%%%
myline = fgetl(fid);
while ( ~strcmp(myline, '$PhysicalNames') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  nb_physical_names = fscanf(fid, '%d\n', 1);
  nb_materials = 0;
  for i = 1:nb_physical_names
    linedata = sscanf(fgetl(fid), '%d %d %s');
    physical_names(i).dimension = linedata(1);
    physical_names(i).number = linedata(2);
    physical_names(i).name = char(linedata(4:end-1))';
    if strcmp(physical_names(i).name, 'Pointsource')
      pointsource_tag = physical_names(i).number;
    end
    if strcmp(physical_names(i).name, 'Receiver')
      receiver_tag = physical_names(i).number;
    end
    % We assume thta .dat files are material files!!!
    [dummy1 dummy2 ext] = fileparts(physical_names(i).name);
    if ( strcmp(ext, '.dat') )
      nb_materials = nb_materials + 1;
      materials{nb_materials} = load_material(physical_names(i).name);
    end
  end
  if ( ~strcmp(fgetl(fid), '$EndPhysicalNames') )
    error('Mesh file has a corrupt $PhysicalNames section.')
  end
end


%%%
%%% Nodes section.
%%%
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$Nodes') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  Nv = fscanf(fid, '%d\n', 1);
  tmp = fscanf(fid, '%d %f %f %f\n', [4 Nv]);
  VX = tmp(2,:);
  VY = tmp(3,:);
  VZ = tmp(4,:);
  if ( ~strcmp(fgetl(fid), '$EndNodes') )
    error('Mesh file has a corrupt $Nodes section.')
  end
end


%%%
%%% Elements section.
%%%
frewind(fid)
myline = fgetl(fid);
while ( ~strcmp(myline, '$Elements')  && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)

  % Read the number of elements in the grid.  These are 'Gmsh elements', which
  % can also be nodes, lines, polygons,...  Currently, we are only interested
  % in elements of type 1 (2-node line, necessary for the boundary conditions)
  % and type 2 (3-node triangle).
  number_of_gmsh_elements = fscanf(fid, '%d\n', 1);
  
  % Read the elements.
  K = 0;
  nb_bc_edges = 0;
  nb_psources = 0;
  nb_receivers = 0;
  for n = 1:number_of_gmsh_elements
  
    linedata = sscanf(fgetl(fid), '%d');
  
    % Read number (index) of the n-th element in the mesh.
    elm_number = linedata(1);
  
    % Read element type that defines the geometrical type of the n-th element.
    elm_type = linedata(2);
  
    % If there are tags, read them.
    %  1st tag: number of the physical entity to which the element belongs.
    %  2nd tag: number of the elementary geometrical entity to which the element
    %           belongs.
    %  3rd tag: number of a mesh partition to which the element belongs.
    number_of_tags = linedata(3);
    tags = linedata(4:3+number_of_tags);
    physical_entity    = tags(1);
    geometrical_entity = tags(2);
    mesh_partition     = tags(3);
  
    % Read the list of the node numbers.
    % For 3-node triangles, node numbers are ordered counterclockwise, see
    % http://www.geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
    % Note bart: the online documentation seems not to match reality: if we
    % look at the .msh file, it seems like the node numbers are ordered clockwise...
    node_number_list = linedata(4+number_of_tags:end);
  

    % Element type 1: 2-node line.
    % -> read the list of edges where Boundary Conditions are specified.
    if (elm_type == 1)
      nb_bc_edges = nb_bc_edges + 1;
      bc_edges(nb_bc_edges,:) = node_number_list;
      bc_codes(nb_bc_edges) = physical_entity;
    end
  
    % Element type 2: 3-node triangle.
    if (elm_type == 2)
  
      K = K + 1;
  
      % Small bug fix because gmsh seems to order clockwise instead
      % of counter clockwise as described at
      % http://geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
      %EToV(K,:) = node_number_list;
      EToV(K,:) = node_number_list([1 3 2]);

      % Store the material for this triangle.
      EToMat(K,:) = physical_entity;
  
    end
  
    % Element type 15: 1-node point (Pointsource or Receiver).
    if (elm_type == 15) 

      if ( exist('pointsource_tag') ...
            && (physical_entity == pointsource_tag) )
        nb_psources = nb_psources + 1;
        fprintf('Pointsource detected at (%f, %f)!\n', ...
                                VX(node_number_list), VY(node_number_list));
        psources(nb_psources).coords = [VX(node_number_list) VY(node_number_list)];
        psources(nb_psources).idx = node_number_list(1);
      end
  
      % Receiver detected!
      if ( exist('receiver_tag') ...
            && (physical_entity == receiver_tag) )
        nb_receivers = nb_receivers + 1;
        fprintf('Receiver detected at (%f, %f)!\n', ...
                                VX(node_number_list), VY(node_number_list));
        receivers(nb_receivers).coords = [VX(node_number_list) VY(node_number_list)];
        receivers(nb_receivers).idx = node_number_list(1);
      end

    end
  
  end

end


% Remove point source nodes from the node list (because they are not mesh
% nodes!) and make sure that indices in EToV matrix stay consistent.
for i = 1:nb_psources
  VX(psources(i).idx) = [];
  VY(psources(i).idx) = [];
  Nv = Nv - 1;
  idx_to_decrement = EToV > psources(i).idx;
  EToV(idx_to_decrement) = EToV(idx_to_decrement) - 1;
  for j = 1:nb_receivers
    if (receivers(j).idx > psources(i).idx)
      receivers(j).idx = receivers(j).idx - 1;
    end
  end
  for j = i:nb_psources
    psources(j).idx = psources(j).idx - 1;
  end
end


% Remove receiver nodes from the node list (because they are not mesh
% nodes!) and make sure that indices in EToV matrix stay consistent.
for i = 1:nb_receivers
  VX(receivers(i).idx) = [];
  VY(receivers(i).idx) = [];
  Nv = Nv - 1;
  idx_to_decrement = EToV > receivers(i).idx;
  EToV(idx_to_decrement) = EToV(idx_to_decrement) - 1;
  for j = i:nb_receivers
    receivers(j).idx = receivers(j).idx - 1;
  end
end


% Set the Boundary Conditions matrix BCType.
%
% Note Bart: I *think* that the convention for the face numbering is
% like this:
%
%    Face 1 is between node 1 and node 2.
%    Face 2 is between node 2 and node 3.
%    Face 3 is between node 3 and node 1.

BCType = zeros(K,3);

if (nb_bc_edges > 0)

  % Face 1: between nodes 1 and 2.
  [dummy idx1 idx2] = intersect(bc_edges, EToV(:,[1 2]), 'rows');
  BCType(idx2,1) = bc_codes(idx1);
  [dummy idx1 idx2] = intersect(fliplr(bc_edges), EToV(:,[1 2]), 'rows');
  BCType(idx2,1) = bc_codes(idx1);
  
  % Face 2: between nodes 2 and 3.
  [dummy idx1 idx2] = intersect(bc_edges, EToV(:,[2 3]), 'rows');
  BCType(idx2,2) = bc_codes(idx1);
  [dummy idx1 idx2] = intersect(fliplr(bc_edges), EToV(:,[2 3]), 'rows');
  BCType(idx2,2) = bc_codes(idx1);
  
  % Face 3: between nodes 3 and 1.
  [dummy idx1 idx2] = intersect(bc_edges, EToV(:,[3 1]), 'rows');
  BCType(idx2,3) = bc_codes(idx1);
  [dummy idx1 idx2] = intersect(fliplr(bc_edges), EToV(:,[3 1]), 'rows');
  BCType(idx2,3) = bc_codes(idx1);

end

fclose(fid);
