function [m] = load_material(filename)
% Read the material properties from the specified file.  Currently, we fill
% in the material structure as follows:
%
%   m.rho               -> scalar
%   m.C                 -> 6x6 matrix
%   m.Cnl               -> 6x6x6 tensor
%   m.elastic_type      -> 'general', 'orthotropic', 'isotropic' (derived from m.C)
%   m.nonlinearity_type -> 'linear', 'nonlinear_classical' (specified in material file)

fid = fopen(filename, 'rt');


% Read the mass density.
myline = fgetl(fid);
while ( ~strcmp(myline, '$Density') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.rho = fscanf(fid, '%g');
  if ( ~strcmp(fgetl(fid), '$EndDensity') )
    error('Material file has a corrupt $Density section.')
  end
end


% Read the permeability.
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$Permeability') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.mu = fscanf(fid, '%g');
  if ( ~strcmp(fgetl(fid), '$EndPermeability') )
    error('Material file has a corrupt $Permeability section.')
  end
end


% Read the permittivity.
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$Permittivity') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.epsilon = fscanf(fid, '%g');
  if ( ~strcmp(fgetl(fid), '$EndPermittivity') )
    error('Material file has a corrupt $Permittivity section.')
  end
end


% Read the (symmetric) stiffness tensor.
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$LinearElasticConstants') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.C = zeros(6,6);
  for row = 1:6
  
    % Read the stars.
    fscanf(fid, '%s', row-1);
  
    % Read the non-stars.
    m.C(row,row:end) = fscanf(fid, '%g', 6-row+1);
  
  end

  % Create symmetric stiffness tensor.
  m.C = m.C + triu(m.C,1)';

  % Set linear elastic material type depending on the values in the stiffness
  % tensor.
  m.elastic_type = 'general';
  
  if ( all(all(m.C(1:3,4:6) == 0)) ...
       && (m.C(4,5) == 0) ...
       && (m.C(4,6) == 0) ...
       && (m.C(5,6) == 0)  )
  
    m.elastic_type = 'orthotropic';
  
  end
  
  if (   (m.C(1,2) == m.C(1,3)) ...
      && (m.C(1.3) == m.C(2,3)) ...
      && (m.C(4,4) == m.C(5,5)) ...
      && (m.C(5,5) == m.C(6,6)) ...
      && (m.C(1,1) == m.C(1,2) + m.C(4,4))  )
  
    m.elastic_type = 'isotropic';
  
  end

  fgetl(fid); % dummy to get newline
  myline = fgetl(fid);
  if ~strcmp(myline, '$EndLinearElasticConstants')
    error('Material file has a corrupt $EndLinearElasticConstants section.')
  end

end


% Read the nonlinearity type.
m.nonlinearity_type = 'linear'; % Default = linear
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$NonlinearityType') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.nonlinearity_type = fgetl(fid);
  if ~strcmp(fgetl(fid), '$EndNonlinearityType')
    error('Material file has a corrupt $EndNonlinearityType section.')
  end
end


% Read the third order elastic constants (for nonlinearity).
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$NonlinearElasticConstants') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)

  m.Cnl = nan(6,6,6);
  for i=1:6
  
    % Skip rows full of stars.
    for row = 1:i-1
      fscanf(fid, '%s', 6);
    end
  
    % Read non-star data.
    for row = i:6
      fscanf(fid, '%s', row-1);
      [mydata count] = fscanf(fid, '%f', 6-row+1);
      if (count ~= 0)
        m.Cnl(i,row,row:end) = reshape(mydata, 1, 1, 6-row+1);
      else
        break
      end
    end
  
  end
  
  %TODO: check this!!!
  % Normally, we should only do 6*6*6-56 = 160 assignments, but here we do
  % 56*5 = 280 assignments.  This means we sometimes assign twice, and we should
  % check that we do not write a *different* value!
  for p = 1:6
    for q = p:6
      for r = q:6
  
         % Test if we don't have conflicts when exploiting the symmetry!
         % (We can't assign two different values to one element).
         if ~isnan(m.Cnl(p,r,q)) && (m.Cnl(p,r,q) ~= m.Cnl(p,q,r))
           error('Changing value!')
         end
         if ~isnan(m.Cnl(r,p,q)) && (m.Cnl(r,p,q) ~= m.Cnl(p,q,r))
           error('Changing value!')
         end
         if ~isnan(m.Cnl(q,p,r)) && (m.Cnl(q,p,r) ~= m.Cnl(p,q,r))
           error('Changing value!')
         end
         if ~isnan(m.Cnl(q,r,p)) && (m.Cnl(q,r,p) ~= m.Cnl(p,q,r))
           error('Changing value!')
         end
         if ~isnan(m.Cnl(r,q,p)) && (m.Cnl(r,q,p) ~= m.Cnl(p,q,r))
           error('Changing value!')
         end
  
         % All 6 permutations of (p,q,r) are equal.
         m.Cnl(p,r,q) = m.Cnl(p,q,r);
         m.Cnl(r,p,q) = m.Cnl(p,q,r);
         m.Cnl(q,p,r) = m.Cnl(p,q,r);
         m.Cnl(q,r,p) = m.Cnl(p,q,r);
         m.Cnl(r,q,p) = m.Cnl(p,q,r);
  
      end
    end
  end

  fgetl(fid); % dummy to get newline
  myline = fgetl(fid);
  if ~strcmp(myline, '$EndNonlinearElasticConstants')
    error('Material file has a corrupt $EndNonlinearElasticConstants section.')
  end

end

fclose(fid);
