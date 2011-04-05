% Note Bart:
%
%   Things that are set after this routine is called:
%
%     um.Tsource The index of the triangle containing this source
%                (used in ElasticsRHS2D.m)
%     um.Vpsi    The value of V times the psi-function, evaluated in the
%                (r,s)-coordinates for this source (used in ElasticsRHS2D.m)
%   
function [um] = init_pointsource2D(um, N)

global EToV
global VX VY
global V

K = size(EToV, 1);

% Find the triangle that contains the source, based on the barycentric
% technique at http://www.blackpawn.com/texts/pointinpoly/default.html
for k=1:K

  % Compute the vectors.
  x_coords = VX(EToV(k,:));
  y_coords = VY(EToV(k,:));
  v0 = [x_coords(2)-x_coords(1) y_coords(2)-y_coords(1)];
  v1 = [x_coords(3)-x_coords(1) y_coords(3)-y_coords(1)];
  v2 = [um.x0-x_coords(1) um.y0-y_coords(1)];

  % Compute the dot products.
  dot00 = dot(v0, v0);
  dot01 = dot(v0, v1);
  dot02 = dot(v0, v2);
  dot11 = dot(v1, v1);
  dot12 = dot(v1, v2);

  % Compute the barycentric coordinates.
  invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  % Check if the source is in the current triangle.
  if ( (u > 0) && (v > 0) && (u + v < 1) )
    um.Tsource = k;
    fprintf('=> The source is located in triangle %d.\n', k)
    break
  end

end


% Determine the (r,s)-coordinates of the source in the reference triangle by
% applying the inverse of the mapping (6.3) on page 172 in the DGFEM book.
% If we write the forward mapping as
%
%   x = A*r + B
%
% then the inverse of this mapping is
%
%   r = A^{-1}*(x-B)
%
% which is what is implemented below.
nA = EToV(um.Tsource,1);
nB = EToV(um.Tsource,2);
nC = EToV(um.Tsource,3);
xA = VX(nA);
xB = VX(nB);
xC = VX(nC);
yA = VY(nA);
yB = VY(nB);
yC = VY(nC);
matA = [(xB-xA)/2 (xC-xA)/2; ...
        (yB-yA)/2 (yC-yA)/2];
vectB = [(xB+xC)/2; ...
         (yB+yC)/2];
vectRS = matA\([um.x0; um.y0]-vectB);

% Once we have the (r,s)-coordinates of the source, evaluate the psi-function
% from formula (2.64) page 81 in YiFeng's PhD thesis at this point.
r_source = vectRS(1);
s_source = vectRS(2);
Np = (N+1)*(N+2)/2;
um.Vpsi = zeros(Np,1);
skip = 1;
for i=0:N
  for j=0:N-i
    if (s_source == 1)
      a = -1;
    else
      a = 2*(1+r_source)/(1-s_source) - 1;
    end
    % We use the Koornwinder-Dubiner polynomials here.  These are normalized
    % versions of the polynomials from formula (6.6) in the book.
    um.Vpsi(skip) = JacobiP(a,0,0,i)*JacobiP(s_source,2*i+1,0,j)*((0.5*(1-s_source))^i)*sqrt((i+0.5)*(i+j+1));
    skip = skip + 1;
  end
end
um.Vpsi = V*um.Vpsi;
