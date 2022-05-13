function [xyz] = ObjectiveCoordTransform(lat, uvw)
%ObjectiveCoordTransform Transforms a lattice vector, uvw, into real space
%   Provide the fully define lattice in the lattice struct (angles in rad):
%
%   lattice.a   lattice.b   lattice.c
%   lattice.alpha lattice.beta  lattice.gamma
%
%   The function will transform a 1x3 vector, uvw, into objective space
%   coordinates xyz, given by a 1x3 vector.  The units of the coordiante
%   space are inherited from the lattice vector

if (isa(lat.a,'int')); lat.a = double(lat.a); end
if (isa(lat.b,'int')); lat.b = double(lat.b); end
if (isa(lat.c,'int')); lat.c = double(lat.c); end

xyz(1) = lat.a*uvw(1)*sin(lat.beta)+lat.b*uvw(2)*((cos(lat.gamma)-...
         (cos(lat.beta)*cos(lat.alpha)))/sin(lat.beta));

xyz(2) = lat.b*uvw(2)*sqrt((sin(lat.alpha)^2)-...
         (cos(lat.gamma)-cos(lat.beta)*cos(lat.alpha))^2/sin(lat.beta)^2);
     
xyz(3) = uvw(3)*lat.c+uvw(2)*lat.b*cos(lat.alpha)+uvw(1)*lat.a*cos(lat.beta);



xyz = round(xyz,10, 'decimals');

end

