function [recip] = ReciprocalTransform(lat)

%A function that reads a lattice described by a structure with fields a, b,
%c, alpha, beta, and gamma in units, and returns a structure with the reciprocal
%lattice vectors in real space, in reciprocal units
    %Used to determine projections of displacement on the momentum transfer
    %of a scattering event

recip.a = 2*pi*cross(ObjectiveCoordTransform(lat,[0 1 0] ), ObjectiveCoordTransform(lat, [0 0 1]))/...
    dot(ObjectiveCoordTransform(lat, [1 0 0]), cross(ObjectiveCoordTransform(lat,[0 1 0] ), ...
    ObjectiveCoordTransform(lat, [0 0 1])));

recip.b = 2*pi*cross(ObjectiveCoordTransform(lat,[0 0 1] ), ObjectiveCoordTransform(lat, [1 0 0]))/...
    dot(ObjectiveCoordTransform(lat, [0 1 0]), cross(ObjectiveCoordTransform(lat,[0 0 1] ), ...
    ObjectiveCoordTransform(lat, [1 0 0])));

recip.c = 2*pi*cross(ObjectiveCoordTransform(lat,[1 0 0] ), ObjectiveCoordTransform(lat, [0 1 0]))/...
    dot(ObjectiveCoordTransform(lat, [0 0 1]), cross(ObjectiveCoordTransform(lat,[1 0 0] ), ...
    ObjectiveCoordTransform(lat, [0 1 0])));
    
end