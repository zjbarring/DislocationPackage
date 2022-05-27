function [dhkl] = InterplanarSpacing(indices, lat, latticeType)

switch(latticeType)
       
    case 'cubic'
        
        dhkl = lat.a/norm(indices);
        
    case 'tetragonal'
        
        term1 = (indices(1).^2+indices(2).^2)/lat.a.^2;
        term2 = indices(3).^2/lat.c.^2;
        dhkl = sqrt(1/(term1+term2));
        
    case 'hexagonal'
        
        term1 = 4/3*(indices(1).^2+indices(1)*indices(2)+indices(2).^2)/lat.a.^2;
        term2 = indices(3).^2/lat.c.^2;
        dhkl = sqrt(1/(term1+term2));
        
    case 'rhombohedral'
        
        term1 = norm(indices).^2*sin(lat.alpha).^2;
        term2 = (sum(indices.'*indices, 'all')- norm(indices).^2).*(cos(lat.alpha).^2-cos(lat.alpha));
        term3 = lat.a.^2.*(1-3.*cos(lat.alpha).^2+2*cos(lat.alpha).^3);
        dhkl = sqrt(term3/(term1+term2));
        
    case 'orthorhombic'
        
        term1 = indices(1).^2/lat.a.^2;
        term2 = indices(2).^2/lat.b.^2;
        term3 = indices(3).^2/lat.c.^2;
        dhkl = sqrt(1/(term1+term2+term3));
        
    case 'monoclinic'
        
        term1 = indices(1).^2/lat.a.^2;
        term2 = indices(2).^2*sin(lat.beta).^2/lat.b.^2;
        term3 = indices(3).^2/lat.c.^2;
        term4 = 2*indices(1)*indices(3)*cos(lat.beta)/lat.a/lat.c;
        dhkl = sqrt(sin(lat.beta).^2/(term1+term2+term3-term4));
        
    case 'triclinic'
        
        error('Triclinic not yet supported');
        
end