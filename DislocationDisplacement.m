function [u] = DislocationDisplacement(theta, b1, b2, b3, nu, phs1, phs2, phs3)
%reads in some polar space, 3 coefficients reperesenting burgers vector
%projections, and 3 phase shifts to align the data.

u = b3/2/pi* ((theta-phs3)+(sin(2*(theta-phs3))/4/(1-nu))) + ...
    b2/2/pi*((cos(2*(theta-phs2))/4/(1-nu))) + b1/2/pi*(theta-phs1);

end

