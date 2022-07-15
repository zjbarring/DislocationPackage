function [u] = DislocationDisplacement2(theta, burgVec,...
        nu, Dsl_Line, G_hkl)
%reads in some polar space, 3 coefficients reperesenting burgers vector
%projections, and 3 phase shifts to align the data.
%Outputs 2d array of displacement



d = Dsl_Line; %d is dislocation line
d = d/norm(d);
b = burgVec; %b is burgers vector

screw = dot(b,d).*d; %Screw component in vector form

edge = b - screw; %Edge component is mixed minus screw

%vector in the direction of u_e_parallel, with length b_e
edge_para = norm(edge)*cross(d,edge)/norm(cross(d,edge)); 

%Calculate projections for the displacement components
screwp = dot(screw, G_hkl)/norm(G_hkl);

edgep_perp = dot(edge, G_hkl)/norm(G_hkl);

edgep_para = dot(edge_para, G_hkl)/norm(G_hkl);


b1 =screwp;
if not(isnan(edgep_para)); b2 = edgep_para; else; b2 = 0; end
b3 = edgep_perp;

            
u =  b3/2/pi* (theta)+ b3/2/pi*(sin(2*(theta))/4/(1-nu)) + ...
    b2/2/pi*(cos(2*(theta))/4/(1-nu)) + b1/2/pi*(theta);

end
