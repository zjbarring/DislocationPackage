function [bestFit] = selectFits(data, params)
%selectFits Check specific Burgers vector
%   Detailed explanation goes here

finished = false;
burgers = input("Please enter the Burgers Vector. \n");

while finished == false

    d = ObjectiveCoordTransform(params.lat,  params.Dsl_Line); %d is dislocation line
    d = d/norm(d);
    b = ObjectiveCoordTransform(params.lat, burgers)/params.prefactor; %b is burgers vector

    screw = dot(b,d).*d; %Screw component in vector form

    edge = b - screw; %Edge component is mixed minus screw

    %vector in the direction of u_e_parallel, with length b_e
    edge_para = norm(edge)*cross(d,edge)/norm(cross(d,edge)); 

    %Calculate projections for the displacement components
    screwp = dot(screw, params.G_hkl)/norm(params.G_hkl);

    edgep_perp = dot(edge, params.G_hkl)/norm(params.G_hkl);

    edgep_para = dot(edge_para, params.G_hkl)/norm(params.G_hkl);


    b1 =screwp;
    if not(isnan(edgep_para)); b2 = edgep_para; else; b2 = 0; end
    b3 = edgep_perp;

    %initial values
    phs1 = 0;
    phs2 = 0;
    phs3 = 0;
    vshift = 0;


    %Set up this curve for fitting
    curve = fittype(@(phs1, phs2, phs3, vshift, theta) DislocationDisplacement(theta, b1, b2, b3, params.nu, phs1, phs2, phs3, vshift)...
         , 'independent','theta', 'dependent', 'u', 'coefficients', {'phs1','phs2','phs3','vshift'});

    options = fitoptions(curve);
    options.StartPoint = [0,0,0,0];
    options.Lower = [-2*pi -2*pi -2*pi -norm(params.G_hkl)/2];
    options.Upper = [2*pi 2*pi 2*pi norm(params.G_hkl)/2];

    [currentfit, gof] = fit(data.angles, data.values, curve, options);

    bestFit.coeffs = coeffvalues(currentfit);
    bestFit.fits = currentfit;
    bestFit.gof = gof;
    bestFit.ijk = burgers;
    
    PlotFits(data, bestFit, params);
    
    contCheck = input("Please enter new for another Burgers Vector, Edit to change parameters or Done to finish. \n", 's');
    
    contCheck = lower(contCheck);
    
    switch contCheck
        
        case "done"
            finished = true;
                  
        case "new"
            burgers = input("Please enter the Burgers Vector. \n");
            
        case "edit"
            disp("Type dbcont to continue");
            keyboard;
            burgers = input("Please enter the Burgers Vector. \n");
            
    end


    
end

end

