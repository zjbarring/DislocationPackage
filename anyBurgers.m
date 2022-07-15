function [bestfits] = anyBurgers(data, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Setup - Most edited variable are above this point
%creates some variables to track iterations
m =0;


%Preallocate variables for later use
b1 = 0;
b2 = 0;
b3 = 0;

bestfits_index = 1;

%% Iterate through all the reasonable Burgers vectors - Don't Edit Except Index Bounds
for i=-params.ijk_range:params.ijk_range

    for j=-params.ijk_range:params.ijk_range
        
        for k = -params.ijk_range:params.ijk_range
            
            m = m+1;%iteration counter
            
            d = ObjectiveCoordTransform(params.lat, params.Dsl_Line); %d is dislocation line
            d = d/norm(d);
            b = ObjectiveCoordTransform(params.lat, [i j k])/params.prefactor; %b is burgers vector
            
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
            
            allfits(m).coeffs = coeffvalues(currentfit);
            allfits(m).fits = currentfit;
            allfits(m).gof = gof;
            allfits(m).ijk = [i j k];
            




        end
    end

end%iterate through all possible Burgers vectors to identify necessary coeeficient

%% Find the best fit based on GOF parameters


currentbest_rmse = 1; %define variables to compare fit quality
currentbest_sse = 1000;

for i=1:m %iterate through all the fits
    
    
    if isempty(allfits(i).gof) %only check ones with goodness of fit data
        continue; 
    end

    tempgof = allfits(i).gof; %retrieve gof data for fit i; kind of a janky work around to access data 
                            %stored in a substructure
    
    
    
    if (tempgof.rmse<currentbest_rmse) %check if iteration is best fit
        
        currentbest_rmse = tempgof.rmse; %save it if it is
        best_rmse_ijk = allfits(i).ijk; %and the burgers vector that it corresponds to
        
    end
    
    if (tempgof.sse<currentbest_sse) %check if iteration is best fit
        
        currentbest_sse = tempgof.sse; %save it if it is
%         best_sse_ijk = allfits(i).ijk; %and the burgers vector that it corresponds to
        
    end

end

for i=1:m %iterate through all the fits
    
    
    if isempty(allfits(i).gof) %only save ones with goodness of fit data
        continue; 
    end

    tempgof = allfits(i).gof; %retrieve gof data for fit i; kind of a janky work around to access data 
                                %stored in a substructure
    
    if (tempgof.rmse - currentbest_rmse) < params.rmse_margin%check if iteration is best fit
        bestfits(bestfits_index)=allfits(i);
        bestfits_index = bestfits_index+1;
    end

end

