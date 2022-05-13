%% Define variables

%First some cleanup
%close all;

%Define lattice and reciprocal lattice
lat.a = 5.467;
lat.b = 5.467;
lat.c = 5.467;
lat.alpha = 0.938;
lat.beta = 0.938;
lat.gamma = 0.938;

recip = ReciprocalTransform(lat);

cstar = recip.c;

%Define variables for later use
b1 = 0;
b2 = 0;
b3 = 0;

%print some reminder text for later outputs
%output = sprintf( 'I\tJ\tK \t Screw-Proj \t  Edge Proj Para  \t  Edge Proj Perp' );

%disp(output);

%This is for some testing. Shortens the Burgers vector by this factor
prefactor = 6;

%load in Data for fitting
load1 = load('Dislocation Fitting\38_27_61_3.mat');
load2 = load('Dislocation Fitting\38_27_61_4.mat');
load3 = load('Dislocation Fitting\38_27_61_5.mat');

%prepare all the data in a single structure for fitting purposes
data.angles = [load1.data.angles; load2.data.angles; load3.data.angles];
data.values = double([load1.data.values; load2.data.values; load3.data.values]);

%creates some variables to track iterations
n=1:125;
m =0;

%% Iterate through all the reasonable Burgers vectors
for i=-2:2

    for j=-2:2
        
        for k = -2:2
            
            m = m+1;%iteration counter
            
            d = ObjectiveCoordTransform(lat, [0 0 1]); %d is dislocation line
            b = ObjectiveCoordTransform(lat, [i j k])/prefactor; %b is burgers vector
            
%             CosTheta = max(min(dot(u,b)/(norm(u)*norm(b)),1),-1);
%             ThetaInDegrees = real(acosd(CosTheta));
%             
            screw = dot(b,d).*[0 0 1]/norm(d);
            
            edge = b - screw;
            
            edge_para = norm(edge)*cross(d,edge)/norm(cross(d,edge)); %vector in the direction of u_e_parallel, with length b_e
            
            screwp = dot(screw, cstar)/norm(cstar);
            
            edgep_perp = dot(edge, cstar)/norm(cstar);
            
            edgep_para = dot(edge_para, cstar)/norm(cstar);
            
%             if(i==-2 && j==1 && k==1) 
%                 pause; 
%             end
        
            %output = sprintf([num2str(i) '\t' num2str(j) '\t' num2str(k) '\t' num2str(screwp)...
             %   '\t\t' num2str(edgep_para) '\t\t' num2str(edgep_perp)]);
           
            %disp(output);
            
            b1 =screwp;
            if not(isnan(b2)); b2 = edgep_para; else; b2 = 0; end
            b3 = edgep_perp;
            nu =0.32;
            phs1 = 0;
            phs2 = 0;
            phs3 = 0;
           
            curve = fittype(@(phs1, phs2, phs3, theta) DislocationDisplacement(theta, b1, b2, b3, nu, phs1, phs2, phs3)...
                 , 'independent','theta', 'dependent', 'u', 'coefficients', {'phs1','phs2','phs3'});
            
            options = fitoptions(curve);
            options.StartPoint = [0,0,0];
            options.Lower = [-2*pi -2*pi -2*pi];
            options.Upper = [2*pi 2*pi 2*pi];
            
            if (i==0 && j==0); continue; end %Don't fit these ones because it breaks
           
            [currentfit, gof] = fit(data.angles, data.values, curve, options);
                
            allfits(m).coeffs = coeffvalues(currentfit);
            allfits(m).fits = currentfit;
            allfits(m).gof = gof;
            allfits(m).ijk = [i j k];
            
%             plot(currentfit, data.angles, data.values);
%             
%             if allfits(m).ijk==[-2,-1,1]
%                 
%                 cftool;
%                 pause;
%                 figure; 
%             
%             end




        end
    end

end%iterate through all possible Burgers vectors to identify necessary coeeficient

%% Find the best fit based on GOF parameters


currentbest_rmse = 1; %define variables
currentbest_sse = 1000;

for i=1:125 %iterate through all the fits
    
    
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
        best_sse_ijk = allfits(i).ijk; %and the burgers vector that it corresponds to
        
    end

end
    
%% plot all the best fits

for i=1:125 %iterate through all the fits
    
    
    if isempty(allfits(i).gof) %only plot ones with goodness of fit data
        continue; 
    end

    tempgof = allfits(i).gof; %retrieve gof data for fit i; kind of a janky work around to access data 
                                %stored in a substructure
    
    
    
    if (tempgof.rmse - currentbest_rmse) < 0.00001%check if iteration is best fit
        
        figure;
        hold on;
        
        %calculate the error propogation based on error in radius,
        %normalized by radius
        error1(1:length(load1.data.angles)) =sqrt( (mean(load1.data.errors) * (1-2*nu)/(2*(1-nu))*...
            norm(ObjectiveCoordTransform(lat, allfits(i).ijk))/2/pi)^2 );
        error2(1:length(load2.data.angles)) =sqrt( (mean(load2.data.errors) * (1-2*nu)/(2*(1-nu))*...
            norm(ObjectiveCoordTransform(lat, allfits(i).ijk))/2/pi)^2 );
        error3(1:length(load3.data.angles)) =sqrt( (mean(load3.data.errors) * (1-2*nu)/(2*(1-nu))*...
            norm(ObjectiveCoordTransform(lat, allfits(i).ijk))/2/pi)^2 );
       
        
        %Plot all the data with error bars
        ebs1 = errorbar(load1.data.angles, load1.data.values, error1, '.b');
        errorbar(load2.data.angles, load2.data.values, error2, '.b');
        errorbar(load3.data.angles, load3.data.values, error3, '.b');
       
        
        %Plot the fit
        fithandle = plot(allfits(i).fits);
        fithandle.LineWidth = 1.5;
        
        xlabel('\theta (rad)','FontSize',18,'FontWeight','bold');
        ylabel('Displacement (ang.)','FontSize',18,'FontWeight','bold');
        prop = gca;
        prop.FontWeight = 'bold';
        prop.FontSize = 18;
        legend([ebs1,fithandle], 'Data', 'Fit', 'Location', 'southeast');
        legend('boxoff');
        ylim( [-0.5 0.5])
        
        %Get length of Burger's vector for reference
        blength = norm(ObjectiveCoordTransform(lat, allfits(i).ijk))/prefactor;
        
        %If short boi save figure to a file before annotating
        if(blength - 1.427)<0.1
            filename = ['Unproc-images/fit2 ' num2str(allfits(i).ijk) '.png'];
            saveas(gcf, filename); 
        end
             
        %Add some formatting for the figure
        title(['[' num2str(allfits(i).ijk) ']']);
        
        %Annotate plot with coefficients and confidence intervals
        coeffs = coeffnames(allfits(i).fits);
        coeffvals= coeffvalues(allfits(i).fits);
        ci = confint(allfits(i).fits,0.95);
        str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
        str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
        str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
        str4 = sprintf('\n Burgers Vector Length = %0.3f', blength);
        str5 = sprintf('\n RMSE of Fit = %0.3f', tempgof.rmse);
        annotation('textbox',[.15 .5 .7 .4],'String',...
            ['Coefficients (with 95% confidence bounds): ', str1, str2, str3, str4, str5], 'EdgeColor', 'none');
%         allfits(i).ijk
%         allfits(i).gof
        
        hold off;
        
        %Close all the ones with long Burgers vectors
        if ~(abs((blength-1.427))<0.1); close; end
        
        continue       
    end
    
    
    if (tempgof.sse == currentbest_sse) %check if iteration is best fit
        
        figure;
        plot(allfits(i).fits, data.angles, data.values)
        
        title(num2str(allfits(i).ijk));
        xlabel('Theta (rad)');
        ylabel('Displacement');
        
        
        %Annotate plot with coefficients and confidence intervals
        coeffs = coeffnames(allfits(i).fits);
        coeffvals= coeffvalues(allfits(i).fits);
        ci = confint(allfits(i).fits,0.95);
        str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
        str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
        str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
        annotation('textbox',[.15 .5 .7 .4],'String',...
            ['Coefficients (with 95% confidence bounds): ', str1, str2, str3], 'EdgeColor', 'none');
        
%         allfits(i).ijk
%         allfits(i).gof
        
    end
    


end
