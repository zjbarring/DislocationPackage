%First some cleanup 
close all; clear

%% Define variables

%integer location of the dislocation core in question. Used to access
%correct fitting data saved in \Dislocation Fitting Data\BWO_333\##_##_##_[2-4].mat
x_loc = 72;
y_loc = 62;
z_loc = 61;

%Shortens the Burgers vector by this factor, ie turning a[111] to a/3[111]
%Graphically the slope is ~1/prefactor
params.prefactor = 4;

%Define lattice and reciprocal lattice - Currently the conventional GaAs cell
%Use the lattice used to index the reflecitons
params.lat.a = 5.49;
params.lat.b = 5.52;
params.lat.c = 17.11;
params.lat.alpha = pi/2;
params.lat.beta = pi/2;
params.lat.gamma = pi/2;

params.nu =0.4; %Poisson's ratio of the material

params.recip = ReciprocalTransform(params.lat);

params.G_hkl = params.recip.c*6;%Reciprocal lattice vector from scattering

params.Dsl_Line = [1 1 0]; %Dislocation line, normalized later in the script


%load in Data for fitting - Use concentric circle plots around a
%dislocation. The script is currently written for phase maps.
load1 = load(strcat('Dislocation Fitting Data\BWO_333\',int2str(x_loc),'_',...
    int2str(y_loc),'_',int2str(z_loc),'_2'));
load2 = load(strcat('Dislocation Fitting Data\BWO_333\',int2str(x_loc),'_',...
    int2str(y_loc),'_',int2str(z_loc),'_3'));
load3 = load(strcat('Dislocation Fitting Data\BWO_333\',int2str(x_loc),'_',...
    int2str(y_loc),'_',int2str(z_loc),'_4'));

%prepare all the data in a single structure for fitting purposes
data.angles = [load1.data.angles; load2.data.angles; load3.data.angles];
data.values = double([load1.data.values; load2.data.values; load3.data.values]);

%Convert phase values to displacement
data.values = data.values/2/pi*2.85;
 
%Paraters for plotting results
params.YRange = [-2,2] ;%Range of the Yaxis of the plots
params.ErrBars = boolean(false); %Enable or disable error bars
params.Annotate = 'true';


params.UnprocDir = 'Dislocation Fitting Data\BWO_333\Unproc-images';

%Params for anyBurgers.m
params.ijk_range = 2; %Range of Burgers indices to iteratively check
params.Burgers_ciel = 10; %Maximum allowed Burger vector length, used to favor low energy dislocations

params.rmse_margin = 0.1; %How much larger the RMSE can be than the best value and still 
                    %display; loosen this if the only fits displayed aren't
                    %physical


%% Call Necessary Functions
% bestFits = anyBurgers(data, params);
% PlotFits(data, bestFits, params);

selectFits;

%%
% %% Setup - Most edited variable are above this point
% %creates some variables to track iterations
% n=1:125;
% m =0;
% 
% 
% %Define variables for later use
% b1 = 0;
% b2 = 0;
% b3 = 0;
% 
% %print some reminder text for later outputs
% %output = sprintf( 'I\tJ\tK \t Screw-Proj \t  Edge Proj Para  \t  Edge Proj Perp' );
% 
% %disp(output);
% 
% 
% 
% 
% %% Iterate through all the reasonable Burgers vectors - Don't Edit Except Index Bounds
% for i=-2:2
% 
%     for j=-2:2
%         
%         for k = -2:2
%             
%             m = m+1;%iteration counter
%             
%             d = ObjectiveCoordTransform(lat, Dsl_Line); %d is dislocation line
%             d = d/norm(d);
%             b = ObjectiveCoordTransform(lat, [i j k])/prefactor; %b is burgers vector
%             
%             screw = dot(b,d).*d; %Screw component in vector form
%             
%             edge = b - screw; %Edge component is mixed minus screw
%             
%             %vector in the direction of u_e_parallel, with length b_e
%             edge_para = norm(edge)*cross(d,edge)/norm(cross(d,edge)); 
%             
%             %Calculate projections for the displacement components
%             screwp = dot(screw, G_hkl)/norm(G_hkl);
%             
%             edgep_perp = dot(edge, G_hkl)/norm(G_hkl);
%             
%             edgep_para = dot(edge_para, G_hkl)/norm(G_hkl);
%             
%            
%             b1 =screwp;
%             if not(isnan(edgep_para)); b2 = edgep_para; else; b2 = 0; end
%             b3 = edgep_perp;
%             
%             %initial values
%             phs1 = 0;
%             phs2 = 0;
%             phs3 = 0;
%             vshift = 0;
%                         
%             
%             %Set up this curve for fitting
%             curve = fittype(@(phs1, phs2, phs3, vshift, theta) DislocationDisplacement(theta, b1, b2, b3, nu, phs1, phs2, phs3, vshift)...
%                  , 'independent','theta', 'dependent', 'u', 'coefficients', {'phs1','phs2','phs3','vshift'});
%             
%             options = fitoptions(curve);
%             options.StartPoint = [0,0,0,0];
%             options.Lower = [-2*pi -2*pi -2*pi -norm(G_hkl)/2];
%             options.Upper = [2*pi 2*pi 2*pi norm(G_hkl)/2];
%             
%             [currentfit, gof] = fit(data.angles, data.values, curve, options);
%             
%             allfits(m).coeffs = coeffvalues(currentfit);
%             allfits(m).fits = currentfit;
%             allfits(m).gof = gof;
%             allfits(m).ijk = [i j k];
%             
% 
% 
% 
% 
%         end
%     end
% 
% end%iterate through all possible Burgers vectors to identify necessary coeeficient
% 
% %% Find the best fit based on GOF parameters
% 
% 
% currentbest_rmse = 1; %define variables to compare fit quality
% currentbest_sse = 1000;
% 
% for i=1:m %iterate through all the fits
%     
%     
%     if isempty(allfits(i).gof) %only check ones with goodness of fit data
%         continue; 
%     end
% 
%     tempgof = allfits(i).gof; %retrieve gof data for fit i; kind of a janky work around to access data 
%                             %stored in a substructure
%     
%     
%     
%     if (tempgof.rmse<currentbest_rmse) %check if iteration is best fit
%         
%         currentbest_rmse = tempgof.rmse; %save it if it is
%         best_rmse_ijk = allfits(i).ijk; %and the burgers vector that it corresponds to
%         
%     end
%     
%     if (tempgof.sse<currentbest_sse) %check if iteration is best fit
%         
%         currentbest_sse = tempgof.sse; %save it if it is
%         best_sse_ijk = allfits(i).ijk; %and the burgers vector that it corresponds to
%         
%     end
% 
% end
%     
% %% plot all the best fits
% 
% for i=1:125 %iterate through all the fits
%     
%     
%     if isempty(allfits(i).gof) %only plot ones with goodness of fit data
%         continue; 
%     end
% 
%     tempgof = allfits(i).gof; %retrieve gof data for fit i; kind of a janky work around to access data 
%                                 %stored in a substructure
%     
%     if (tempgof.rmse - currentbest_rmse) < rmse_margin%check if iteration is best fit
%         
%         figure;
%         hold on;
%         
%        
%         
%         %Plot all the data with error bars
%         if ErrBars
%             ebs1 = errorbar(data.angles, data.values, error1, '.b');
% %             errorbar(load2.data.angles, load2.data.values, error2, '.b');
% %             errorbar(load3.data.angles, load3.data.values, error3, '.b');
%             
%         else
%             ebs1 = plot(data.angles, data.values, '.b');
% %             plot(load2.data.angles, load2.data.values, '.b');
% %             plot(load3.data.angles, load3.data.values, '.b');
%         end
%         
%         %Plot the fit
%         fithandle = plot(allfits(i).fits);
%         fithandle.LineWidth = 1.5;
%         
%         xlabel('\theta (rad)','FontSize',18,'FontWeight','bold');
%         ylabel('Displacement (ang.)','FontSize',18,'FontWeight','bold');
%         prop = gca;
%         prop.FontWeight = 'bold';
%         prop.FontSize = 18;
%         legend([ebs1,fithandle], 'Data', 'Fit', 'Location', 'southeast');
%         legend('boxoff');
%         ylim( YRange )
%         
%         %Get length of Burger's vector for reference
%         blength = norm(ObjectiveCoordTransform(lat, allfits(i).ijk))/prefactor;
%         
%         %If short boi save figure to a file before annotating
%         if(blength - Burgers_ciel)<0.1
%             filename = [UnprocDir '\fit2 ' num2str(allfits(i).ijk) '.png'];
%             saveas(gcf, filename); 
%         end
%              
%         %Add some formatting for the figure
%         title(['[' num2str(allfits(i).ijk) ']']);
%         
%         %Annotate plot with coefficients and confidence intervals
%         coeffs = coeffnames(allfits(i).fits);
%         coeffvals= coeffvalues(allfits(i).fits);
%         ci = confint(allfits(i).fits,0.95);
%         str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
%         str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
%         str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
%         str4 = sprintf('\n Burgers Vector Length = %0.3f', blength);
%         str5 = sprintf('\n RMSE of Fit = %0.3f', tempgof.rmse);
%         annotation('textbox',[.15 .5 .7 .4],'String',...
%             ['Coefficients (with 95% confidence bounds): ', str1, str2, str3, str4, str5], 'EdgeColor', 'none');
% %         allfits(i).ijk
% %         allfits(i).gof
%         
%         hold off;
%         
%         %Close all the ones with long Burgers vectors
%         if (Burgers_ciel-blength)<0; close; end
%         
%         continue       
%     end
%     
%     
%     if (tempgof.sse == currentbest_sse) %check if iteration is best fit
%         
%         figure;
%         plot(allfits(i).fits, data.angles, data.values)
%         
%         title(num2str(allfits(i).ijk));
%         xlabel('Theta (rad)');
%         ylabel('Displacement');
%         
%         
%         %Annotate plot with coefficients and confidence intervals
%         coeffs = coeffnames(allfits(i).fits);
%         coeffvals= coeffvalues(allfits(i).fits);
%         ci = confint(allfits(i).fits,0.95);
%         str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
%         str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
%         str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
%         annotation('textbox',[.15 .5 .7 .4],'String',...
%             ['Coefficients (with 95% confidence bounds): ', str1, str2, str3], 'EdgeColor', 'none');
%         
% %         allfits(i).ijk
% %         allfits(i).gof
%         
%     end
%     
% 
% 
% end
