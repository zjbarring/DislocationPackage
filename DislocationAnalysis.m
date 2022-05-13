clear; close all;

%Data preparation for analyzing dilocations in Bi2WO6

Dataset = 'BWO_333';

%Load Data

%%%%%%%%%%%Saved as complex variable
% load('complexRho_GaAs.mat');
% 
% Phase = angle(complexRho_GaAs);

%%%%%%%%%%%Saved as a tiff file
Phase = Tiff2Var(strcat(Dataset,'.tiff'));

%Vector orthogonal to the radii you wish to take (pixel coords)
Pointing = [0, 0, 1];

%Image Coordinates from ImageJ
x_loc = 72;
y_loc = 62;
z_loc = 61;

clusterNum = 2; %cluster param for clusterdata

%Dislocation Core
[data, fig1, graph1] = circular_profile_3d(Phase, Pointing , x_loc, y_loc, z_loc, 3, mode(Phase,'all')); %Grab data
Prepare4Fit(data, strcat('Dislocation Fitting Data\',Dataset,'\',int2str(x_loc),'_',...
    int2str(y_loc),'_',int2str(z_loc),'_3'), clusterNum); %Data formatting and file save

[data, fig2, graph2] = circular_profile_3d(Phase, Pointing , x_loc, y_loc, z_loc, 2, mode(Phase,'all'));
Prepare4Fit(data, strcat('Dislocation Fitting Data\',Dataset,'\',int2str(x_loc),'_',...
    int2str(y_loc),'_',int2str(z_loc),'_2'), clusterNum);

[data, fig3, graph3] = circular_profile_3d(Phase, Pointing , x_loc, y_loc, z_loc, 4, mode(Phase,'all'));
Prepare4Fit(data, strcat('Dislocation Fitting Data\',Dataset,'\',int2str(x_loc),'_',...
    int2str(y_loc),'_',int2str(z_loc),'_4'), clusterNum);

%display(fig1);
%display(fig2);
%display(fig3);

% [gradx, grady, gradz] = truegradient(Phase);
% options.overwrite=1;
% saveastiff(gradx, 'true_gradx.tiff',options);
% saveastiff(grady, 'true_grady.tiff',options);
% saveastiff(gradz, 'true_gradz.tiff',options);

