clear; close all;

threshold = 0.000001;
filename = 'CBNO/June22 Beamtime Measurements/L0_after.mat';
exp_filename = 'CBNO/June22 Beamtime Measurements/L0_after_coretrace.tiff';

% %temporary
% load('Peer Edits and revisions/d33 paper/d33_displacement_E123.mat');
% mask = uint8(uu1_aa~=0);

[val] = bruteforce3(filename, threshold);

%Load reconstruction into the variable 'complexCrystal'
S = load(filename);
varName = fieldnames(S);
complexCrystal = S.(varName{1});

%% Group Points & Save results
% disl_flags = val>4;
% gridsize = size(val);
% [xx,yy,zz] = ndgrid(1:1:gridsize(1), 1:1:gridsize(2), 1:1:gridsize(3));
% 
% disl_points = [val(disl_flags), xx(disl_flags), yy(disl_flags), zz(disl_flags)];

%ind=clusterdata(disl_points, 3);



options.overwrite = 'true';

saveastiff(val, exp_filename , options);



%% Visualization:

%Whole network of dislocations in one color
figure; p = patch(isosurface(val, 4), 'FaceColor', 'yellow', 'EdgeColor', 'blue');
axis equal; axis off;

camlight;

axis equal; axis off;

p = patch(isosurface(abs(complexCrystal),threshold));%alpha(0.4);
set(p,'facecolor','red','edgecolor','none','facealpha',0.2);

% hold on
% sfs = Tiff2Var('true_gradz.tiff');
% patch(isosurface(sfs, 0.5), 'FaceColor', 'yellow', 'EdgeColor', 'black')
% hold off



%Dislocation points grouped by proximity

% Dislocations = IsolateDislocation(val, 7);

% figure;
% hold on;

% p1 = patch(isosurface(Dislocations(1).dislocation, 4), 'FaceColor', 'yellow', 'EdgeColor', 'yellow');
% p2 = patch(isosurface(Dislocations(2).dislocation, 4), 'FaceColor', 'red', 'EdgeColor', 'red');
% p3 = patch(isosurface(Dislocations(3).dislocation, 4), 'FaceColor', 'blue', 'EdgeColor', 'blue');
% p4 = patch(isosurface(Dislocations(4).dislocation, 4), 'FaceColor', 'green', 'EdgeColor', 'green');
% p5 = patch(isosurface(Dislocations(5).dislocation, 4), 'FaceColor', 'white', 'EdgeColor', 'white');
% p6 = patch(isosurface(Dislocations(7).dislocation, 4), 'FaceColor', 'black', 'EdgeColor', 'black');
% 
% axis equal; axis off;
% 
% camlight;
% 
% p = patch(isosurface(abs(complexRho_GaAs),.1));%alpha(0.4);
% set(p,'facecolor','red','edgecolor','none','facealpha',0.2);
% 
% %Create color coded tiff stack to crossreference traced dislocation and
% %phase maps; use royal color map in FIJI
% dis_ref = zeros(size(val));
% dis_ref = dis_ref + (Dislocations(1).dislocation>4).*11; %Yellow
% dis_ref = dis_ref + (Dislocations(2).dislocation>4).*15; %Red
% dis_ref = dis_ref + (Dislocations(3).dislocation>4).*4; %Blue
% dis_ref = dis_ref + (Dislocations(4).dislocation>4).*9; %Green
% dis_ref = dis_ref + (Dislocations(5).dislocation>4).*20; %White
% dis_ref = dis_ref + (Dislocations(7).dislocation>4).*1; %Black

% option.overwrite = 'true';
% saveastiff(dis_ref, 'Dislocation Cross Reference.tiff',options);

