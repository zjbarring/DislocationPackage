function [results, bin_display, graph] = circular_profile_3d(data, vector, location_x, location_y, location_z, radius, bkg)
%Circular plot around a point orthogonal to some vector
%   Given a dataset, vector, point in space, and a radius this function
%   will select the data in a circle around the given point, with the given
%   radius. The circle will be tilted to be normal to the given vector. Any
%   values exactly equal to zero in the dataset will be discarded. The
%   data will be returned in a structure with fields 'angles' & 'values'
%
%   Written by Zachary Barringer 10/18/2021
%
%   See also 'circular_profile.mat' by Elijah Schold for parent program


i3=data;
lx=location_x;
ly=location_y;
lz=location_z;
rr3=radius;
vec = vector/norm(vector); %normalize the input vector

%origin not allowed outside of image:
if lx>size(i3,2) || ly>size(i3,1) || lz>size(i3,3)
    error('The origin chosen is outside the image... try again!')
end

%create meshgrid of image
[xx, yy, zz]=meshgrid(1:size(i3,1),1:size(i3,2), 1:size(i3,3));

%recenter origin on point of interest, determined by lx, ly, and lz
xx=xx-lx;
yy=yy-ly;
zz=zz-lz;

%Calculate rotation about z to change coordinate frame so vector is in yz
%plane
rot_th = atan(vec(1)/vec(2));
if isnan(rot_th); rot_th=pi/2; end
%create rotation matrix to rotate about z
rot_mz = [cos(rot_th) -sin(rot_th) 0;...
          sin(rot_th) cos(rot_th) 0;...
          0 0 1];
      
%Calculat rotation about x to change coord frame so vector is along z
rot_phi = atan(vec(2)/vec(3));
if isnan(rot_phi); rot_phi=pi/2; end
%create rotation matrix to rotate about x (rotating -phi)
rot_mx = [1 0 0;...
          0 cos(rot_phi) sin(rot_phi);...
          0 -sin(rot_phi) cos(rot_phi)];      
      
xyz = [xx(:) yy(:) zz(:)];

rot_xyz = xyz*rot_mz*rot_mx;

xx_r = reshape(rot_xyz(:,2), size(xx)); 
yy_r = reshape(rot_xyz(:,1), size(xx));
zz_r = reshape(rot_xyz(:,3), size(xx));

% %convert to spherical:
% [th, phi, rr]=cart2sph(xx_r,yy_r,zz_r);

%convert to cylindrical:
[th, rr, zz_p]=cart2pol(xx_r,yy_r,zz_r);


%reset theta range to 0 to 2pi (instead of -pi to pi) with x=1, y=0 as th=0
%alternatively, you can have th=mod(th+pi()/4,2*pi()) to shift th=0 to -pi/4
th=mod(th,2*pi());


%extract radius of interest
rbin=zeros(size(rr),'logical');
rbin(round(rr)==rr3 & abs(zz_p)<1 & i3~=bkg)=true;



%% Extracting values and assigning angles

%Get the theta and phi of the binned points
th_b = th(rbin);
val_b = i3(rbin);

[~,B] = sort(th_b,'ascend');

 



results.angles = th_b(B);
results.values = val_b(B);

results.binning = rbin;
results.errors = ((radius-radius).^2+0.5)/radius;

%show what we got: feel free to comment out!(just delete one % in next line
%%{


bin_display = figure(223);
xsection = location_z;
imagesc(i3(:,:,xsection));axis image xy

%go back to image coordinates on meshgrid 
xx=xx+lx;yy=yy+ly;zz=zz+lz;
hold on;
plot(xx(rbin(:,:,xsection)),yy(rbin(:,:,xsection)),'rx','LineWidth',2);

slider = uicontrol('Style', 'slider', 'Callback', @sliderCallback);
slider.Max = size(i3,3);
slider.Min = 1;
slider.Value = xsection;
slider.SliderStep = [1 1]/(slider.Max - slider.Min);

function sliderCallback(slider, evt)
    fprintf('Slider value is: %d\n', get(slider, 'Value') );
    hold on;
    imagesc(i3(:,:,slider.Value)); axis tight;
    plot(xx(rbin(:,:,slider.Value)),yy(rbin(:,:,slider.Value)),'rx','LineWidth',2);
end
hold off;

graph = figure;
plot(results.angles*180/pi(),results.values,'+','LineWidth',2)
xlabel(['Angles(' char(176) ')'],'FontSize',16);
ylabel('circplot.values(units)','FontSize',16);
%}



end

