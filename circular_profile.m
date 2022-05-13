function [circplot]=circular_profile(image,location_x,location_y,radius)
%plots circular line profile about point given by location_x and location_y
%with radius given in pixels. Uses angles in mathmatically positive
%direction from 0 to 2pi, with x,y=1,0 corresponding to 0 angle. Have fun!

%written by Elijah Schold, RPI, 4-1-20

i3=image;
lx=location_x;
ly=location_y;
rr3=radius;

%origin not allowed outside of image:
if lx>size(i3,2) || ly>size(i3,1)
    error('The origin chosen is outside the image... try again!')
end

%create meshgrid of image
[xx, yy]=meshgrid(1:size(i3,2),1:size(i3,1));

%recenter origin on point of interest, determined by lx, ly
xx=xx-lx;
yy=yy-ly;

%convert to polar:
[th, rr]=cart2pol(xx,yy);


%reset range to 0 to 2pi (instead of -pi to pi) with x=1, y=0 as th=0
%alternatively, you can have th=mod(th+pi()/4,2*pi()) to shift th=0 to -pi/4
th=mod(th,2*pi());

%extract radius of interest
rbin=zeros(size(rr),'logical');
rbin(round(rr)==rr3)=true;%can try ceil or floor(rr) instead of round, but
                          %round makes the nicest circle
%% What we really came here for:

%get theta circplot.values
aa=th(rbin);
[~, B]=sort(aa);
circplot.angles=aa(B);

%equivelantly get intensity circplot.values:
ii3=i3(rbin);
circplot.values=ii3(B);

%while we're here, may as well grab radii as well: (can make additional
%output variable if you want them
rr3=rr(rbin);
radii=rr3(B);

%show what we got: feel free to comment out!(just delete one % in next line
%%{
figure(223);
imagesc(i3);axis image xy

%go back to image coordinates on meshgrid 
xx=xx+lx;yy=yy+ly;
hold on;
plot(xx(rbin),yy(rbin),'rx','LineWidth',2);
hold off;

figure(777777);
plot(circplot.angles*180/pi(),circplot.values,'+','LineWidth',2)
xlabel(['Angles(' char(176) ')'],'FontSize',16);
ylabel('circplot.values(units)','FontSize',16);
%}


circplot.errors = ((radii-radius).^2+0.5)/radius;






end