function  [output, cuboid] = cropCDI(input, padding)
% Remove Slices of a 3D array that are 0 or NaN, leaving padding in each
% dimention. Data is returned in output, the cuboid is return so the same
% cropping can be applied elsewhere.


input(input==0) = NaN;
input(input==input(1,1,1)) = NaN;

[x, y, z] = ind2sub(size(input),find(~isnan(input)));

cuboid = [min(y)-padding, min(x)-padding, min(z)-padding,...
            max(y)-min(y)+2*padding, max(x)-min(x)+2*padding, max(z)-min(z)+2*padding];

output = imcrop3(input, cuboid);


end