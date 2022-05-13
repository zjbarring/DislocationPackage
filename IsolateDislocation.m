function [Sep_Dislocations] = IsolateDislocation(Crystal,NumUnique)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Sep_Dislocations = struct;

[xCoord, yCoord, zCoord] = size(Crystal);

Threshold = 4; %Threshold fo consider a point part of a dislocation. Inherited from previous display
Tolerance = 3; %How far appart can two points be and still be a part of the same dislocation

for n=1:NumUnique
    
    DislocationFlag = 0; %set flag to track whether we have found a uniques dislocation to track
    Sep_Dislocations(n).dislocation = zeros(size(Crystal)); %define starting values for active dislocation
    
    for i=1:xCoord
        for j=1:yCoord
            for k=1:zCoord
                
                if(Crystal(i,j,k) > Threshold)%value is big enough to be part of dislocation
                    
                    %If there is no active dislocation
                    if(DislocationFlag == 0)
                        Sep_Dislocations(n).dislocation(i,j,k) = Crystal(i,j,k);
                        Crystal(i,j,k) = 1234; %Value of 1234 indicates the index is a part of the current dislocaiton
                        DislocationFlag = 1;
                        
                    
                    %If the point is adjacent to the active dislocation    
                    elseif max(Crystal((i-Tolerance):(i+Tolerance),...
                            (j-Tolerance):(j+Tolerance), (k-Tolerance):(k+Tolerance)),[], 'all') == 1234
                            Sep_Dislocations(n).dislocation(i,j,k) = Crystal(i,j,k);
                            Crystal(i,j,k) = 1234;
                                  
                    end%end DislocationFlag if
                        
                end%end for thresholding
                
                
            end%end for zCoord
        end%end for yCoord
    end%end for zCoord
        
            
    %after cycling through all the points, set all points flagged for the active dislocation to zero
    Crystal = Crystal.*(Crystal~=1234);
    

    
end%end for NumUnique

end

