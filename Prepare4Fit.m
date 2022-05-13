function [data] = Prepare4Fit(data, name, clustering)
%Prepares data gathered in a circle around a dislocation for fitting
%Clusters data with points nearby and aligns them in either a decreasing or
%increasing order depending on the +/- trend of the points within the
%clusters. Adds 2pi phase to angles as needed to assemble the clusters in
%a single trend, "unwrapping" them in y by using the periodicity in x.
%Saves edited data in '[name].mat' and '[name].csv' for use in other
%DislocationFitting2.m and other programs. 

dataMatrix = [data.angles, data.values];

if(clustering==0) %So we can easily skip if the as grabbed data is periodic and doesn't need prep
    csvname = [name '.csv'];
    mname = [name '.mat'];

    writematrix(dataMatrix, csvname); %save CSV
    save(mname, 'data'); %save MATLAB variable

    figure;%Plot output for quick checks
    plot(data.angles, data.values, '.b');
    
    return %we don't need to perform clustering if data is already good, so exit
end

idx = clusterdata(dataMatrix, 2);

numClusters = max(idx);

%Preallocate varialables
clusterMeanVal = zeros(numClusters,1);
clusterMeanAng = zeros(numClusters,1);
trend = zeros(numClusters,1);
outputs.angles = [];
outputs.values = [];

for i=1:numClusters
   
    clusterMeanVal(i) = mean(data.values(idx==i));
    clusterMeanAng(i) = mean(data.angles(idx==i));
    trend(i) = (mean(data.angles(idx==i).*data.values(idx==i))-(clusterMeanAng(i)*clusterMeanVal(i)))...
                           /(mean(data.angles(idx==i).^2)-clusterMeanAng(i)^2);
    if isinf(trend(i)); trend(i)=0; end
    if isnan(trend(i)); trend(i)=0; end
                    
end

totalTrend = sum(trend);

if totalTrend>0
    [~,B] = sort(clusterMeanVal, 'ascend');
    
    for i=1:numClusters
        
        outputs.values = [outputs.values; data.values(idx==B(i))];
        
        if i~=1 && clusterMeanAng(B(i))<clusterMeanAng(B(i-1))
            outputs.angles = [outputs.angles; (data.angles(idx==B(i))+2*pi)];
        else
            outputs.angles = [outputs.angles; data.angles(idx==B(i))];
        end
    end
end
    
if totalTrend<=0
    [~,B] = sort(clusterMeanVal, 'descend');
    
    for i=1:numClusters
        
        outputs.values = [outputs.values; data.values(idx==B(i))];
        
        if i~=1 && clusterMeanAng(B(i))>clusterMeanAng(B(i-1))
            outputs.angles = [outputs.angles; (data.angles(idx==B(i))+2*pi)];
            data.angles(idx==B(i)) = data.angles(idx==B(i))+2*pi;
        else
            outputs.angles = [outputs.angles; data.angles(idx==B(i))];
        end
    end
end

%% Cluster again to break up inconsistent clusters

% data = outputs;
% dataMatrix = [data.angles, data.values];
% 
% idx = clusterdata(dataMatrix, 1);
% 
% numClusters = max(idx);
% 
% %Preallocate varialables
% clusterMeanVal = zeros(numClusters,1);
% clusterMeanAng = zeros(numClusters,1);
% trend = zeros(numClusters,1);
% outputs.angles = [];
% outputs.values = [];
% 
% for i=1:numClusters
%    
%     clusterMeanVal(i) = mean(data.values(idx==i));
%     clusterMeanAng(i) = mean(data.angles(idx==i));
%     trend(i) = (mean(data.angles(idx==i).*data.values(idx==i))-(clusterMeanAng(i)*clusterMeanVal(i)))...
%                            /(mean(data.angles(idx==i).^2)-clusterMeanAng(i)^2);
%     if isinf(trend(i)); trend(i)=0; end
%     if isnan(trend(i)); trend(i)=0; end
%                     
% end
% 
% totalTrend = sum(trend);
% 
% if totalTrend>0
%     [~,B] = sort(clusterMeanVal, 'ascend');
%     
%     for i=1:numClusters
%         
%         outputs.values = [outputs.values; data.values(idx==B(i))];
%         
%         if i~=1 && clusterMeanAng(B(i))<clusterMeanAng(B(i-1))
%             outputs.angles = [outputs.angles; (data.angles(idx==B(i))+2*pi)];
%         else
%             outputs.angles = [outputs.angles; data.angles(idx==B(i))];
%         end
%     end
% end
%     
% if totalTrend<=0
%     [~,B] = sort(clusterMeanVal, 'descend');
%     
%     for i=1:numClusters
%         
%         outputs.values = [outputs.values; data.values(idx==B(i))];
%         
%         if i~=1 && clusterMeanAng(B(i))>clusterMeanAng(B(i-1))
%             outputs.angles = [outputs.angles; (data.angles(idx==B(i))+2*pi)];
%             data.angles(idx==B(i)) = data.angles(idx==B(i))+2*pi;
%         else
%             outputs.angles = [outputs.angles; data.angles(idx==B(i))];
%         end
%     end
% end

%% Algorithm before clustering aka the Graveyard of Failure
%if length(data.values)~=length(unique(data.values)); throw('all values are not unique. Perform by hand'); end

% %Ensure there is at least 1 zero so dataset can be aligned
% data.angles = data.angles-min(data.angles);
% 
% %Check number of zeros at beginning of angle data, to use in check data
% %repetition
% numZero = sum(data.angles==0);
% %Check for repetion in the dataset; since all data from circular profile
% %and circular profile is sorted by angle on a range [0,2pi) all of the
% %0 angles should be in the first numZero positions of data.angles if the
% %data hasn't been cloned
% isRepeated = sum(data.angles(1:numZero)~=0);
% 
% origLength = length(data.angles);
% 
% 
% 
% if(~isRepeated)%check if the cloning has been done already
%     
%     data.angles = [data.angles-2*pi; data.angles; data.angles+2*pi];
%     data.values = [data.values; data.values; data.values];
%     data.errors = [data.errors; data.errors; data.errors];
%     
% end
% 
% %Determine whether the data is trending positive or negative to seperate
% %overlapping peridic boundaries
% meanValue = mean(data.values);
% meanAngle = mean(data.angles);
% 
% %formula for slope of trendline
% trendslope = sum(data.angles-meanAngle)*sum(data.values-meanValue)/sum(data.angles-meanAngle)^2;
% 
% %Round to +/-1 for
% if(trendslope>=0);trend=1;end
% if(trendslope<0);trend=-1;end
% 
% %Create offset arrays for convenience
% index = [data.values; 0];
% indexplus1 = [0; data.values];
% 
% %Find the largest jump between two values to set as the start of the loop
% [maxchange,maxind] = max(index-indexplus1);
% [minchange,minind] = min(index-indexplus1);
% 
% 
% if(abs(minchange)>abs(maxchange)) %detemine index to start at
%     loopstart=minind; 
%     maxchange=abs(minchange);
%     %if(trend==1); loopstart=loopstart-1;end %If immediate slope is opposite trend get previous point to start
%     
% else
%     loopstart = maxind; 
%     %if(trend==-1); loopstart=loopstart-1;end %If immediate slop is opposite trend get  previous point to start
% end
% 
% 
% 
% 
% %define some variables
% holdthese = data;
% clear data;
% j = 1; %indexing variable for new data structure
% 
% %Time to clean the data
% for i = 1:length(holdthese.angles) %Iterate through the cloned data set
%     
%     %Skip through determined loop start index   
%     if (i == loopstart)
%         data.values(j) = holdthese.values(i);
%         data.angles(j) = holdthese.angles(i);
%         data.errors = holdthese.errors;
%         
%           
%         %Check for overlap of the periodic boundary by comparing to
%         %limiting value: max or min depending on trend
%         if abs(data.values(j)-min(holdthese.values))> maxchange*0.3 && trend==1 %for uptrending
%            data.angles(j) = data.angles(j)+(2*pi);  %move early overlap to the end
%            holdthese.values(i-1) = min(holdthese.values); %clone min data to point before start to check overlap for point 2
%            j = j+1;
%            continue %continues ensure each point is only moved once
%         end
%         
%         if abs(data.values(j)-max(holdthese.values))> maxchange*0.3 && trend==-1 %for downtrending
%            data.angles(j) = data.angles(j)+(2*pi);  %move early overlap to the end
%            holdthese.values(i-1) = max(holdthese.values); %clone max data to point before start to check overlap for point 2
%            j = j+1;
%            continue
%         end
%         
%         j = j+1;
%         
%                 
%     elseif i>loopstart && i<loopstart+origLength %Now go until we have a length equal to input
%         
% %         plot(data.angles, data.values, '.b') %debug
%         data.values(j) = holdthese.values(i);
%         data.angles(j) = holdthese.angles(i);
%         data.errors = holdthese.errors;
% 
%         
%         %Check for overlap of the periodic boundary
%         
%         %This case finds overlapping points by comparing to previous point
%         if abs(data.values(j) - data.values(j-1)) > maxchange*0.5 && abs(data.angles(j) - data.angles(j-1)) < pi%abs(data.values(j) - holdthese.values(i-2)) > maxchange*0.5
%            if(j<origLength/2)
%                data.angles(j) = data.angles(j)+(2*pi);  %move early overlap to the end
%                %holdthese.angles(i) = data.angles(j)+(2*pi);
%                j = j+1;
%                continue %continues to ensure each point is only moved once
%            elseif(j>=origLength/2)
%                data.angles(j) = data.angles(j)-(2*pi);  %move late overlap to the beginning
%                %holdthese.angles(i) = data.angles(j)-(2*pi);
%                j = j+1;
%                continue
%            end
%         end
%         
%         %This case finds overlapping points by comparing closest value's
%         %angle        
%         [~, minindex] = min(abs(data.values(1:j-1)-data.values(j))); 
%         while abs(data.angles(j)-data.angles(minindex))>pi &&...       %if the nearest value is farther than pi
%                 abs(data.values(minindex)-data.values(j)) < maxchange*0.5   %and the values are close
%             
%             if data.angles(j)>data.angles(minindex) %move big angles down
%                 data.angles(j)=data.angles(j)-2*pi; 
%                 %j = j+1; continue; 
%                 if abs(data.angles(j)-data.angles(minindex))<pi; break; end
%             end 
%             
%             if data.angles(j)<data.angles(minindex) %move small angles up
%                 data.angles(j)=data.angles(j)+2*pi; 
%                 %j = j+1; continue; 
%                 if abs(data.angles(j)-data.angles(minindex))<pi; break; end
%             end 
%         end
%         
%         if abs(data.values(j) - data.values(j-1)) > maxchange*0.5 && abs(data.angles(j) - data.angles(j-1)) > pi
%             %Compare slope to trendline until correct
%             while (data.values(j)-data.values(j-1))/(data.angles(j)-data.angles(j-1)) < 0 && trend==1
%                 data.angles(j) = data.angles(j)+2*pi;
%             end
%             while (data.values(j)-data.values(j-1))/(data.angles(j)-data.angles(j-1)) > 0 && trend==-1
%                 data.angles(j) = data.angles(j)+2*pi;
%             end
%         
%         end
%         
%         j = j+1;
%         
% 
%         
%            
%         
%         
%     end
%     
% end
% 
% if min(data.angles(2:origLength))-data.angles(1) > 2*pi
%     while min(data.angles(2:origLength))-data.angles(1) > 2*pi; data.angles(1) = data.angles(1)+2*pi; end
% elseif max(data.angles(2:origLength))-data.angles(1) < -2*pi
%     while max(data.angles(2:origLength))-data.angles(1) < -2*pi; data.angles(1) = data.angles(1)-2*pi; end
% end
%     
% while(min(data.angles)>0 || max(data.angles)<0)%check that the data is centered around 0
%     
%     switch(sum(data.angles)>0)
%         
%         case true 
%             
%             data.angles = data.angles-2*pi;
%             
%         case false
%             
%             data.angles = data.angles+2*pi;
%         
%     end
%     
% end

 %% Export results
 
data = outputs;
export_matrix = [outputs.angles; outputs.values]; %Data in matrix form for csv

csvname = [name '.csv'];
mname = [name '.mat'];

writematrix(export_matrix, csvname);
save(mname, 'data');

figure;
plot(outputs.angles, outputs.values, '.b');
end
