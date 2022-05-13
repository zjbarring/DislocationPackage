function [pixel_objective] = bruteforce3(dataset, threshold)

%% Written by A. Ulvestad, M. Menickelly, S. M. Wild


global ampmask phasedata uncertainty_width

%complexCrystal = load(dataset);
complexCrystal = (dataset);
uncertainty_width = 1;

crystal = abs(complexCrystal);

[lenx,leny,numSlides] = size(crystal);

%ampmask = (abs(crystal)>=threshold);
ampmask = shrink_wrap(crystal, threshold, 1,'gauss');
phasedata = angle(complexCrystal).*ampmask;

pixel_objective = zeros(lenx,leny,numSlides);

for i = 1:lenx
    for j = 1:leny
        for k = 1:numSlides
            if ampmask(i,j,k) ~= 0
                x = i;
                y = j;
                z = k;
                
                nhood_phases = [];

                lbx = max(1,floor(x - uncertainty_width)); ubx = min(lenx,ceil(x + uncertainty_width));
                lby = max(1,floor(y - uncertainty_width)); uby = min(leny,ceil(y + uncertainty_width));
                lbz = max(1,floor(z - uncertainty_width)); ubz = min(numSlides,ceil(z + uncertainty_width));
                for kk = lbx:ubx
                    for ll = lby:uby
                        for mm = lbz:ubz
                            if ampmask(kk,ll,mm) == 1
                                nhood_phases = [nhood_phases phasedata(kk,ll,mm)];
                            end
                        end
                    end
                end

                sorted_phases = sort(nhood_phases + pi);
                numPhases = length(nhood_phases);
                diffs = zeros(1,numPhases);
                for p = 1:(numPhases-1)
                    diffs(p) = 2*pi - (sorted_phases(p+1) - sorted_phases(p));
                end
                diffs(numPhases) = sorted_phases(numPhases) - sorted_phases(1);
                pixel_objective(i,j,k) = min(diffs);
                
                
%                 shift = pi;
%                 if any(sorted_phases>=pi)
%                     negshift = min(sorted_phases(sorted_phases>=pi));
%                     shift = 2*pi-negshift + .01; %.01 is just a small offset
%                 end
%                 pixel_objective(i,j,k) = max_nhood_offset3([i;j;k],shift);
            end % end if ampmask(i,j) ~= 0  
        end % end for k = 1:numSlides
    end % end for j = 1:leny
end % end for i = 1:lenx

end
