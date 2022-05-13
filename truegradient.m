function [mgradx, mgrady, mgradz] = truegradient(phaseData, range)
%truegradient Averages gradient components over several phase shifts
%   Use this to get better estimate of phase gradients in data with
%   discontinuities in the phase due to periodic data. Works by applying
%   phasae shifts and averaging the results. Only applies to -pi to pi
%   data, not from spikes due to physical defects. 
%
%   phaseData is the dataset of the phase to be gradient. range is a 2
%   element vector containing the maximum and miniumum of the range,
%   defaults to -pi and pi.

if ~exist('range', 'var')
    range = [-pi, pi];
end

phaseData = phaseData + min(range);

upperLim = max(range)-min(range);

mgradx = zeros(size(phaseData));
mgrady = zeros(size(phaseData));
mgradz = zeros(size(phaseData));

N = 4; %# of averages

for i=1:N
    
    phaseShift = mod(phaseData + pi * (i-1)/N, upperLim);
    
    [gradx, grady, gradz] = gradient(phaseShift);
    
    mgradx = mgradx + gradx.*(abs(gradx)<upperLim/2);
    mgrady = mgrady + grady.*(abs(grady)<upperLim/2);
    mgradz = mgradz + gradz.*(abs(gradz)<upperLim/2);
    
end

mgradx = mgradx/N;
mgrady = mgrady/N;
mgradz = mgradz/N;




end

