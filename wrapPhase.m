function wrapped = wrapPhase(data, wrap)

shift = 0;
if min(data)<0 
    data = data+pi;
    shift = 1;
end

data = data+wrap;
data = mod(data, 2*pi);

if shift; wrapped = data-pi; 
else; wrapped = data; end

end