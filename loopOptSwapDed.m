function loopOptSwapDed

kos = [0];
swaps = [1,2,3,4,5];

for i=1:length(kos)
    for j=1:length(swaps)
        runOptSwapDed(kos(i), swaps(j));
    end
end

end