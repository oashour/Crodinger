function [suprema] = regions_f0(PSI, t)

cutoff = 0.95;

row = find(PSI >= cutoff);
results = zeros(length(row)+1, 1);

for i = 1:length(row)
    results(i) = [t(row(i)) PSI(row(i))];
end

end

