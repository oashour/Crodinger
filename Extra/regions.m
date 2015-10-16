function [suprema] = regions(PSI, x, t, q)

t_s = 0;
x_s = x;
lambda = 2*sqrt(2*q*(1-2*q));
Omega = 2*sqrt(1-2*q);
psi_s = ((1-4*q)*cosh(lambda*t_s)+sqrt(2*q)*cos(Omega*x_s)+1i*lambda*sinh(lambda*t_s))./(sqrt(2*q)*cos(Omega*x_s)-cosh(lambda*t_s));
maximum = max(abs(psi_s)).^2;
cutoff = maximum - 0.2*maximum;

[row,col] = find(abs(PSI).^2' >= cutoff);
results = zeros(length(row)+1, 3);
for i = 1:length(row)
    results(i, :) = [t(col(i)) x(row(i)) abs(PSI(col(i), row(i))).^2];
end
results(length(row)+1, :) = [9999999999999 0 0];

pivot = 1;
suprema = [];
for i = 1:length(row)+1;
    if results(i, 1) - results(pivot, 1) < 1 
        continue
    else
        region = results(pivot:i-1, :);
        [~, index] = max(region(:, 3));
        suprema = [suprema; region(index, :)];
        pivot = i;
    end
end
end

