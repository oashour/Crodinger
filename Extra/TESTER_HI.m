close all

dir = '../output';

fileP = 'param.txt';
fileID = fopen(fullfile(dir, fileP), 'r');
A = fscanf(fileID, '%f');
dt = A(1);                              % Temporal step size
Nx = A(2);                              % Number of fourier modes/spatial nodes
Tmax = A(3);                            % Maximum time to run simulation
A1 = A(4);                           % Maximum time to run simulation
q = A(5);
order = A(6);
A = fscanf(fileID, '%c');
type = A(1);
Nt = Tmax/dt;                           % Number of temporal nodes

fileR = sprintf('psi_f.bin');
fileID = fopen(fullfile(dir, fileR));
PSI = fread(fileID, [1 Nt], 'double')';
fclose(fileID);
t = 0:dt:(Tmax-dt);

plot(t, log(PSI))

hold on
[~,~,xmin,imin] = extrema(log(PSI));
plot(t(imin), xmin, 'ro')

tn = sort(t(imin));

h = zeros(1, length(tn)-1);
for i=1:length(tn)-1
    h(i) = tn(i+1)-tn(i);
end

figure
plot(h, '--')
figure
plot(abs(fft(h))/length(h))