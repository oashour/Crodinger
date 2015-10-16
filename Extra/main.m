clear all
close all

dir = 'output';

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

lambda = 2*sqrt(2*q*(1-2*q));
Omega = 2*sqrt(1-2*q);

fileR = sprintf('psi_r.bin');
fileI = sprintf('psi_i.bin');
fileX = sprintf('x.bin');
fileID = fopen(fullfile(dir, fileR));
psi_r = fread(fileID, [Nx Nt], 'double')';
fclose(fileID);
fileID = fopen(fullfile(dir, fileI));
psi_i = fread(fileID, [Nx Nt], 'double')';
fclose(fileID);
fileID = fopen(fullfile(dir, fileX));
x = fread(fileID, 'double');
fclose(fileID);
t = 0:dt:(Tmax-dt);
PSI = psi_r + 1i*psi_i;
Lx = pi/sqrt(1-2*q);
k = 2*(-Nx/2:1:Nx/2-1)'*pi/Lx;          % Wave number

% % Plot results 
%densityPlot(PSI, x, t, dt, 200);
 
maxima = regions(PSI, x, t, q);
%a_plot(PSI, Nx, t, 6);
%b_plot(PSI, Nx, t, 6, maxima(1,1), q);
%c_plot(PSI, Nx, t, 6, maxima(1,1), q);
% %A = recon(Nx, 5000, max(t), maxima(1,1));
% [shift] = ab(PSI, x, t, q);      
% energy(PSI, t, k2, Nx, gamma, dt);

% %T_F = maxima(1,1);
% %T_R = maxima(2,1)-T_F;
% %if ~(isempty(maxima))
% %    fid = fopen(fileRes, 'a+');
% %    fprintf(fid, '%2.1f    %2.16f    %2.16f\n', w, T_F, T_R);
% %    fclose(fid);
% %end

t0 = -log(A1/Omega)/lambda + log(lambda)/lambda;
tf = maxima(1,1);

fid = fopen(fullfile(dir, 'results_q.txt'), 'a+');
fprintf(fid, '%2.2f     %2.9f    %2.16f\n', A1, q, tf);
