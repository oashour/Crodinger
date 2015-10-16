% Copyright 2015, Omar Ashour.
% This sourcecode is available from <https://github.com/oashour/HighNLSE/>
%
% This file is part of HighNLSE.
% 
% HighNLSE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HighNLSE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HighNLSE.  If not, see <http://www.gnu.org/licenses/>.

function [A_0, time] = a_plot(PSI, Nx, t, num)
% FUNCTION: Plots analytical function of r, c, d, e vs actual data
% INPUT:
%       PSI: full spatiotemporal wave function
%       Nx: Number of Fourier modes/spatial nodes
%       t: time array
%       num: number of elements to plot. 5 is recommended.
% See the file recon for more info on the analytical fit.

CO = [   0         0    1.0000; %BLUE
         0    0.5000         0; %GREEN
    1.0000         0         0; %RED
    0.5843    0.3882    0.3882; %BROWN
    0.7500         0    0.7500; %VIOLET
    0.0784    0.1686    0.5490; %NAVY
    0.2500    0.2500    0.2500;];%?

% Prepare the actual data
J = 2500;                                    % Scaling down of number of points
PSI_k = abs(fft(PSI'))/Nx;                  % Absolute normalized fft
data = log(PSI_k(1:num, 1:J:end));          % Our data
time = t;
r = sqrt(3)/2;                               % Parameter
[tg, kg] = meshgrid(time, 1:num-1);            % mesh t and k for 2D calculations
c = 10.6352;
A_k = exp(r*kg.*(tg-c));
A_0 = sqrt(1 - 2*sum(A_k.^2));              % A_0

fit = log([A_0; A_k]);
% Plot
figure
hFit = plot(time, fit, '-', 'LineWidth', 1.5);
hold on
hData = plot(t(1:J:end), data, 'o', 'MarkerSize', 6); grid on;

set(hFit                      , ...
  'LineStyle'           , 'none' , ...
  'LineWidth'           , 1.5 );

set(hData                                , ...
  'LineStyle'           , '-'       , ...
  'Marker'              , 'none'          , ...
  'MarkerSize'          , 6      );

hXLabel = xlabel('t');
hYLabel = ylabel('log(|A_k|)');


for i=1:num
    set(hFit(i)                       , ...
        'Color'   , CO(i, :) );
    set(hData(i)                      , ...
        'Color'   , CO(i, :) );
end

%hLegend = legend( ...
%  [hFit, hData], ...
%  'Fit' , ...
%  'Numerical Data'      , ...
%  'location', 'SouthWest');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel         ], ...
    'FontName'   , 'Helvetica');
%set([hLegend, gca]             , ...
%    'FontSize'   , 9           );
set([hXLabel, hYLabel       ]  , ...
    'FontSize'   , 10          );

set(gca, ...
  'Box'         , 'off'        , ...
  'TickDir'     , 'out'        , ...
  'TickLength'  , [.02 .02]    , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'XGrid'       , 'on'         , ...
  'YGrid'       , 'on'         , ...
  'XColor'      , [.3 .3 .3]   , ...
  'YColor'      , [.3 .3 .3]   , ...
  'LineWidth'   , 1            , ...
  'YLim'        , [-6 1]      , ... 
  'XLim'        , [0 40]       , ...
  'XTick'       , [0:5:40]     );

get(gca, 'ColorOrder')

delete(hFit(1));
end