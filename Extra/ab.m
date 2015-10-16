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

function [t_shift] = ab(psi, x, t, q)
% FUNCTION: [spatial temporal] = ab(psi, x, t, peak)
%           Compares a full time evolution of a wave function psi to an Akhmediev
%           breather using supplied value of a
% INPUT:    psi:  Matrix representing wafe function with rows as time and
%                 columns as space.
%           x:    spatial variable used to create psi
%           t:    temporal variable used to create psi
%           peak: actual peak to compare to theoretical peak
%           a:    parameter used for akhmediev breather
% OUTPUT:   spatial: spatial profile array
%           temporal: temporal profile array

% Find row and column of the rogue wave peak
[~, maxIndex] = max(psi(:));                    %
[row, col] = ind2sub(size(psi), maxIndex);      %

x_shift = x(col);                               % x-value of the peak
t_shift = t(row);                               % t-value of the peak

% Spatial profile
spatial = psi(row,:);                           % spatial profile corresponds to 
                                                % a row of psi
t_s = 0;                                        % Use t = 0 for analytical AB 

Nx = 2^10;
dx = pi/sqrt(1-2*q)/Nx;                         % Spatial step size
x_2 = (-Nx/2:1:Nx/2-1)'*dx;                     % Spatial grid points
x_s = x_2 - x_shift;                            % Shift space for analytical AB

%b = sqrt(8*a*(1-2*a));                          % Paramter
%omega = 2*sqrt(1-2*a);                          % Parameter
% psi_s = (1 + (2*(1-2*a)*cosh(b*t_s) + 1i*b*sinh(b*t_s))./ ... % Analytical
%         (sqrt(2*a)*cos(omega*x_s)-cosh(b*t_s))).*exp(1i*t_s); % Spatial profile

lambda = 2*sqrt(2*q*(1-2*q));
Omega = 2*sqrt(1-2*q);
psi_s = ((1-4*q)*cosh(lambda*t_s)+sqrt(2*q)*cos(Omega*x_s)+1i*lambda*sinh(lambda*t_s))./(sqrt(2*q)*cos(Omega*x_s)-cosh(lambda*t_s));

% Temporal profile
%temporal = psi(:,col);                          % temporal profile corresponds 
                                                % to column of psi
%t_t = t - t_shift;                              % Shift time
%x_t = 0;                                        % Use x = 0 for analytical AB
%psi_t = (1 + (2*(1-2*a)*cosh(b*t_t) + 1i*b*sinh(b*t_t))./ ... % Analytical
%        (sqrt(2*a)*cos(omega*x_t)-cosh(b*t_t))).*exp(1i*t_t); % Temporal profile

% Plot the results
figure
hAnal = plot(x_2, abs(psi_s).^2); % Analytical spatial
hold on;                             % Grid
hNum = plot(x, abs(spatial).^2);
%xlim([-min(x), max(x)]);                              % Limits
%xlabel('x'); ylabel('|\psi|^2');             % Labels
%legend('Analytical', 'Numerical', 0); 
          
% figure
% plot(t, abs(psi_t).^2, '-r', 'LineWidth', 2); % Analytical temporal
% grid on; hold on;                             % Grid
% plot(t(1:40:end), abs(temporal(1:40:end)).^2, 'b+', ...           % Actual temporal
%      'MarkerFaceColor', 'b', 'MarkerSize', 8);
% xlim([t_shift-5, t_shift+5]);                 % Limits
% xlabel('t'); ylabel('|\psi|^2|');             % Labels
% legend('Analytical', 'Calculation', 0);          

set(hNum                        , ...
  'LineStyle'           , 'none' , ...
  'Color'               , 'b'    , ...
  'Marker'              , '+'    , ...
  'MarkerSize'          , 8      );

set(hAnal                      , ...
  'LineStyle'           , '-'  , ...
  'Color'               , 'r'  , ...
  'LineWidth'           , 1.5  );


hXLabel = xlabel('x');
hYLabel = ylabel('|\psi|^2');

hLegend = legend( ...
  [hAnal, hNum], ...
  'Analytical' , ...
  'Numerical'      , ...
  'location', 'NorthWest');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel         ], ...
    'FontName'   , 'Helvetica');
set([hLegend, gca]             , ...
    'FontSize'   , 9           );
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
  'XLim'        , [min(x) max(x)]);
  %'XTickLabel'       , sprintf('%.1f |',get(gca, 'XTick')+0.5)
      
end