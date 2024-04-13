% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) using the TMM approach.
%
% by Gary P. Scavone, McGill University, 2021-2022.

clear; clf;
lossy = 0;   % 0 = lossless, 1 = traditional losses, 2 = Zwikker-Kosten; 3 = Bessel function
endType = 3; % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open

% Evaluation frequencies
fmax = 6000;          % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = eps:finc:fmax;
T = 20;   % temperature (C)

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Get geometry data
figure(1)
drawBore 'sevenCylinders'
fingering = 0;
[boreData, holeData] = sevenCylinders( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations
Zin = tmm( boreData, holeData, endType, f, lossy, T ); % ideally open end

% Plot result using rzplot script
figure(2)
rzplot( f, Zin, 1, true);
title('Seven cylinder structure.')
ylim([-100 100])
