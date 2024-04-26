% Matlab script to test the digital waveguide class implementation.
%
% Gary Scavone, McGill University, 2022-2024.

clear; clf
N = 2000;  % number of samples to compute
fs = 48000; % sample rate
T = 24;     % temperature

% Cylinder geometry
L = 0.6;    % length (m)
r = 0.0075; % radius (m)

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Create the digital waveguide class instance and add cylindrical segment
mydwg = dwg( fs, T );        % create the class for given sample rate and temperature
mydwg.addSegment( L, r );    % add the segment to the structure

mydwg.setLossFlag( 0 );      % 0 = no losses, 1 = loss filtering
mydwg.setFracDelayFlag( 0 ); % 0 = no fractional delay, 1 = fractional delay

mydwg.setOutputEnd( 1 ); % 0 = closed; 1 = ideally open; 2 = open unflanged, 3 = open flanged
mydwg.setInputEnd( 1 );  % 0 = closed input; 1 = anechoic input (reflectance)
%mydwg.drawShape();       % draw a 3D representation of the geometry
%disp('paused ... hit key to continue'); pause

x = [1, zeros(1, N-1)];    % pressure traveling-wave unit impulse vector
p = mydwg.processInput(x); % compute outputs for length of input vector

t = (0:N-1) * 1000 / fs;
plot(t, p)
grid
xlabel('Time (ms)');
ylabel('Reflection function gain')
ylim([-1.1 0.1])
title( 'Reflection function of 60 cm pipe (lossless, Z_L = 0)')
