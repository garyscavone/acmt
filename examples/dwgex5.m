% Matlab script to test the digital waveguide class implementation.
%
% Gary Scavone, McGill University, 2022-2024.

clear; clf
N = 4000;  % number of samples to compute
fs = 48000; % sample rate
T = 24;     % temperature

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Load a geometry file
fingering = 0;
[ boreData, holeData ] = keefeFlute();

% Create the digital waveguide class instance and add cylindrical segment
mydwg = dwg( fs, T );
mydwg.setDefaults(struct('fracType', 'lagrange', 'fracOrder', 5, 'lossType', 'shelf', ...
  'lossOrder', 5, 'toneholeType', 'wdf'));
mydwg.setGeometry( boreData, holeData );
mydwg.setLossFlag( 1 );      % 0 = no losses, 1 = loss filtering
mydwg.setFracDelayFlag( 1 ); % 0 = no fractional delay, 1 = fractional delay

% Specify an unflanged open end modeled by a 2-zero, 1-pole digital filter
mydwg.setOutputEnd( 2, 2, 1 ); % 0 = closed; 1 = ideally open; 2 = open unflanged, 3 = open flanged
mydwg.setInputEnd( 1 );   % 0 = closed input; 1 = anechoic input (reflectance)
%mydwg.drawShape();       % draw a 3D representation of the geometry
%disp('paused ... hit key to continue'); pause

x = [1, zeros(1, N-1)];     % pressure traveling-wave impulse
r = mydwg.processInput(x);  % reflection function (anechoic input)
t = (0:N-1) * 1000 / fs;
subplot(2, 1, 1)
plot(t, r, 'b');
xlim([0 15]);

% Compute reflectance and then input impedance
R = fft( r );
M = (N / 2) + 1;
R = R(1:M);
Zin = (1 + R) ./ (1 - R); % convert reflectance to impedance for cylindrical input
f = fs*(0:M-1)/(M-1)/2;
subplot(2, 1, 2)
plot(f, 20*log10(abs(Zin(1:M))), 'b');

plotTypes = [11 1];

% Compute TMM result for comparison
lossType = 1;  % 0 = lossless, 1 = traditional losses, 2 = Zwikker-Kosten; 3 = Bessel function
endType = 1;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open
f(1) = eps;
Zin = tmm( boreData, holeData, endType, f, lossType, T );
rzplot( f, Zin, plotTypes, true, true, [], 'r', true );
ylim( [-50 50] );
subplot(2, 1, 1)
legend('DWG', 'TMM')
title( 'Impulse response and input impedance of Keefe flute (lossy, unflanged Z_L)')
