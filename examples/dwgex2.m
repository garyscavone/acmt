% Matlab script to test the digital waveguide class implementation.
%
% Gary Scavone, McGill University, 2022-2024.

clear; clf
N = 4000;  % number of samples to compute
fs = 48000; % sample rate
T = 24;     % temperature

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Load a geometry file for a 0.6 m pipe of radius 0.0075 m
boreData = pipe();

% Create the digital waveguide class instance and add cylindrical segment
mydwg = dwg( fs, T );
%mydwg.setDebugFlag( true );
%mydwg.setPlotFlag( true );
mydwg.setGeometry( boreData );

mydwg.setLossFlag( 0 );      % 0 = no losses, 1 = loss filtering
mydwg.setFracDelayFlag( 0 ); % 0 = no fractional delay, 1 = fractional delay

% Specify an unflanged open end modeled by a 2-zero, 1-pole digital filter
mydwg.setOutputEnd( 2, 2, 1 ); % 0 = closed; 1 = ideally open; 2 = open unflanged, 3 = open flanged
mydwg.setInputEnd( 1 );  % 0 = closed input; 1 = anechoic input (reflectance)

x = [1, zeros(1, N-1)];    % pressure traveling-wave unit impulse vector
p = mydwg.processInput(x); % compute outputs for length of input vector

R = fft( p );              % compute reflectance from reflection function
Zin = (1 + R) ./ (1 - R);  % convert reflectance to impedance for cylindrical input
M = (N / 2) + 1;
f = fs*(0:M-1)/(M-1)/2;
t = (0:N-1) * 1000 / fs;

plotTypes = [11 1];
subplot(2, 1, 1)
plot(t, p, 'b');
subplot(2, 1, 2)
plot(f, 20*log10(abs(Zin(1:M))), 'b');
% rzplot( f, Zin(1:M), plotTypes, true, false, [], 'b' ); use rzplot instead

% Compute TMM result for comparison
lossType = 0;  % 0 = lossless, 1 = traditional losses, 2 = Zwikker-Kosten; 3 = Bessel function
endType = 1;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open
Zin = tmm( boreData, [], endType, f, lossType, T );
rzplot( f, Zin, plotTypes, true, true, [], 'r' );
ylim( [-50 50] );
subplot(2, 1, 1)
legend('DWG', 'TMM')
title( 'Reflection function and input impedance of 60 cm pipe (unflanged Z_L)')
