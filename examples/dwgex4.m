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
boreData = conicSection();

% Specify whether the input is closed (use 0) to compute an impulse
% response or anechoically terminated (use 1) to compute a reflection
% function.
inputEndType = 0;

% Create the digital waveguide class instance and add cylindrical segment
mydwg = dwg( fs, T );
mydwg.setDefaults(struct('fracType', 'lagrange', 'fracOrder', 5, 'lossType', 'shelf', ...
  'lossOrder', 5, 'toneholeType', 'twoport'));
mydwg.setGeometry( boreData );
mydwg.setLossFlag( 1 );      % 0 = no losses, 1 = loss filtering
mydwg.setFracDelayFlag( 1 ); % 0 = no fractional delay, 1 = fractional delay

% Specify an unflanged open end modeled by a 2-zero, 2-pole digital filter
mydwg.setOutputEnd( 2, 2, 2 ); % 0 = closed; 1 = ideally open; 2 = open unflanged, 3 = open flanged
mydwg.setInputEnd( inputEndType );  % 0 = closed input; 1 = anechoic input (reflectance)
%mydwg.drawShape();       % draw a 3D representation of the geometry
%disp('paused ... hit key to continue'); pause

% Note that we use a unit impulse pressure traveling-wave input signal but
% for an input conic segment with a closed input end, this will be passed
% through a spherical-wave characteristic impedance filter to convert it to
% a corresponding volume flow traveling-wave signal.
x = [1, zeros(1, N-1)];    % pressure traveling-wave unit impulse vector
p = mydwg.processInput(x);

M = (N / 2) + 1;
f = fs*(0:M-1)/(M-1)/2;
f(1) = eps;
t = (0:N-1) * 1000 / fs;

[c, rho] = thermoConstants( T );
Z0 = rho * c / (pi * boreData(2, 1)^2);
x0 = boreData(2, 1) * diff(boreData(1, 1:2)) / diff(boreData(2, 1:2));
jkx0 = 1i * 2 * pi * f * x0 / c;
Zc = Z0 * ( jkx0 ./ ( jkx0 + 1 )); % spherical-wave characteristic impedance in +x direction

if inputEndType == 0  % closed ... impulse response
  Zin = fft( p );     % compute input impedance from impulse response
  Zin = Zin(1:M);
else % anechoic ... reflection function
  Rs = fft( p );         % compute spherical reflectance from reflection function
  Rs = Rs(1:M);
  Zin = Zc .* (Rs + 1) ./ (1 - Rs.*Zc./conj(Zc)) / Z0;
end

% Choose the signals to plot: % 1 = input impedance magnitude, 6 =
% reflectance magnitude; 10 = impulse response; 11 = reflection function
plotTypes = [10 1];
rzplot( f, Zin, plotTypes, true, false, [], 'b', 0, Z0, Zc );

% Compute TMM result for comparison
lossType = 1;  % 0 = lossless, 1 = traditional losses, 2 = Zwikker-Kosten; 3 = Bessel function
endType = 1;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open
Zin = tmm( boreData, [], endType, f, lossType, T );
rzplot( f, Zin, plotTypes, true, true, [], 'r', 0, Z0, Zc );
if sum(plotTypes == 1) % input impedance magnitude
  ylim( [-50 50] );
end
subplot(2, 1, 1)
legend('DWG', 'TMM')
title( 'Impulse response and input impedance of conic segment (unflanged Z_L)')
