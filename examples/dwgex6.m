% Matlab script to test the digital waveguide class implementation with a
% pseudo-clarinet-like reed function.
%
% Gary Scavone, McGill University, 2022-2024.

clear; clf
N = 100000;  % number of samples to compute
fs = 44100;  % sample rate
T = 24;      % temperature

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Load a geometry file
boreData = pipe();

% Create the digital waveguide class instance and add cylindrical segment
mydwg = dwg( fs, T );
mydwg.setDefaults(struct('fracType', 'lagrange', 'fracOrder', 5, 'lossType', 'shelf', ...
  'lossOrder', 5, 'toneholeType', 'wdf'));
mydwg.setGeometry( boreData );
mydwg.setLossFlag( 1 );      % 0 = no losses, 1 = loss filtering
mydwg.setFracDelayFlag( 1 ); % 0 = no fractional delay, 1 = fractional delay

% Specify an unflanged open end modeled by a 2-zero, 1-pole digital filter
mydwg.setOutputEnd( 2, 2, 1 ); % 0 = closed; 1 = ideally open; 2 = open unflanged, 3 = open flanged
mydwg.setInputEnd( 1 );        % 0 = closed input; 1 = anechoic input (reflectance)

y = zeros(1, N);        % initialize output vector
pm = 0;                 % mouth pressure
pinc = 0.01;            % mouth pressure increment

for i = 1:N

	pminus = mydwg.getNextPminus();

  % Pressure-dependent reflection coefficient (a la STK)
	p_delta = 0.5*pm - pminus;
  rc = 0.8 + 0.3*p_delta;
  if rc > 1.0
    rc = 1.0;
  end
  y(i) = mydwg.processInput( 0.5*pm - p_delta*rc );
	
	if pm < 1    % increment mouth pressure
		pm = pm + pinc;
	end
end

t = (0:N-1) / fs;
plot(t, y, 'b');
xlabel('Time (s)')
ylabel('Mouthpiece pressure')
soundsc(y, fs)
title('Clarinet-like mouthpiece pressure output')
