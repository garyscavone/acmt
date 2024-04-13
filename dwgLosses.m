function [b, a] = dwgLosses(r1, r2, L, fs, nz, nb, na, T, type, doplot)
% DWGLOSSES: Discrete IIR filter approximation to the attenuation and phase
%            characteristic of viscous and thermal losses in a rigid duct.
%
% [B, A] = DWGLOSSES(R1, R2, L, FS, NZ, NB, NA, T, TYPE, DOPLOT) computes
% real numerator and denominator coefficients B and A for an IIR filter
% that approximates the losses in a duct of input radius R1, output radius
% R2 and length L (all lengths are in meters). FS is the system sampling
% rate (in Hz) and NZ is the number of unit sample delays for which the
% losses are modeled (which can be fractional). NB and NA are optional,
% with default values of NB=3 and NA=2 for the INVFREQZ fit and NB=5 for
% the SHELF fit. If the temperature T is not specified, a reference
% temperature of 20 degrees celsius is assumed. TYPE can be either
% 'invfreqz' (default, least squares) or 'shelf'. For the 'invfreqz' fit,
% this function calls INVFREQZ with an omega^(-2) weighting function for
% the filter design process. For the 'shelf' fit, a cascade of NB
% first-order shelf filters is computed as described in Abel, Smyth and
% Smith, 2003). If the optional value of DOPLOT != 0 (default = 0), the
% continuous-time and discrete-time reflectance magnitude and phase
% responses will be plotted.
%
% References:
%
% Scavone, "An Acoustic Analysis of Single-Reed Woodwind Instruments with
% an Emphasis on Design and Performance Issues and Digital Waveguide
% Modeling Techniques, Ph.D. thesis, Stanford University, 1997.
%
% Abel, Smyth and Smith, "A Simple, Accurate Wall Loss Filter for
% Acoustic Tubes," DAFx 2003.
%
% By Gary P. Scavone, CCRMA, Stanford University, March 1997 - 2024.
% Updated in 2021 for Abel, Smyth and Smith shelf method.

if nargin < 5 || nargin > 10
  error('dwgLosses: Number of arguments is incorrect.');
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'type', 'var')
  type = 'invfreqz';
end
if ~exist( 'doplot', 'var')
  doplot = 0;
end

% Physical constants and evaluation frequencies
N = 1024;                % Number of evaluation frequencies
omegaof = pi/N:pi/N:pi;  % normalized omega
f = omegaof*fs/(2*pi);
c = thermoConstants( T );

% Note that L is the physical duct length but nz may represent the number
% of delays for propagation over two times L (which can be more efficient
% to implement), so we don't use L in the calculation of H. But L is used
% in sectionLosses() to determine an appropriate equivalent radius for
% conical segments.
[G, ~] = sectionLosses( r1, r2, L, f, 20, 1 );
H = exp(-G*nz*c/fs + 1i*omegaof*nz);

if strcmp(type, 'invfreqz')
  
  % Design digital filter using method described in JOS thesis,
  % pp. 47-50 & 101-103 and implemented in invfreqz() function.

  % Weighting to help fit at low frequencies
  wt = (omegaof).^(-2);
  
  % Using an iteration limit should guarantee a stable result but the
  % value used here is arbitrary.
  if nargin < 4
    nb = 3;
    na = 2;
  end
  [b, a] = invfreqz(H, omegaof, nb, na, wt, 30);

  % Check filter stability
  ps = roots(a);
  if sum( abs(ps) >= 1.0) > 0
    disp('dwgLosses: Filter is unstable ... change design parameters!');
  end
  
elseif strcmp(type, 'shelf')

  % Shelf-filtering approach
  if nargin < 4
    nb = 5;
  end
  M = nb;
  if M < 2
    M = 2;
  end
  gpi = zeros(1, M);
  ft = zeros(1, M);
  b0 = zeros(1, M);
  b1 = zeros(1, M);
  a1 = zeros(1, M);
  b = 1;
  a = 1;
  tmp = sum(sqrt(((1:M)-0.5)/M));
  for m = 1:M
    gpi(m) = exp(sqrt((m-0.5)/M)/tmp * log(abs(H(end))));
    ft(m) = ((m - 0.5)/M)^3;
    eta = (gpi(m) + 1)/(gpi(m) - 1);
    if eta == 1
      alpha1 = 0;
    else
      alpha1 = eta - sign(eta)*sqrt(eta^2 - 1);
    end
    rhow = sin(pi*ft(m)/2 - pi/4)/sin(pi*ft(m)/2 + pi/4);
    beta0 = (1 + gpi(m))/2 + alpha1*(1 - gpi(m))/2;
    beta1 = (1 - gpi(m))/2 + alpha1*(1 + gpi(m))/2;
    b0(m) = (beta0 + rhow*beta1)/(1 + rhow*alpha1);
    b1(m) = (beta1 + rhow*beta0)/(1 + rhow*alpha1);
    a1(m) = (rhow + alpha1)/(1 + rhow*alpha1);
    b = conv(b, [b0(m) b1(m)]);
    a = conv(a, [1 a1(m)]);
  end
  
else
  error('dwgLosses: Type argument unrecognized.');
end

% Plot responses
if doplot
  h = freqz(b, a, omegaof);
  clf
  subplot(2, 1, 1)
  plot(fs*omegaof/(2000*pi), 20*log10(abs([H; h])).')
  legend('Continuous-Time Response','Discrete-Time Response');
  title(['Boundary-Layer Losses (type = ', type, ', order = ', num2str(max([na nb])), ')']);
  ylabel('Magnitude (dB)');
  xlabel('Frequency (kHz)')
  grid
  xlim([0 fs/2000])
  subplot(2, 1, 2)
  plot(fs*omegaof/(2000*pi), -unwrap(angle([H; h].'))./omegaof.')
  ylabel('Phase Delay (samples)');
  xlabel('Frequency (kHz)')
  grid
  xlim([0 fs/2000])
  disp('Paused ... hit any key to continue.')
  pause
end