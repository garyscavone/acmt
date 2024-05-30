function [br, ar, bt, at] = dwgTonehole( delta, rb, height, state, fs, T, chimney, rpad, hpad, w, doplot )
% DWGTONEHOLE:  Determine digital filter coefficients for tonehole
%               junction reflectance and transmittance.
%
%   [BR, AR, BT, AT] = DWGTONEHOLE( DELTA, RB, HEIGHT, STATE, FS, T, ...
%                                   CHIMNEY, RPAD, HPAD, W, DOPLOT)
%   gives real numerator and denominator reflectance coefficients BR and AR
%   and transmittance coefficients BT and AT, where DELTA is the ratio of
%   tonehole to main bore radii, RB is the tonehole radius, HEIGHT is the
%   tonehole wall height, T is the temperature (celsius), STATE is an
%   openness value between 0 (closed) and 1 (fully open), FS is the model
%   sampling rate (in Hz), CHIMNEY is the length the tonehole extends out
%   from the bore, RPAD is the radius of a hanging pad over the hole (use 0
%   for no pad), HPAD is the distance of the pad from the hole and W is the
%   wall thickness of the tonehole (default values = 0 if not otherwise
%   specified). All length values should be given in meters. The
%   discrete-time filters are designed using INVFREQZ with an omega^(-2)
%   weighting function. The discrete-time open hole reflectance filter is
%   designed by Kopec's method, in conjunction with INVFREQZ, and is a
%   two-zero/one-pole filter. All other filters are two-zero/two-pole. If
%   the value of DOPLOT is 1, the continuous-time and discrete-time
%   immitance magnitude and phase responses will be plotted.
%
%   By Gary P. Scavone, CAML, McGill University, 2022-2024.

if nargin < 5 || nargin > 11
  error( 'dwgTonehole: Incorrect number of parameters.' );
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'chimney', 'var')
  chimney = 0;
end
if ~exist( 'rPad', 'var')
  rPad = 0;
end
if ~exist( 'hPad', 'var')
  hPad = 0;
end
if ~exist( 'w', 'var')
  w = 0;
end
if ~exist( 'doplot', 'var')
  doplot = 0;
end

N = 4096;
f = fs * (1:N)/N/2;
lossType = 1; % standard tmm first-order losses
type = 'Lefebvre2012'; %'Keefe1990'

Gamma = sectionLosses( rb, rb, 0, f, T, lossType );
[A, B, C, D] = tmmTonehole( delta, rb, height, state, Gamma, type, T, chimney, rPad, hPad, w );

[c, rho] = thermoConstants( T );
ra = rb / delta;
Zo = rho * c / (pi * ra^2 );
Zo2 = Zo^2;
den = (A*Zo + B + C*Zo2 + D*Zo);
R = (A*Zo + B - C*Zo2 - D*Zo) ./ den;
T = 2 * (A.*D*Zo - B.*C*Zo) ./ den;

% Design digital filter using method described in Julius O. Smith thesis,
% pp. 47-50 & 101-103 and implemented in invfreqz() function.

omegaof = 2 * pi * f / fs; % normalized omega

% Weighting to help fit at low frequencies
wt = (omegaof).^(-2);

% Numerator and denominator filter orders
nb = 2;
na = 2;
iter = 40;

if state
  [b1, a1] = invfreqz(R, omegaof, 0, 1, wt, iter);
  temp = freqz(b1, a1, omegaof);
  [b2, a2] = invfreqz(temp./R, omegaof, 0, 2, wt, iter);
  br = b1 * a2 / b2;
  ar = a1;
  br = -abs(R(1))*sum(ar)*br/sum(br);
else
  [br, ar] = invfreqz(R, omegaof, nb, na, wt, iter);
end
[bt, at] = invfreqz(T, omegaof, nb, na, wt, iter);

% Check filter stability. Using an iteration limit (last parameter in
% invfreqz() function above) should guarantee a stable result but the
% value we are using is arbitrary.
psr = roots(ar);
pst = roots(at);
if sum( abs(psr) >= 1.0) > 0 || sum( abs(pst) >= 1.0 ) > 0
  disp('dwgTonehole: Filter is unstable ... change design parameters!');
end

if doplot
  hr = freqz(br, ar, omegaof);
  ht = freqz(bt, at, omegaof);
  clf
  subplot(2, 1, 1)
  plot(fs*omegaof/(2000*pi), 20*log10(abs([R; hr; T; ht])).')
  legend('Continuous-Time Response','Discrete-Time Response');
  title(['Tonehole Immitances (type = ', type, ')']);
  ylabel('Magnitude (dB)');
  xlabel('Frequency (kHz)')
  grid
  xlim([0 fs/2000])
  subplot(2, 1, 2)
  plot(fs*omegaof/(2000*pi), unwrap(angle([R; hr; T; ht].')))
  ylabel('Phase (degrees)');
  xlabel('Frequency (kHz)')
  grid
  xlim([0 fs/2000])
  disp('Paused ... hit any key to continue.')
  pause
end
