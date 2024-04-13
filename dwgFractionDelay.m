function [b, a, m] = dwgFractionDelay(nZ, type, doplot, maxN)
% DWGFRACTIONDELAY: Discrete filter approximation for fractional delay
%                   length.
%
% [B, A, M] = DWGFRACTIONDELAY(NZ, TYPE, DOPLOT, MAXN) computes real numerator
% and denominator coefficients B and A for a digital filter that
% approximates the fractional delay length of NZ - floor(NZ) + M, where NZ
% is the total delay length available (including the fractional part) and M
% is a number of unit delays that are incorporated into the filter design
% (which should be subtracted from NZ before implementing the integer delay
% length). TYPE can be either 'lagrange' (default) or 'thiran'. If the
% optional value of DOPLOT > 0 (default = 0), the continuous-time and
% discrete-time magnitude and phase responses will be plotted. MAXN
% specifies the maximum filter order (default = 5 for Lagrange and 4 for
% Thiran). Thiran filters are allpass, while Lagrange filters are FIR. Only
% odd-order Lagrange filters are designed.
% 
% References:
%
% - https://ccrma.stanford.edu/~jos/pasp/Matlab_Code_Lagrange_Interpolation.html
% 
% - https://ccrma.stanford.edu/~jos/pasp/Thiran_Allpass_Interpolation_Matlab.html
%
% Valimaki and Laakso, "Principles of Fractional Delay Filters," ICASSP
% 2000.
%
% By Gary P. Scavone, McGill University, 2022-2024.

if nargin < 1
  error('dwgFractionDelay: At least one argument is necessary.');
end

if nargin < 2, type = 'lagrange'; end
if nargin < 3, doplot = 0; end
if nargin < 4, maxN = 5; end

if rem(maxN, 1)
  error('dwgFractionDelay: maxN argument must be an integer.');
end

delta = nZ - floor(nZ);
D = delta;
N = 1;
m = 0;

if strcmp(type, 'lagrange')

  if ~rem(maxN, 2), maxN = maxN - 1; end % force odd
  for n = maxN:-2:1
    m = floor(n/2);
    if nZ > m
      D = delta + m;
      N = n;
      break;
    end
  end
  b = lagrange(N, D);
  a = 1;

elseif strcmp(type, 'thiran')
  % Thiran: stable for D > N - 1, optimal range is N-0.5 to N+0.5
  % - if D between 0 - 1.5, N=1
  % - if D between 1.5 - 2.5, N=2
  % - if D between 2.5 - 3.5, N=3
  % - if D between 3.5 - 4.5, N=4
  for n = maxN:-1:1
    if nZ > n + 1
      m = n;
      D = delta + m;
      N = round(D);
      break;
    end
  end
  [b, a] = thiran(N, D);
  if sum(isnan(a)) || sum(abs(roots(a)) > 1)
    fprintf(1, 'Unstable filter for delay = %f and order = %f\n', D, N);
  end
else
  error('dwgFractionDelay: Type argument unrecognized.');
end

if doplot
  [h, w] = freqz(b, a);
  subplot(2, 1, 1)
  plot(w/pi, 20*log10(abs(h)));
  ylabel('Magnitude (dB)')
  title(['Fractional Delay Response (type = ', type, ', order = ', num2str(N), ', delay = ', num2str(delta), ')']);
  subplot(2, 1, 2)
  plot(w/pi, -unwrap(angle(h)) ./ w)
  ylabel('Phase Delay (samples)')
  xlabel('Normalized Frequency')
  disp('Paused ... hit any key to continue.')
  pause
end
