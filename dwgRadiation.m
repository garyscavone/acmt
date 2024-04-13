function [b, a] = dwgRadiation(ra, fs, nb, na, T, type, doplot)
% OPENPIPE: Cylindrical pipe open-end reflectance filter least squares fit
%           to either the unflanged solution by Levine & Schwinger (1948),
%           the unflanged approximation of Dalmont et al. (2001), an
%           unflanged approximation by Causse (1984) or the flanged
%           solution of Norris and Sheng (1989). The theoretical results
%           are only strictly valid out to k*a = 3.83 though this function
%           will design a filter out to fs / 2.
%
% [B, A] = DWGRADIATION(RA, FS, NB, NA, T, TYPE, DOPLOT) gives real
% numerator and denominator coefficients B and A of orders NB and NA
% respectively, where RA is the radius of the cylinder (in meters), FS is
% the desired sampling rate (in Hz), T is an optional air temperature in
% degrees Celsius (default = 20 C), TYPE is an optional parameter
% specifying a particular condition or formula, and DOPLOT is an optional
% parameter to plot the continuous-time and discrete-time reflectance
% magnitude and phase responses (0 = no plotting (default), otherwise
% plotting). For TYPE, the default ('dalmont') is an unflanged
% approximation provided in [1]. Other options include the unflanged
% ('unflanged') solution by Levine & Schwinger (1948), an unflanged
% approximation by Causse ('causse') and the flanged ('flanged') solution
% of Norris and Sheng (1989). This function calls INVFREQZ with an f^(-2)
% weighting function for the filter design process.  This function calls
% radiation.m to get the appropriate continuous-time frequency response.
%
% 1. J. Dalmont and C.J. Nederveen, "Radiation impedance of tubes with
%      different flanges: numerical and experimental investigations,"
%      Journal of Sound and Vibration,  Vol. 244, pp. 505-534, 2001.
%
% 2. R. Caussé, J. Kergomard, and X. Lurton, "Input impedance of brass
%      musical instruments - Comparaison between experimental and numerical
%      models," J. Acoust. Soc. Am.,  Vol. 75, pp. 241-254, 1984.
%
% 3. H. Levine and J. Schwinger, "On the radiation of sound from an
%      unflanged circular pipe," Phys. Rev., 73(4), pp. 383-406, 1948.
%
% 4. A.N. Norris and I.C. Sheng. "Acoustic radiation from a circular pipe
%      with an infinite flange." Journal of Sound and Vibration, Vol. 135,
%      pp. 85-93, 1989.
%
% By Gary P. Scavone, McGill University, 2021-2024.

if nargin < 4 || nargin > 7
  error('dwgRadiation: Number of arguments is incorrect.');
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'type', 'var')
  type = 'dalmont';
end
if ~exist( 'doplot', 'var')
  doplot = 0;
end
if sum( type == ["unflanged"; "dalmont"; "causse"; "flanged"] ) == 0
  error('dwgRadiation: Type argument unrecognized.');
end

% Physical constants and evaluation frequencies
n = 1024;           % number of evaluation frequencies
omega = pi/n:pi/n:pi;

[~, R] = radiation( ra, omega*fs/(2*pi), T, type);

% Design digital filter using method described in Julius O. Smith thesis,
% pp. 47-50 & 101-103.

% Weighting to help fit at low frequencies
wt = 1./omega.^2;
[b, a] = invfreqz(R, omega, nb, na, wt, 30); % iter parameter forces stability

% Check filter stability
ps = roots(a);
if ~isempty(ps)
  for n = 1:length(ps)
    if abs(ps(n)) >= 1.0
      disp('dwgRadiation: Filter is unstable ... change design parameters!');
    end
  end
end

% Plot responses
if doplot
  h = freqz(b, a, omega);
  clf
  subplot(2,1,1)
  plot(omega/pi, 20*log10(abs([R; h].')));
  legend('Continuous-Time Response','Discrete-Time Response');
  title(['Open Pipe Reflectance (type = ', type, ', order = ', num2str(max([na nb])), ')']);
  ylabel('Log Magnitude (dB)');
  xlabel('Normalized Discrete-Time Radian Frequency')
  grid
  subplot(2,1,2)
  plot(omega/pi, angle([R; h].'));
  ylabel('Phase (radians)');
  xlabel('Normalized Discrete-Time Radian Frequency')
  grid
  disp('Paused ... hit any key to continue.')
  pause
end
