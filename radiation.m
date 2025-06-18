function [Zr, R] = radiation( a, f, T, type, b )
% [ZR, R] = RADIATION( A, F, TYPE, T ) computes the radiation impedance ZR
%      (not normalized by Zc) and reflectance (or reflection coefficient) R
%      of a cylindrical pipe at frequencies specified in the 1D vector F
%      (Hertz) and for pipe radius A (in meters). T is an optional air
%      temperature in degrees Celsius (default = 20 C). TYPE is an optional
%      parameter specifying a particular condition or formula. The default
%      ('dalmont') is an unflanged approximation provided in [1]. Other
%      options include the unflanged ('unflanged') solution by Levine &
%      Schwinger (1948), an unflanged approximation by Causse ('causse'),
%      the flanged ('flanged') solution of Norris and Sheng (1989) and a
%      thick pipe ('thickpipe') approximation provided in [1]. The
%      parameter B (in meters), which specifies the outer wall thickness,
%      is required with the 'thickpipe' type.
%
% by Gary P. Scavone, McGill University, 2013-2025.
% Based in part on functions from WIAT by Antoine Lefebvre.
%
% References:
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
% Generally, Zr = Zc*(1+R)/(1-R), where R is the reflection coefficient for
% the open end of a pipe and Zc is the characteristic impedance.  R can
% also be expressed as:
%     R = -|R0|exp(-2j*k*a*delta),
% where delta is the frequency dependant length correction.  See
% reference [1], pg. 509.  A low-frequency approximation for an unflanged
% Zr is:
%     Zr = 0.25*ka^2 + 0.61j*ka

if nargin < 2 || nargin > 5
  error( 'radiation: Invalid number of arguments.');
end
if ~isvector(f)
  error( 'radiation: f should be a 1D vector of frequencies in Hertz.' );
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'type', 'var')
  type = 'dalmont';
end
if strcmp( type, 'thickpipe' ) && nargin < 5
  error( 'radiation: b is a required parameter for thickpipe type.' );
end

[c, rho] = thermoConstants( T );
ka = 2 * pi * f * a / c;
ka2 = ka.^2;

if strcmp( type, 'unflanged' )
  % The Levine & Schwinger results are calculated by numerical
  % integrations. Eq. (VI.5) of Levine and Schwinger is used to calculate
  % the reflection coefficient magnitude.
  intvi5 = @(x, z) atan(besselk(1, x)./(pi*besseli(1, x))) ...
    .*(1 - z./sqrt(z^2 + x.^2)).*(1./x);

  if ka(1) == 0, ka(1) = eps; end % avoid division by 0
  sum1 = zeros(size(ka));
  for n=1:length(ka)
    sum1(n) = integral(@(x)intvi5(x, ka(n)), 0, 20); % upper limit determined empirically
  end
  r = (pi*ka).^(0.5).*exp(-ka+(sum1./pi));
  
  % First and second integrals in Eq. (VI.4) of Levine and Schwinger, used
  % to calculate the length correction.
  intvi4a = @(x, z) log(pi*abs(besselj(1,x)).*sqrt(besselj(1,x).^2 + ...
    bessely(1,x).^2))./(x.*sqrt(z^2 - x.^2));
  intvi4b = @(x, z) log(1./(2*besseli(1,x).*besselk(1,x))) ...
    ./(x.*sqrt(x.^2 + z^2));

  warning('off'); % turn off warning about integration interval if ka(1) close to zero.
  sum2 = zeros(size(ka));
  for n=1:length(ka)
    sum2(n) = integral(@(x)intvi4a(x, ka(n)), 0, ka(n));
  end
  sum3 = zeros(size(ka));
  for n=1:length(ka)
    sum3(n) = integral(@(x)intvi4b(x, ka(n)), 0, 600); % upper limit determined empirically
  end
  warning('on')
  loa = (sum2+sum3)./pi;
  R = -r.*exp(-2*1i*ka.*loa);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'dalmont' )
  
  % Unflanged pipe radiation impedance approximation (ka < 1.5) from
  % reference [1].
  [R0, delta] = unflanged(ka, ka2);
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'causse' )
  
  % Unflanged pipe radiation impedance approximation (ka < 1.5) from
  % reference [2], formula taken from reference [1].
  Zr = 1j*0.6113*ka - 1j*ka.^3 .* (0.036-0.034*log(ka) + 0.0187*ka2) + ...
    0.25*ka2 + ka.^4.*(0.0127+0.082*log(ka) - 0.023*ka2);
  R = (Zr -1) ./ (Zr + 1);
  
elseif strcmp( type, 'flanged' )

  % Infinite flange fit formula by Norris and Sheng for ka < 3.5, taken
  % from reference [1].
  %
  % A.N. Norris and I.C. Sheng. "Acoustic radiation from a circular pipe
  % with an infinite flange." Journal of Sound and Vibration, Vol. 135,
  % pp. 85-93, 1989.
  [R0, delta] = flanged(ka, ka2);
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'thickpipe' )
  
  kb = 2 * pi * f * b / c;
  [R0_uf, delta_uf] = unflanged(ka, ka2);
  [R0_f, delta_f] = flanged(ka, ka2);
  
  % [1] Combine radiation length corrections Eq. (41) Due to phase wrapping
  % issues, compute delta_c and R_c separately Use Eq. (41) in [1],
  % together with frequency-dependent length corrections for the infinite
  % flanged and unflanged pipes to get real valued thick pipe (or circular
  % flanged) length correction.
  delta_c = delta_f + (a/b) * ( delta_uf - delta_f) ...
    + 0.057 * (a/b) * ( 1 - ( a / b)^5 );

  % Real circular flange reflection coefficient as a weighted interpolation
  % between infinite flanged and unflanged coefficients.
  R0_c = R0_f.^((b - a)/b).*R0_uf.^(a/b);

  % Complex flange reflectance, but without edge effect
  R_flange = -abs(R0_c) .* exp( -2j * ka.* delta_c);

  % Eq. 42b from [1], magnitude of "edge" reflection coefficient
  R0_edge = -0.43 * ( ( b - a ) * a )/( b^2 ) * sin( kb / ( 1.85 - a/b ) ).^2;

  % Complex edge reflection coefficient
  R_edge = R0_edge .* exp( -1j*kb .* (1 + a/b * (2.3 - a/b - 0.3*(ka2))));

  % Eq. 42a from [1], combined reflectance from interpolated flange plus
  % edge effect
  R = R_flange + R_edge;
  Zr = (1 + R) ./ (1 - R);

else
  error( 'radiation: Unknown type argument.' );
end

Zr = Zr * rho * c / (pi * a.^2); % scale by Zc

end

% Dalmont unflanged approximation from [1]
function [R0, delta] = unflanged(ka, ka2)
  R0 = (1 + 0.2*ka - 0.084*ka2) ./ (1 + 0.2*ka + (0.5-0.084)*ka2); % ka<3.5
  delta = 0.6133*((1 + 0.044*ka2)./(1 + 0.19*ka2) - 0.02*sin(2*ka).^2);
end

% Norris and Sheng flanged approximation from [4]
function [R0, delta] = flanged(ka, ka2)
  R0 = (1 + 0.323*ka - 0.077*ka2) ./ (1 + 0.323*ka + (1-0.077)*ka2);
  delta = 0.8216*(1 + (0.77*ka).^2./(1 + 0.77*ka)).^(-1);
end