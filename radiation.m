function [Zr, R] = radiation( a, f, T, type, b, RH, CO2 )
% [ZR, R] = RADIATION( A, F, TYPE, T, B ) computes the radiation impedance
%      ZR (not normalized by Zc) and reflectance (or reflection
%      coefficient) R of a cylindrical pipe at frequencies specified in the
%      1D vector F (Hertz) and for pipe radius A (in meters). T is an
%      optional air temperature in degrees Celsius (default = 20 C). TYPE
%      is an optional parameter specifying a particular condition or
%      formula. The default ('unflanged-silva-noncausal') is an unflanged
%      approximation provided in [5]. Other unflanged options include: a
%      causal unflanged solution by Silva ('unflanged-silva-causal'),
%      approximate solutions by Dalmont [1] ('unflanged-dalmont') and
%      Causse [2] ('unflanged-causse-') and a computationally expensive
%      numerical integration of a solution from [3]
%      ('unflanged-ls-integration'). Infinite flanged tube options include:
%      causal and noncausal formulae provided by Silva [5]
%      ('flanged-silva-noncausal' and 'flanged-silva-causal'), and a
%      solution from [4] ('flanged-norris'). A normalized flange
%      ('normalizedflange') and a thick pipe ('thickpipe') approximation
%      are provided from [1]. The parameter B (in meters), which specifies
%      the outer radius of the pipe wall at the open end, is required with
%      the 'thickpipe' type. The other optional parameters are relative
%      humidity RH (default = 50%) and carbon dioxide percentage CO2
%      (default = 0.042%).
%
% by Gary P. Scavone, McGill University, 2013-2026.
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
% 4. A.N. Norris and I.C. Sheng, "Acoustic radiation from a circular pipe
%      with an infinite flange." Journal of Sound and Vibration, Vol. 135,
%      pp. 85-93, 1989.
%
% 5. Silva, F., Guillemain, P., Kergomard, J., Mallaroni, B., Norris, A.N.,
%      "Approximation formulae for the acoustic radiation impedance of a
%      cylindrical pipe," Journal of Sound and Vibration, Vol. 322, pp.
%      255-262, 2009.
% 
%
% Generally, Zr = Zc*(1+R)/(1-R), where R is the reflection coefficient for
% the open end of a pipe and Zc is the characteristic impedance.  R can
% also be expressed as:
%     R = -|R0|exp(-2j*k*a*delta),
% where delta is the frequency dependant length correction.  See
% reference [1], pg. 509.  A low-frequency approximation for an unflanged
% Zr is:
%     Zr = 0.25*ka^2 + 0.61j*ka

if nargin < 2 || nargin > 7
  error( 'radiation: Invalid number of arguments.');
end
if ~isvector(f)
  error( 'radiation: f should be a 1D vector of frequencies in Hertz.' );
end
if ~exist( 'T', 'var') || isempty(T)
  T = 20;
end
if ~exist( 'type', 'var') || isempty(type)
  type = 'dalmont';
end
if strcmp( type, 'thickpipe' ) && nargin < 5
  error( 'radiation: b is a required parameter for thickpipe type.' );
end
if ~exist( 'RH', 'var') || isempty(RH)
  RH = 50; % percent
end
if ~exist( 'CO2', 'var') || isempty(CO2)
  CO2 = 0.042; % percent
end

[c, rho] = thermoConstants( T, RH, CO2 );
ka = 2 * pi * f * a / c;
ka2 = ka.^2;

if strcmp( type, 'unflanged-ls-integration' )
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

elseif strcmp( type, 'unflanged-dalmont' )
  
  % Unflanged pipe radiation impedance approximation (ka < 1.5) from
  % reference [1].
  R0 = (1 + 0.2*ka - 0.084*ka2) ./ (1 + 0.2*ka + (0.5-0.084)*ka2); % ka<3.5
  delta = 0.6133*((1 + 0.044*ka2)./(1 + 0.19*ka2) - 0.02*sin(2*ka).^2);
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'unflanged-silva-noncausal')

  [absR, LoA] = unflanged(ka);
  Zr = 1j * tan( ka.*LoA + 0.5*1j*log( absR ) );
  R = -absR.*exp(-2*1i*ka.*LoA);

elseif strcmp( type, 'flanged-silva-noncausal')

  [absR, LoA] = flanged(ka);
  Zr = 1j * tan( ka.*LoA + 0.5*1j*log( absR ) );
  R = -absR.*exp(-2*1i*ka.*LoA);

elseif strcmp ( type, 'unflanged-silva-causal')
  beta = 0.5;
  eta = 0.6133;
  d2 = 0.457;
  d1 = 1.393;
  n1 = 0.167;

  % Compute d1, n1 from constraints [5] Eq. (17). This ensures eta length
  % correction is more accurate.
  d1 = eta + (d2 + beta)/(2*eta);
  n1 = (d2 + beta)/(2*eta) - eta;

  % Note Zr ~= direct computation from R in [5]. A correction is made to
  % make the equation passive.
  jka = 1j * ka;
  Zr = ( (d1 - n1)*jka + d2 * jka.^2 ) ./ (2 + (d1 + n1)*jka + d2*jka.^2);
  %R = (1 - n1*jka) ./ (1 - d1*jka + d2*jka.*jka);
  R = (Zr - 1) ./ (Zr + 1);

elseif strcmp ( type, 'flanged-silva-causal')
  beta = 1.0;
  eta = 0.8216;
  d2 = 0.649;
  %d1 = 1.825;
  %n1 = 0.182;

  % Compute d1, n1 from constraints [1] Eq. (17). This ensures eta length
  % correction is more accurate
  d1 = eta + (d2 + beta)/(2*eta);
  n1 = (d2 + beta)/(2*eta) - eta;

  % Note Zr ~= direct computation from R in [1]. A correction is made to
  % make the equation passive.
  jka = 1j * ka;
  Zr = ( -(n1 - d1)*jka + d2 * jka.^2 )./(2 + (d1 + n1)*jka + d2*jka.^2);
  %R = -(1 - n1*jka) ./ (1 - d1*jka + d2*jka.*jka);
  R = (Zr - 1) ./ (Zr + 1);

elseif strcmp( type, 'unflanged-causse' )
  
  % Unflanged pipe radiation impedance approximation (ka < 1.5) from
  % reference [2], formula taken from reference [1].
  Zr = 1j*0.6113*ka - 1j*ka.^3 .* (0.036-0.034*log(ka) + 0.0187*ka2) + ...
    0.25*ka2 + ka.^4.*(0.0127+0.082*log(ka) - 0.023*ka2);
  R = (Zr -1) ./ (Zr + 1);
  
elseif strcmp( type, 'flanged-norris' )

  % Infinite flange fit formula by Norris and Sheng for ka < 3.5, taken
  % from reference [1].
  %
  % A.N. Norris and I.C. Sheng. "Acoustic radiation from a circular pipe
  % with an infinite flange." Journal of Sound and Vibration, Vol. 135,
  % pp. 85-93, 1989.
  %[R0, delta] = flanged(ka, ka2);
  R0 = (1 + 0.323*ka - 0.077*ka2) ./ (1 + 0.323*ka + (1-0.077)*ka2);
  delta = 0.8216*(1 + (0.77*ka).^2./(1 + 0.77*ka)).^(-1);
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'normalizedflange' )
  % From [1, Eq. 44], the normalized flange is very similar to the infinite flange,
  % with a slightly shorter end correction.
  [R0, delta] = flanged(ka);
  delta = 0.975 * delta;
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'thickpipe' )
  
  kb = 2 * pi * f * b / c;
  [R0_uf, delta_uf] = unflanged(ka);
  [R0_f, delta_f] = flanged(ka);
  
  % [1] Combine radiation length corrections Eq. (41). Due to phase
  % wrapping issues, compute delta_c and R_c separately. Use Eq. (41) in
  % [1], together with frequency-dependent length corrections for the
  % infinite flanged and unflanged pipes to get real valued thick pipe (or
  % circular flanged) length correction.
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

% Silva non-causal unflanged approximation from [5]
function [absR, LoA] = unflanged(ka)
  %R0 = (1 + 0.2*ka - 0.084*ka2) ./ (1 + 0.2*ka + (0.5-0.084)*ka2); % ka<3.5
  %delta = 0.6133*((1 + 0.044*ka2)./(1 + 0.19*ka2) - 0.02*sin(2*ka).^2);
  beta = 0.5;
  eta = 0.6133;
  a1 = 0.8; a2 = 0.266; a3 = 0.0263;
  b1 = 0.0599; b2 = 0.238; b3 = -0.0153; b4 = 0.0015;
  absR = ( 1 + a1 * ka.^2 ) ./ ...
    (1 + ( beta + a1 )*ka.^2 + a2*ka.^4 + a3*ka.^6);
  LoA = eta * (1 + b1*ka.^2 ) ./ ...
    (1 + b2*ka.^2 + b3*ka.^4 + b4*ka.^6 );
end

% Silva non-causal flanged approximation from [5]
function [absR, LoA] = flanged(ka)
  %R0 = (1 + 0.323*ka - 0.077*ka2) ./ (1 + 0.323*ka + (1-0.077)*ka2);
  %delta = 0.8216*(1 + (0.77*ka).^2./(1 + 0.77*ka)).^(-1);
  beta = 1.0;
  eta = 0.8216;
  a1 = 0.730; a2 = 0.372; a3 = 0.0231;
  b1 = 0.244; b2 = 0.723; b3 = -0.0198; b4 = 0.00366;
  absR = ( 1 + a1 * ka.^2 ) ./ ...
    (1 + ( beta + a1 )*ka.^2 + a2*ka.^4 + a3*ka.^6);
  LoA = eta * (1 + b1*ka.^2 ) ./ ...
    (1 + b2*ka.^2 + b3*ka.^4 + b4*ka.^6 );
end