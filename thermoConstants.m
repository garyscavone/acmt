function [c, rho, gamma, lv, Pr] = thermoConstants( T, RH, CO2 )
% THERMOCONSTANTS: Determine various thermodynamic constants for a
% specified temperature, relative humidity and CO2 level.
%
% [C, RHO, GAMMA, LV, PR] = THERMOCONSTANTS( T, RH, CO2 ) returns the speed
% of sound C, the density of air RHO, the ratio of specific heats GAMMA,
% the viscous characteristic length LV, and the Prandtl number PR, based on
% a temperature T in degrees Celsius (default T = 20 degrees Celsius),
% relative humidity RH in percent (default RH = 50%) and a percentage
% carbon dioxide level CO2 (default CO2 = 0.042%).
%
% Reference:
%
% 1. D. Keefe, "Acoustical wave propagation in cylindrical ducts:
%      Transmission line parameter approximations for isothermal and
%      nonisothermal boundary conditions," Journal of the Acoustical
%      Society of America, Vol. 75, No. 1, pp. 58-62, 1984.
%
% 2. A. Ernoult, "Effect of air humidity and carbon dioxide in the sound
%    propagation for the modeling of wind musical instruments," Research
%    Report RR-9500, Inria, 2023. URL: https://hal.inria.fr/hal-04008847.
%
% by Gary P. Scavone and Champ Darabundit, McGill University, 2013-2026.

if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'RH', 'var')
  RH = 50; % percent
end
if ~exist( 'CO2', 'var')
  CO2 = 0.042; % percent
end

% Convert T in Celsius to Kelvin, RH and CO2 in percentage to decimal
T = T + 273.15;
RH = RH / 100;
CO2 = CO2 / 100;

% Reference values
Tref = 293.15;
xcref = 4.2e-4;
xvref = 1.157e-2;

% Compute molar fraction of water vapor [1] Eq. (53)
xv = RH * 10^( 5.21899 - 5.8294*( Tref / T ) - 1.0252 * ( Tref / T )^2 );

% Compute difference values for xv, xc and T from reference values
xvDiff = xv - xvref;
xcDiff = CO2 - xcref;
TDiff = T/Tref - 1;

% Compute specific heat capacity in J/(K*kg) [1] Eq. (54)
Cp = 1012.25 * ( 1 + 0.5438 * xvDiff + 0.638 * xvDiff^2 ...
     - 0.1594 * xcDiff + 0.075 * xcDiff^2 ...
     + 0.00952 * TDiff + 0.0406 * TDiff^2 + 0.3976 * xcDiff * TDiff );

% Ratio of specific heats [1] Eq. (55)
gamma = 1.40108 * ( 1 - 0.060 * xvDiff - 0.104 * xcDiff ...
        - 0.0087 * TDiff - 0.154 * xcDiff * TDiff );

% Speed of sound in m/s [1] Eq. (56)
c = 343.986 * sqrt( T / Tref ) * sqrt( 1 + 0.314 * xvDiff - 0.520 * xcDiff ...
    + 0.25 * xcDiff^2 - 0.16 * xcDiff * TDiff );

% Air density in kg / m^3 [1] Eq. (57)
rho = 1.19930 * ( Tref / T ) * ( 1 - 0.3767 * xvDiff + 0.4162 * xcDiff ...
      - 0.0029 * TDiff );

% Compute viscosity and thermal conductivity in kg/(m*s) [1] Eq. (58)
% The results using Eq. (58) seem to deviate significantly from other
% sources, including the values given in Table 6 of [1], so we are
% currently using the Table 6 for dry air.
%mu = 1.8206e-5 * (1 + 0.77013 * TDiff );
mu = (-9.8601 + 0.90801*T - 1.1764e-3*(T.^2) + 1.2350e-6*(T.^3) - ...
     5.7971e-10*(T.^4))*10^(-7);

% Compute thermal conductivity in kg/(m*s) [1] Eq. (59)
kappa = 2.5562e-2 * (1 + 0.8490 * TDiff ); % W / (m*K)

Pr = mu * Cp / kappa;
lv = mu / (rho * c);

% Values from Keefe (1984)
%deltaT = T - 26.85;
%c = 347.23 * ( 1 + 0.00166 * deltaT );        % speed of sound in air (m/s)
%rho = 1.1769 * ( 1 - 0.00335 * deltaT );      % density of air (kg/m^3)
%mu = 1.846*10^(-5) * ( 1 + 0.0025 * deltaT ); % shear viscosity coefficient (kg/m s)
%gamma = 1.4017 * ( 1 - 0.00002 * deltaT );    % ratio of specific heats
%Pr = (0.8410 * ( 1 - 0.00002 * deltaT ))^2;   % Prandtl number
