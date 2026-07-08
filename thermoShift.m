function shiftedZ = thermoShift( Z, f, T, deltaT )
% THERMOSHIFT: Reinterpolate impedance values based on a temperature
% difference.
%
% SHIFTEDZ = THERMOSHIFT( Z, F, T, deltaT ) takes a 1D vector of complex
% impedance values Z sampled at the frequencies specified in the 1D vector
% F, measured or calculated at a temperature T (in degrees Celsius), and
% returns a new complex impedance SHIFTEDZ (sampled at the frequencies in
% F) that (approximately) corresponds to what Z is expected to be at
% temperature T + DELTAT (where DELTAT can be negative or positive).
%
% by Gary P. Scavone, McGill University, 2026.

if nargin ~= 4
  error( 'thermoshift: Invalid number of arguments.');
end
if (~isvector(f) || ~isvector(Z)) | size(f) ~= size(Z)
  error( 'thermoshift: Z and f should be 1D vectors of the same size.' );
end

c1 = thermoConstants( T );
c2 = thermoConstants( T + deltaT );
f2 = f * c2/c1;
shiftedZ = interp1(f2, Z, f, 'spline'); 