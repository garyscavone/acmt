function shiftedZ = zShift( Z, f, torf1, torf2, type )
% ZSHIFT: Reinterpolate impedance values based on either a temperature
% or a frequency difference.
%
% SHIFTEDZ = ZSHIFT( Z, F, TORF1, TORF2, TYPE ) takes a 1D vector of
% complex impedance values Z sampled at the frequencies specified in the 1D
% vector F, and returns a new complex impedance SHIFTEDZ (sampled at the
% frequencies in F). If TYPE is equal to 'temperature' (default), the shift
% is calculated as corresponding to what Z is expected to be at temperature
% TORF2 if it was originally measured or calculated at TORF1 (in Celsius).
% If TYPE is equal to 'frequency', the shift is calculated to shift a peak
%
% by Gary P. Scavone, McGill University, 2026.

if nargin <= 3
  error( 'zshift: Invalid number of arguments.');
end
if nargin < 5, type = 'temperature'; end

if (~isvector(f) || ~isvector(Z)) | size(f) ~= size(Z)
  error( 'thermoshift: Z and f should be 1D vectors of the same size.' );
end

if strcmp( type, 'temperature' )
  c1 = thermoConstants( torf1 );
  c2 = thermoConstants( torf2 );
  f2 = f * c2/c1;
elseif strcmp( type, 'frequency' )
  f2 = f * torf2 / torf1;
else
  error( 'zshift" Invalid type argument.');
end

shiftedZ = interp1(f2, Z, f, 'cubic'); 