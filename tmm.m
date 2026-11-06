function Zin = tmm( boreData, holeData, endType, f, lossType, T, b, RH, CO2 )
% TMM: Compute the normalized input impedance of a system using the
%      transfer matrix method.
%
% ZIN = TMM( BOREDATA, HOLEDATA, ENDTYPE, F, LOSSTYPE, T, B, RH, CO2 )
% returns the input impedance of a system defined by BOREDATA and HOLEDATA,
% normalized by the characteristic impedance at the input, at frequencies
% specified in the 1D vector F, given an optional air temperature T in
% degrees Celsius (default = 20 C), an optional relative humidity RH
% (default = 50%) and an optional carbon dioxide percentage CO2 (default =
% 0.042%).
% 
% The parameter ENDTYPE can either specify a particular bore end case
% number [0 = rigidly closed; 1 = unflanged open; 2 = flanged open; 3 =
% ideally open (Zl = 0)], or a more specific end formula ['closed',
% 'ideally-open', 'unflanged-silva-noncausal' (default for case 1),
% 'unflanged-silva-causal', 'flanged-silva-noncausal' (default for case 2),
% 'flanged-silva-causal', 'unflanged-dalmont', 'unflanged-ls-integration',
% 'flanged-norris', 'normalizedflange', 'thickpipe'] or it can be a 1D
% vector representing a pre-computed load impedance (which should have the
% same dimensions as F and should not be normalized by a characteristic
% impedance). If ENDTYPE = 'thickpipe', the parameter B, which specifies
% the outer wall radius at the end of the object, must be provided. Details
% of the various formula can be found in radiation.m.
% 
% The optional parameter LOSSTYPE specifies how losses are approximated [0
% = no losses; 1 = lowest order losses (previous tmm method, default); 2 =
% Zwikker-Kosten; 3 = full Bessel function computations]. 
%
% BOREDATA is a 2D matrix, with values in the first row corresponding to
% positions along the center axis of a specified geometry, from input to
% output ends, and values in the second row corresponding to radii at those
% positions (all values in meters).
%
% HOLEDATA is a 2D matrix specifying information about holes along a
% geometry. HOLEDATA can be empty ([]) or given by zeros(6, 0) if no holes
% exist. If holes do exist, the first row specifies positions along the
% center axis and each subsequent row specifies corresponding hole radii,
% hole heights, hole protrusion lengths, hole states (open or closed), pad
% states, pad radii, pad heights and wall thicknesses (all values, other
% than states, are in meters).
%
% Initially by Gary P. Scavone, McGill University, 2013-2024, updates
% provided by Champ Darabundit, 2026.

if nargin < 4 || nargin > 9
  error( 'tmm: Invalid number of arguments.');
end
if ~isvector(f)
  error( 'tmm: f should be a 1D vector of frequencies in Hertz.' );
end
if ~isscalar(endType) && ~ischar(endType)
  if length(endType) ~= length(f)
    error( 'tmm: endType impedance vector must be the same length as F.' );
  elseif size(endType, 1) ~= size(f, 1)
    % Transpose if not same dimension
    endType = endType';
  end
elseif isscalar(endType)
  if endType < 0 || endType > 3
    error('tmm: scalar endType must be between 0 - 3.');
  end
end
if ~exist( 'T', 'var') || isempty(T)
  T = 20;
end
if ~exist( 'b', 'var') || isempty(b)
  b = 0;
end
if ~exist( 'RH', 'var') || isempty(RH)
  RH = 50; % percent
end
if ~exist( 'CO2', 'var') || isempty(CO2)
  CO2 = 0.042; % percent
end
if ~exist( 'lossType', 'var')
  lossType = 1;
end
if isempty( holeData )
  holeData = zeros(6, 0);
end

% Bore dimensions
idx = find(diff(boreData(1,:)) == 0);
boreData(1, idx+1) = boreData(1, idx+1) + eps; % avoid double values
x = sort( [boreData(1,:) holeData(1,:)] ); % segment positions along x-axis
L = diff( x );                             % lengths of segments

% Interpolate bore radii at x values
ra = interp1(boreData(1,:), boreData(2,:), x, 'linear');

isHole = zeros(size(x));                   % is x value at a tonehole?
for n = 1:length(x)
  isHole(n) = 1 - isempty(find(x(n)==holeData(1,:), 1));
end
nHole = sum(isHole);

if nHole
  % Tonehole dimensions and states
  rb = holeData(2,:).';            % tonehole radii
  t = holeData(3,:).';             % tonehole heights
  chimney = holeData(4,:);         % tonehole chimney height
  states = holeData(5,:);          % tonehole states
  [n, m] = size(holeData);
  padr = zeros(1, m);
  padt = zeros(1, m);
  holew = zeros(1, m);
  if n > 6, padr = holeData(7,:); end % tonehole pad radii
  if n > 7, padt = holeData(8,:); end % tonehole pad heights
  if n > 8, holew = holeData(9,:); end % tonehole wall thickness
end

if f(1) == 0 % avoid zero frequency calculations
  f(1) = eps;
end

% Work our way back from the load impedance at the end.
if isscalar(endType) || ischar(endType)
  switch endType
    case {0, 'closed'}
      Zl = Inf;
    case {1, 'unflanged-silva-noncausal'}
      Zl = radiation( ra(end), f, T, 'unflanged-silva-noncausal', [], RH, CO2 );
    case {2, 'flanged-silva-noncausal'}
      Zl = radiation( ra(end), f, T, 'flanged-silva-noncausal', [], RH, CO2 );
    case {3, 'ideally-open'}
      Zl = 0;
    case {'thickpipe'}
      Zl = radiation( ra(end), f, T, 'thickpipe', b, RH, CO2 );
    otherwise
      Zl = radiation( ra(end), f, T, endType, [], RH, CO2 );
  end
else
  Zl = endType;
end

for n = length(L):-1:1
  if L(n) > eps
    [Gamma, Zc] = sectionLosses( ra(n), ra(n+1), L(n), f, T, lossType, RH, CO2 );
    [A, B, C, D] = tmmCylCone( ra(n), ra(n+1), L(n), Gamma, Zc );
    if Zl == 0
      Zl = B ./ D;
    else
      Zl = (A + B./Zl) ./ (C + D./Zl);
    end
  end
  
  if isHole(n)
    Gamma = sectionLosses( rb(nHole), rb(nHole), 0, f, T, lossType, RH, CO2 );
    [A, B, C, D] = tmmTonehole( rb(nHole)/ra(n), rb(nHole), t(nHole), ...
      states(nHole), Gamma, '', T, chimney(nHole), padr(nHole), ...
      padt(nHole), holew(nHole) );
    nHole = nHole - 1;
    Zl = (A + B./Zl) ./ (C + D./Zl);
  end
end

% Calculate real wave impedance at entrance for scaling
[c, rho] = thermoConstants( T, RH, CO2 );
Zc = rho * c / ( pi * ra(1) * ra(1) );
Zin = Zl ./ Zc;
