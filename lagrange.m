function b = lagrange( N, D )
% LAGRANGE: Returns numerator coefficients of a Lagrange interpolation
%           filter of order N for the fractional delay D (in samples).
%
% by Julius Smith, CCRMA.
% https://ccrma.stanford.edu/~jos/pasp/Matlab_Code_Lagrange_Interpolation.html

n = 0:N;
b = ones(1,N+1);
for k = 0:N
  index = find(n ~= k);
  b(index) = b(index) *  (D-k)./ (n(index)-k);
end