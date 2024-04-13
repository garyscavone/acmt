function [B, A] = thiran( N, D )
% THIRAN: Returns coefficients of a Thiran allpass interpolation filter of
%         order N for the fractional delay D (samples).
%
% by Julius Smith, CCRMA.
% https://ccrma.stanford.edu/~jos/pasp/Thiran_Allpass_Interpolation_Matlab.html

A = zeros(1,N+1);
for k=0:N
  Ak = 1;
  for n=0:N
    Ak = Ak * (D-N+n)/(D-N+k+n);
  end
  A(k+1) = (-1)^k * nchoosek(N,k) * Ak;
end

B = A(N+1:-1:1);