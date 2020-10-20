function [utau,lstar]=comp_re_channel(u,Re)
%
% Computes the various Reynolds numbers based on a
% Chebyshev-distributed velocity profile
%

[N,M]=size(u);

[y,D] = chebdif(N,1);

du=D*u;

utau = sqrt(abs(du(1,:))/Re);
lstar = 1./utau/Re;

