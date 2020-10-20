function [re_theta,re_deltastar,utau,lstar,H12]=comp_re(Ly,u,Re)
%
% Computes the various Reynolds numbers based on a
% Chebyshev-distributed velocity profile
%

[N,M]=size(u);


[y,D] = chebdif(N,1);
D=D*(-2/Ly);
y=(y(N:-1:1)+1)*(Ly/2);

du=D*u;

utau = sqrt(du(1,:)/Re);
lstar = 1./utau/Re;
cf = 2*(utau./u(N,:)).^2;

DD=D;
DD(1,:)=0;
DD(1,1)=1;

for i=1:M
  d1=inv(DD)*(1-u(N:-1:1,i)/u(N,i));
  d1=d1(N);
  d2=inv(DD)*((1-u(N:-1:1,i)/u(N,i)).*u(N:-1:1,i)/u(N,i));
  d2=d2(N);
  H12(i)=d1/d2;
  re_deltastar(i) = d1*Re*u(N);
  re_theta(i) = d2*Re*u(N);
  delta99(i) = interp1(u(:,i),y,0.99);
  delta1(i) = d1;
  delta2(i) = d2;
end

delta99;
delta1;
delta2;

