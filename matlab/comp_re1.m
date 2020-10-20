function [re,delta,utau,lstar,H12]=comp_re1(Ly,u,Re)
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

  delta(i,1) = interp1(u(:,i),y,0.99);
  delta(i,4) = interp1(u(:,i),y,0.95);

  
  d1=inv(DD)*(1-u(:,i)/u(N,i));
  d1 = d1-d1(1);
  
  
%  figure(100)
%  hold on
%  plot(y,d1,'b');

  
  
  d2=inv(DD)*((1-u(:,i)/u(N,i)).*u(:,i)/u(N,i));
  d2 = d2-d2(1);

%  plot(y,d2,'r');
%  plot(y,d1./d2,'g')
  
  d199 = interp1(y,d1,delta(i,1));
  d299 = interp1(y,d2,delta(i,1));


  d1=d1(N);
  d2=d2(N);
  H12(i)=d1/d2;
  delta(i,2) = d1;
  delta(i,3) = d2;
  delta(i,5) = d199;
  delta(i,6) = d299;
  re(i,1) = delta(i,1)*Re*u(N);
  re(i,2) = delta(i,2)*Re*u(N);
  re(i,3) = delta(i,3)*Re*u(N);
  re(i,4) = delta(i,4)*Re*u(N);
end

