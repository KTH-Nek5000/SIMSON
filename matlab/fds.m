% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [M1] = fds(CDmode, n, grid, alpha,xx)

% close all

% This file serves as a test file for different compact finite difference
% schemes. The Helmholtz equation (Laplacian + lambda^2) u = 0 is solved for several
% schemes stemming from "Compact Finite Difference Schemes with
% Spectral-like Resolution" by Lele (1991).
% 
% The CDmode-option designates which discretisation scheme to use:
% CDmode = 1 : tridiagonal fourth order scheme
% CDmode = 2 : tridiagonal sixth order scheme
% CDmode = 3 : pentadiagonal tenth order scheme
%
% n is the numbers of points
% alpha is a parameter defining the gridspacing
% grid is a parameter deciding to use Gauss-Lobato (grid ~=1) or a
% stretched grid

% Created by Peter Lenaers, Modified by Liang Wei (2012).


%global h x k

% Grid definition:
a = -1;
b = 1;
h = (b-a)/(n-1);    % average grid spacing

x = zeros(n,1);

if grid == 1
    temp = 1. / (n-1);
    for i = 1:n
        x(i) = (i-1)*temp;
        x(i) = tanh(alpha * (x(i) - 0.5)) / tanh(0.5 * alpha);
    end
elseif grid == 2 % Gauss-Lobato points
    temp = pi / (n-1);
    for i = 1:n
        x(i) = 1+ cos((i-1)*temp);
    end
elseif grid == 3
    for i = 1:n
        x(i) = a + (i-1)*h;
    end
elseif grid ==4 % specified grid points
    for i=1:n
        x(i) = xx(i);
    end
    
end


% wave number
k = 0:pi/(n-1):pi;

% In this function the coefficients of the compact difference matrix for the
% first derivative are calculated. If the spacing is
% equidistant, these coefficients can be found in the article by Lele.
%
% The code differentiates between schemes with different accuracy. The
% schemes are tri- or pentadiagonal on the LHS, and penta- or heptadiagonal
% in the RHS.
dX = diff(x);
P1 = sparse(n,n);
Q1 = sparse(n,n);

if CDmode == 1 % using the implicit tridiagonal, explicit tridiagonal in the function, fourth order method
    % We calculate the coefficients for
    % a_(i-1)*f'_(i-1) + f'_i + a_(i+1)*f'_(i+1) = b_(i-1)*f_(i-1) + ... + b_(i+1)*f_(i+1)
    
    % Firstly, the coefficients of the interior nodes are calculated
    for i = 2:n-1
       Acoef = sparse(5,5);
       bcoef = sparse(5,1);
       bcoef(2) = 1;
       h1 = -dX(i-1);
       h2 =  dX(i);
       
       Acoef(1,3) = 1;
       Acoef(1,4) = 1;
       Acoef(1,5) = 1;
       
       for j = 1:4
          Acoef(j+1,1) = -h1^(j-1) * factorial(j);
          Acoef(j+1,2) = -h2^(j-1) * factorial(j);
          Acoef(j+1,3) =  h1^ j    * factorial(j-1);
          Acoef(j+1,5) =  h2^ j    * factorial(j-1);
       end
       coef = Acoef \ bcoef;
       P1(i,i-1) = coef(1);
       P1(i,i)   = 1;
       P1(i,i+1) = coef(2);
       Q1(i,i-1) = coef(3);
       Q1(i,i)   = coef(4);
       Q1(i,i+1) = coef(5);
    end
    
    % The points near the boundary need to be treated separately since the full
    % scheme uses gridpoints outside the domain. To assure stability of the
    % scheme, the order near the boundary needs to be lowered.
    %
    % i = 1
    % We calculate the coefficients for f'_1 + a_2*f'_2 = b_1*f_1 + b_2*f_2 + b_3*f_3
    Acoef = sparse(4,4);
    bcoef = sparse(4,1);
    bcoef(2) = 1;
    
    h1 = dX(1);
    h2 = h1 + dX(2);
    
    Acoef(1,2) = 1;
    Acoef(1,3) = 1;
    Acoef(1,4) = 1;
    
    for i = 1:3
       Acoef(i+1,1) = -h1^(i-1) * factorial(i);
       Acoef(i+1,3) =  h1^i     * factorial(i-1);
       Acoef(i+1,4) =  h2^i     * factorial(i-1);
    end
    coef = Acoef \ bcoef;
    P1(1,1) = 1;
    P1(1,2) = coef(1);
    Q1(1,1) = coef(2);
    Q1(1,2) = coef(3);
    Q1(1,3) = coef(4);
    
    % i = n
    % We calculate the coefficients for a_(n-1)*f'_(n-1) + f'_n = b_(n-2)*f_(n-2) + ... + b_n*f_n
    Acoef = sparse(4,4);
    bcoef = sparse(4,1);
    bcoef(2) = 1;
    
    h2 = - dX(n-1);
    h1 = h2 - dX(n-2);
    
    Acoef(1,2) = 1;
    Acoef(1,3) = 1;
    Acoef(1,4) = 1;
    
    for i = 1:3
       Acoef(i+1,1) = -h2^(i-1) * factorial(i);
       Acoef(i+1,2) =  h1^i     * factorial(i-1);
       Acoef(i+1,3) =  h2^i     * factorial(i-1);
    end
    coef = Acoef \ bcoef;
    P1(n,n)   = 1;
    P1(n,n-1) = coef(1);
    Q1(n,n-2) = coef(2);
    Q1(n,n-1) = coef(3);
    Q1(n,n)   = coef(4);
    
    I = sqrt(-1);
    num = Q1(2,2) + Q1(2,1)*exp(I*k*(x(1) - x(2))/h) + Q1(2,3)*exp(I*k*(x(3) - x(2))/h);
    denom = 1 + P1(2,1)*exp(I*k*(x(1) - x(2))/h) + P1(2,3)*exp(I*k*(x(3) - x(2))/h);
    k1 = -h * I * num ./ denom;
    
elseif CDmode == 2 % using the implicit tridiagonal, explicit pentadiagonal, sixth order method
   % We calculate the coefficients for
   % a_(i-1)*f'_(i-1) + f'_i + a_(i+1)*f'_(i+1) = b_(i-2)*f_(i-2) + ... + b_(i+2)*f_(i+2)
   
   % Firstly, the coefficients of the interior nodes are calculated
   for i = 3:n-2
      Acoef = sparse(7,7);
      bcoef = sparse(7,1);
      bcoef(2) = 1;
      
      h1 = -dX(i-1) - dX(i-2);
      h2 = -dX(i-1);
      h3 =  dX(i);
      h4 =  dX(i) + dX(i+1);
      
      Acoef(1,3) = 1;
      Acoef(1,4) = 1;
      Acoef(1,5) = 1;
      Acoef(1,6) = 1;
      Acoef(1,7) = 1;
      
      for j = 1:6
         Acoef(j+1,1) = -h2^(j-1) * factorial(j);
         Acoef(j+1,2) = -h3^(j-1) * factorial(j);
         Acoef(j+1,3) =  h1^ j    * factorial(j-1);
         Acoef(j+1,4) =  h2^ j    * factorial(j-1);
         Acoef(j+1,6) =  h3^ j    * factorial(j-1);
         Acoef(j+1,7) =  h4^ j    * factorial(j-1);
      end
      coef = Acoef \ bcoef;
      P1(i,i-1) = coef(1);
      P1(i,i)   = 1;
      P1(i,i+1) = coef(2);
      Q1(i,i-2) = coef(3);
      Q1(i,i-1) = coef(4);
      Q1(i,i)   = coef(5);
      Q1(i,i+1) = coef(6);
      Q1(i,i+2) = coef(7);
   end
   
   % The points near the boundary need to be treated separately since the full
   % scheme uses gridpoints outside the domain
   %
   % i = 1
   % We calculate the coefficients for f'_1 + a_2*f'_2 = b_1*f_1 + ... + b_3*f_3
   Acoef = sparse(4,4);
   bcoef = sparse(4,1);
   bcoef(2) = 1;
   
   h1 = dX(1);
   h2 = h1 + dX(2);
   
   Acoef(1,2) = 1;
   Acoef(1,3) = 1;
   Acoef(1,4) = 1;
   
   for j = 1:3
      Acoef(j+1,1) = -h1^(j-1) * factorial(j);
      Acoef(j+1,3) =  h1^j     * factorial(j-1);
      Acoef(j+1,4) =  h2^j     * factorial(j-1);
   end
   coef = Acoef \ bcoef;
   P1(1,1) = 1;
   P1(1,2) = coef(1);
   Q1(1,1) = coef(2);
   Q1(1,2) = coef(3);
   Q1(1,3) = coef(4);
   
   % i = n
   % We calculate the coefficients for a_(n-1)*f'_(n-1) + f'_n = b_(n-2)*f_(n-2) + ... + b_n*f_n
   Acoef = sparse(4,4);
   bcoef = sparse(4,1);
   bcoef(2) = 1;
   
   h1 = - dX(n-2) - dX(n-1);
   h2 = - dX(n-1);
   
   Acoef(1,2) = 1;
   Acoef(1,3) = 1;
   Acoef(1,4) = 1;
   
   for j = 1:3
      Acoef(j+1,1) = -h2^(j-1) * factorial(j);
      Acoef(j+1,2) =  h1^j     * factorial(j-1);
      Acoef(j+1,3) =  h2^j     * factorial(j-1);
   end
   coef = Acoef \ bcoef;
   P1(n,n-1) = coef(1);
   P1(n,n)   = 1;   
   Q1(n,n-2) = coef(2);
   Q1(n,n-1) = coef(3);
   Q1(n,n)   = coef(4);
   
   % i = 2
   % We calculate the coefficients for a_1*f'_1 + f'_2 + a_3*f'_3 = b_1*f_1 + ... + b_4*f_4
   Acoef = sparse(6,6);
   bcoef = sparse(6,1);
   bcoef(2) = 1;
   
   h1 = -dX(1);
   h2 =  dX(2);
   h3 = h2 + dX(3);
   
   Acoef(1,3) = 1;
   Acoef(1,4) = 1;
   Acoef(1,5) = 1;
   Acoef(1,6) = 1;
   
   for j = 1:5
      Acoef(j+1,1) = -h1^(j-1) * factorial(j);
      Acoef(j+1,2) = -h2^(j-1) * factorial(j);
      Acoef(j+1,3) =  h1^j     * factorial(j-1);
      Acoef(j+1,5) =  h2^j     * factorial(j-1);
      Acoef(j+1,6) =  h3^j     * factorial(j-1);
   end
   coef = Acoef \ bcoef;
   P1(2,1) = coef(1);
   P1(2,2) = 1;
   P1(2,3) = coef(2);
   Q1(2,1) = coef(3);
   Q1(2,2) = coef(4);
   Q1(2,3) = coef(5);
   Q1(2,4) = coef(6);
   
   % i = n-1
   % We calculate the coefficients for a_(n-2)*f'_(n-2) + f'_(n-1) + a_n*f'_n = 
   % b_(n-3)*f_(n-3) + ... + b_n*f_n
   Acoef = sparse(6,6);
   bcoef = sparse(6,1);
   bcoef(2) = 1;
   
   h1 = - dX(n-2) - dX(n-3);
   h2 = - dX(n-2);
   h3 = dX(n-1);
   
   Acoef(1,3) = 1;
   Acoef(1,4) = 1;
   Acoef(1,5) = 1;
   Acoef(1,6) = 1;
   
   for j = 1:5
      Acoef(j+1,1) = -h2^(j-1) * factorial(j);
      Acoef(j+1,2) = -h3^(j-1) * factorial(j);
      Acoef(j+1,3) =  h1^j     * factorial(j-1);
      Acoef(j+1,4) =  h2^j     * factorial(j-1);
      Acoef(j+1,6) =  h3^j     * factorial(j-1);
   end
   coef = Acoef \ bcoef;
   P1(n-1,n-2) = coef(1);
   P1(n-1,n-1) = 1;
   P1(n-1,n)   = coef(2);
   Q1(n-1,n-3) = coef(3);
   Q1(n-1,n-2) = coef(4);
   Q1(n-1,n-1) = coef(5);
   Q1(n-1,n)   = coef(6);
   
end

M1=P1\Q1;

