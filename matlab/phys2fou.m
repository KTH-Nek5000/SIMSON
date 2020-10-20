% ***********************************************************************
%
% $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/matlab/fou2phys.m $
% $LastChangedDate: 2006-09-22 11:01:19 +0200 (Fri, 22 Sep 2006) $
% $LastChangedBy: mattias@MECH.KTH.SE $
% $LastChangedRevision: 147 $
%
% ***********************************************************************
function [fou]=phys2fou(phys);
% Fourier transform

s=size(phys);
Nx=s(1);
Nz=s(2);
if length(s)==3; Ny=s(3); else; Ny=1;end;

fou=zeros(Nx/2,Nz-1,Ny);

for indy=1:Ny
  % To Fourier space
  x=phys(:,:,indy);

  %x=fftshift(x);
  X=fft2(x)/(Nx*Nz);

  % Remove conjugate part and odd balls
  fou(:,:,indy)=[X(1:Nx/2,1:Nz/2) ,  X(1:Nx/2,Nz/2+2:Nz)];
end
