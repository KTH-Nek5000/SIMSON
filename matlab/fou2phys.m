% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [phys,NNx,NNy,NNz]=fou2phys(fou,ex,ez,varargin);
% Fourier transform
% fou Two-dimensional plane in Fourier space
% ex  Number of extra zeros in x-direction (0 normally) to
%     improve interpolation to physical space
% ez  Number of extra zeros in z-direction (0 normally)
%     improve interpolation to physical space

s=size(fou);
Nx=2*s(1);
Nz=s(2)+1;

if length(s)>=3; Ny=s(3); else; Ny=1;end;
if length(s)>=4; Ns=s(4); else; Ns=1;end;
if length(s)>=5; Nr=s(5); else; Nr=1;end;

if size(varargin) > 0
  ncomp = varargin{1}
else
  ncomp=3
end
NNx=Nx+ex;
NNy=Ny/ncomp;
NNz=Nz+ez;

phys=zeros(Nx+ex,Nz+ez,Ny,Ns,Nr);

for inds=1:Ns
  for indy=1:Ny
    for indr=1:Nr
      % The z odd ball
      tp=fou(:,:,indy,inds,indr);
      tp=[tp(:,1:Nz/2) , zeros(Nx/2,1+ez) , tp(:,Nz/2+1:Nz-1)];

      % The conjugate and x odd ball
      X=[ tp ;
      zeros(1+ex,Nz+ez) ; 
      conj(flipud(tp(2:Nx/2,1))) , conj(flipud(fliplr(tp(2:Nx/2,2:Nz+ez))))];

      % To physical space
      x=ifft2(X)*(Nx+ex)*(Nz+ez);
%      x=fftshift(x); 

      phys(:,:,indy,inds,indr)=real(x);

    end
  end
end

