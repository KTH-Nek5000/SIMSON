% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [dxfou,dyfou,dzfou]=gradfield(fou,yF,kxvec,kzvec);
%
% Differentiate field in Fourier space
% Return data in Fourier space
%
% fou Two-dimensional plane in Fourier space
%

s=size(fou);
Nx=2*s(1);
Nz=s(2)+1;

if length(s)>=3; Ny=s(3); else; Ny=1;end;
if length(s)>=4; Ns=s(4); else; Ns=1;end;
if length(s)>=5; Nr=s(5); else; Nr=1;end;

NNy=Ny/3;

%
% Set up differentiation matrix including boundaries
%
[xx,DM]=chebdif(NNy,1);

dxfou=0.0*fou;
dyfou=0.0*fou;
dzfou=0.0*fou;

%
% Compute streamwise and spanwise derivatives
%
for indy=1:NNy
  for indz=1:Nz-1
    % d/dz
    dzfou(:,indz,indy)       = sqrt(-1)*kzvec(indz)*fou(:,indz,indy);
    dzfou(:,indz,indy+NNy)   = sqrt(-1)*kzvec(indz)*fou(:,indz,indy+NNy);
    dzfou(:,indz,indy+2*NNy) = sqrt(-1)*kzvec(indz)*fou(:,indz,indy+2*NNy);
  end
  for indx=1:Nx/2
    % d/dx
    dxfou(indx,:,indy)       = sqrt(-1)*kxvec(indx)*fou(indx,:,indy);
    dxfou(indx,:,indy+NNy)   = sqrt(-1)*kxvec(indx)*fou(indx,:,indy+NNy);
    dxfou(indx,:,indy+2*NNy) = sqrt(-1)*kxvec(indx)*fou(indx,:,indy+2*NNy);
  end
end

%
% Wall-normal derivatives are taken in physical space
%
foup=fou2phys(fou,0,0);
for indx=1:Nx
  for indz=1:Nz
    tmpfou(1:NNy)=foup(indx,indz,1:NNy);
    dyfoup(indx,indz,1:NNy)=(DM*tmpfou');
    
    tmpfou(1:NNy)=foup(indx,indz,NNy+1:2*NNy);
    dyfoup(indx,indz,NNy+1:2*NNy)=(DM*tmpfou');

    tmpfou(1:NNy)=foup(indx,indz,2*NNy+1:3*NNy);
    dyfoup(indx,indz,2*NNy+1:3*NNy)=(DM*tmpfou');
  end
end
dyfou=phys2fou(dyfoup);

%
% Compute wall-normal derivatives in Fourier space
% Does not work at the moment???
%
%for indx=1:Nx/2
%  for indz=1:Nz-1
    % d/dy

%    tmpfou(1:NNy)=fou(indx,indz,1:NNy);
%    dyfou(indx,indz,1:NNy)=(DM*tmpfou');
    
%    tmpfou(1:NNy)=fou(indx,indz,NNy+1:2*NNy);
%    dyfou(indx,indz,NNy+1:2*NNy)=(DM*tmpfou');

%    tmpfou(1:NNy)=fou(indx,indz,2*NNy+1:3*NNy);
%    dyfou(indx,indz,2*NNy+1:3*NNy)=(DM*tmpfou');
%  end
%end
