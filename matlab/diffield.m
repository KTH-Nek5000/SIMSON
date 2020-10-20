% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [dfou]=diffield(fou,yF,kxvec,kzvec,ex,ez);
% Differentiate field in Fourier space
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

NNx=Nx+ex
NNy=Ny/3
NNz=Nz+ez

[xx,DM]=chebdif(NNy,1);

dfou=0.0*fou;

for indy=1:NNy

  for indz=1:NNz-1
    % du/dz
    dfou(:,indz,indy+NNy)=sqrt(-1)*kzvec(indz)*fou(:,indz,indy);
    % dv/dz
    dfou(:,indz,indy)=-sqrt(-1)*kzvec(indz)*fou(:,indz,indy+NNy);
  end

  for indx=1:NNx/2
    % dw/dx
    dfou(indx,:,indy+NNy)=dfou(indx,:,indy+NNy)-sqrt(-1)*kxvec(indx)*fou(indx,:,indy+2*NNy);
    % dv/dx
    dfou(indx,:,indy+2*NNy)=sqrt(-1)*kxvec(indx)*fou(indx,:,indy+NNy);
  end

end

for indx=1:NNx/2
  for indz=1:NNz-1
    % dw/dy
    tmpfou(1:NNy)=fou(indx,indz,2*NNy+1:3*NNy);
    tmpdfou(1:NNy)=dfou(indx,indz,1:NNy);
    dfou(indx,indz,1:NNy)=tmpdfou'+(DM*tmpfou');

    % du/dy
    tmpfou2(1:NNy)=fou(indx,indz,1:NNy);
    tmpdfou2(1:NNy)=dfou(indx,indz,2*NNy+1:3*NNy);
    dfou(indx,indz,2*NNy+1:3*NNy)=tmpdfou2'-(DM*tmpfou2');
  end
end
