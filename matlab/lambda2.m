% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [l2,dxfoup,dyfoup,dzfoup]=lambda2(fou,yF,kxvec,kzvec,ex,ez);
%
% Function that computes the second eigenvalue of the velocity
% gradient tensor. Lambda2 is returned in physical space
%

%
% Compute velocity gradient tensor
%
[dxfou,dyfou,dzfou]=gradfield(fou,yF,kxvec,kzvec);

%
% Transform velocity gradient tensor to physical space
%
[dxfoup,NNx,NNy,NNz]=fou2phys(dxfou,ex,ez);
dyfoup=fou2phys(dyfou,ex,ez);
dzfoup=fou2phys(dzfou,ex,ez);

%
% Construct symmetric and skew-symmetric local tensors
%
for k=1:NNz
  for j=1:NNy
    for i=1:NNx

      sr(1,1) = dxfoup(i,k,j);
      sr(2,1) = .5*(dxfoup(i,k,NNy+j)+dyfoup(i,k,j));
      sr(3,1) = .5*(dxfoup(i,k,2*NNy+j)+dzfoup(i,k,j));

      sr(1,2) = .5*(dyfoup(i,k,j)+dxfoup(i,k,NNy+j));
      sr(2,2) = dyfoup(i,k,NNy+j);
      sr(3,2) = .5*(dyfoup(i,k,2*NNy+j)+dzfoup(i,k,NNy+j));

      sr(1,3) = .5*(dzfoup(i,k,j)+dxfoup(i,k,2*NNy+j));
      sr(2,3) = .5*(dzfoup(i,k,NNy+j)+dyfoup(i,k,2*NNy+j));
      sr(3,3) = dzfoup(i,k,2*NNy+j);

      or(1,1) = 0.;
      or(2,1) = .5*(dxfoup(i,k,NNy+j)-dyfoup(i,k,j));
      or(3,1) = .5*(dxfoup(i,k,2*NNy+j)-dzfoup(i,k,j));

      or(1,2) = .5*(dyfoup(i,k,j)-dxfoup(i,k,NNy+j));
      or(2,2) = 0.;
      or(3,2) = .5*(dyfoup(i,k,2*NNy+j)-dzfoup(i,k,NNy+j));

      or(1,3) = .5*(dzfoup(i,k,j)-dxfoup(i,k,2*NNy+j));
      or(2,3) = .5*(dzfoup(i,k,NNy+j)-dyfoup(i,k,2*NNy+j));
      or(3,3) = 0.;
               
      for ii=1:3
	for kk=1:3
	  a(ii,kk) = 0.;
	  for jj=1:3
	    a(ii,kk) = a(ii,kk) + sr(ii,jj)*sr(jj,kk) + or(ii,jj)*or(jj,kk);
	  end
	end
      end
      %
      % Solve eigenvalue problem for each local matrix
      %
      w=eig(a);

      sort(w);

      l2(i,k,j) = w(2);
      
    end
  end
end
