% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [Retau]=retauf(foup,yF,kxvec,kzvec,Re);
%
%Compute Retau
%
sp=size(foup);
NNx=sp(1);
NNz=sp(2);
NNy=sp(3);
%
% Set up differentiation matrix including boundaries
%
%[xx,DM]=chebdif(NNy,1);

   yy=yF(1:end);
   my= fds(1, NNy, 4, 0,yy);
%    size(my)
%    mx1=mx1*2.0/(xF(end)-xF(1));
% dxfoup=zeros(NNx+1,NNz,3*NNy);

%
% Wall-normal derivatives are taken in physical space
%
%foup=fou2phys(fou(:,:,1:NNy),0,0);
for indx=1:NNx
  for indz=1:NNz
    tmpfou=reshape(foup(indx,indz,1:NNy),NNy,1);
%    dyfoup(indx,indz,1:NNy)=(DM*tmpfou');
    
%     dudyw1(indx,indz)=(DM(1,:)*tmpfou');
%     dudyw2(indx,indz)=(DM(NNy,:)*tmpfou');
    dudyw1(indx,indz)=(my(1,:)*tmpfou);
    dudyw2(indx,indz)=(my(NNy,:)*tmpfou);

  end
end
dudyw=sum(sum(abs(dudyw1(:,:)-dudyw2(:,:))))/(NNx*NNz*2);
Retau=sqrt(Re*dudyw);


