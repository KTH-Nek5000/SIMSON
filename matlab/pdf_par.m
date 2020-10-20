% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [pc,p2,ps,pdfylcx,pdfyl2x,pdfylsx,NNy, Q, R, QWmean,pcy,p2y,psy]=pdf_par(u,v,w, xF,yF,zF,kxvec,kzvec,ex,ez,lb,ub,step,nproc,xfour,xfd,zfd,dim);
%
% Function that computes the p.d.f. of the cosine of the angle between
% the vorticity vector and three eigenvectors of the strain rate tensor,
% Q R of velocity gradient tensor and mean Q of the rotation tensor. 
% Compute the p.d.f. of the cosine of the angle between
% three eigenvectors of the strain rate tensor and the "dim" direction
% 
%

% %
% % Compute velocity gradient tensor
% %
% [dxfou,dyfou,dzfou]=gradfield(fou,yF,kxvec,kzvec);
% 
% %
% % Transform velocity gradient tensor to physical space
% %
% [dxfoup,NNx,NNy,NNz]=fou2phys(dxfou,ex,ez);
% dyfoup=fou2phys(dyfou,ex,ez);
% dzfoup=fou2phys(dzfou,ex,ez);
% 

%    [phys,NNx,NNy,NNz]=fou2phys(fou,0,0);
%sp=size(phys);
%NNx=sp(1);
%NNz=sp(2);
%NNy=sp(3)/3;

sp=size(u);
NNx=sp(1);
NNz=sp(2);
NNy=sp(3);

nytest=floor(NNy);
%
% Compute velocity gradient tensor
%

%dxfoup=zeros(NNx,NNz,3*NNy);
dyfoup=zeros(NNx,NNz,3*NNy);
dzfoup=zeros(NNx,NNz,3*NNy);

yy=yF(1:nytest);
zz=zF(1:NNz);


%xfd=1;
%foudif=0;
%matgrad=0;
dxfoup=zeros(NNx,NNz,3*NNy);
if     xfd==1
%
% Use finite difference scheme on x direction
%
	xx=xF(1:NNx);
	[mx1]= fds(1, NNx, 4, 0,xx);
%    mx1=mx1*2.0/(xF(end)-xF(1));
%	 dxfoup=zeros(NNx,NNz,3*NNy);
 
 
%	for j=1:nytest
%	    for k=1:NNz
%        u(NNx+1,k,j)=u(1,k,j);
%        v(NNx+1,k,j)=v(1,k,j);
%        w(NNx+1,k,j)=w(1,k,j);
%	         dxfoup(:,k,j)=mx1*u(:,k,j);
%	         dxfoup(:,k,j+NNy)=mx1*v(:,k,j);
%	         dxfoup(:,k,j+NNy*2)=mx1*w(:,k,j);
%
%	    end
%	end


else
%
% Use fourier difference scheme on x direction
%
	xx=xF(1:end-1);

	[xx mx1]=fourdif(NNx,1);
	clearvars xx
	mx1=mx1*(2.0*pi/abs(xF(end)-xF(1)));
end
%	 dxfoup=zeros(NNx,NNz,3*NNy);
    

for j=1:nytest
    for k=1:NNz
         dxfoup(:,k,j)=mx1*u(:,k,j);
         dxfoup(:,k,j+NNy)=mx1*v(:,k,j);
         dxfoup(:,k,j+NNy*2)=mx1*w(:,k,j);
    end
end

%end



if zfd==1
    	zz=zF(1:NNz);
	[mz1]= fds(1, NNz, 4, 0,zz);
else

	[zz mz1]=fourdif(NNz,1);
	clearvars zz
	mz1=mz1*(2.0*pi/abs(zF(end)-zF(1)));
end
for j=1:nytest
	for i=1:NNx
         dzfoup(i,:,j)=mz1*u(i,:,j)';
         dzfoup(i,:,j+NNy)=mz1*v(i,:,j)';
         dzfoup(i,:,j+NNy*2)=mz1*w(i,:,j)';
	end
end


[my1]= fds(1, nytest, 4, 0,yy);
for k=1:NNz
    for i=1:NNx
        dyfoup(i,k,1:nytest)=my1*reshape(u(i,k,1:nytest),nytest,1);
        dyfoup(i,k,NNy+1:NNy+nytest)=my1*reshape(v(i,k,1:nytest),nytest,1);

        dyfoup(i,k,NNy*2+1:NNy*2+nytest)=my1*reshape(w(i,k,1:nytest),nytest,1);
    end
end


QWtemp=zeros(NNy,1);
Q=zeros(NNx,NNz,NNy);
R=zeros(NNx,NNz,NNy);


nbins=floor((ub-lb)/step+1);
pylcx=zeros(nbins,1);
 pylsx=zeros(nbins,1);
 pyl2x=zeros(nbins,1);

 pdfylcx=zeros(nbins,NNy);
 pdfyl2x=zeros(nbins,NNy);
 pdfylsx=zeros(nbins,NNy);

 pcy=zeros(nbins,NNy);
 p2y=zeros(nbins,NNy);
 psy=zeros(nbins,NNy);
 


 try
     matlabpool(nproc)
 catch
     matlabpool close
     matlabpool(nproc)
 end
for j=1:NNy

   for k=1:NNz
    for i=1:NNx
        


         vgt =[dxfoup(i,k,j) dyfoup(i,k,j) dzfoup(i,k,j);
            dxfoup(i,k,NNy+j) dyfoup(i,k,NNy+j) dzfoup(i,k,NNy+j);
            dxfoup(i,k,2*NNy+j) dyfoup(i,k,2*NNy+j) dzfoup(i,k,2*NNy+j)];
      or =[0. .5*(vgt(1,2)-vgt(2,1)) .5*(vgt(1,3)-vgt(3,1));
          .5*(vgt(2,1)-vgt(1,2)) 0. .5*(vgt(2,3)-vgt(3,2));
          .5*(vgt(3,1)-vgt(1,3)) .5*(vgt(3,2)-vgt(2,3)) 0.];
      
        for ll=1:3
            for mm=1:3
                Q(i,k,j) = Q(i,k,j) + (-0.5)*vgt(ll,mm)*vgt(mm,ll);
                QWtemp(j) = QWtemp(j) + (-0.5)*or(ll,mm)*or(mm,ll);
                for nn=1:3
                    R(i,k,j) = R(i,k,j) + (-1./3.)*vgt(ll,mm)*vgt(mm,nn)*vgt(nn,ll);
                end
            end
        end
        
    
         
      sr =[dxfoup(i,k,j) .5*(dyfoup(i,k,j)+dxfoup(i,k,NNy+j)) .5*(dzfoup(i,k,j)+dxfoup(i,k,2*NNy+j));
          .5*(dxfoup(i,k,NNy+j)+dyfoup(i,k,j)) dyfoup(i,k,NNy+j) .5*(dzfoup(i,k,NNy+j)+dyfoup(i,k,2*NNy+j));
          .5*(dxfoup(i,k,2*NNy+j)+dzfoup(i,k,j)) .5*(dyfoup(i,k,2*NNy+j)+dzfoup(i,k,NNy+j)) dzfoup(i,k,2*NNy+j)];
      
      omr =[(dyfoup(i,k,2*NNy+j)-dzfoup(i,k,NNy+j)) (dzfoup(i,k,j)-dxfoup(i,k,2*NNy+j)) (dxfoup(i,k,NNy+j)-dyfoup(i,k,j))];

      [eigvec,eigval]=eig(sr);
%      eigval
      lc=eigvec(:,1);
      l2=eigvec(:,2);
      ls=eigvec(:,3);
	lc=lc/norm(lc);
	l2=l2/norm(l2);
	ls=ls/norm(ls);
	omr=omr/norm(omr);
      
      l2=l2*sign(dot(l2,omr));
	lc=lc*sign(lc(1));
    	ls=ls*sign(dot(ls,cross(l2,lc)));  
            
      pdftempc(i,k,j)=dot(lc,omr);
      pdftemp2(i,k,j)=dot(l2,omr);
      pdftemps(i,k,j)=dot(ls,omr);
      
         
   pdfytemplcx(i,k,j)=lc(dim);
   pdfytempl2x(i,k,j)=l2(dim);
   pdfytemplsx(i,k,j)=ls(dim);
    
    end
   end
end
matlabpool close




% dudyw=sum(sum(abs(dyfoup(:,:,1)-dyfoup(:,:,NNy))))/(NNx*NNz*2);
% Retau=sqrt(Re*dudyw);



%matlabpool(nproc) 
for j=1:NNy
        
            pylcx=hist(reshape(pdfytemplcx(:,:,j),NNx*NNz,1),lb:step:ub);
            pyl2x=hist(reshape(pdfytempl2x(:,:,j),NNx*NNz,1),lb:step:ub);
            pylsx=hist(reshape(pdfytemplsx(:,:,j),NNx*NNz,1),lb:step:ub);
            
            pdfylcx(:,j)=pylcx(:)/((ub-lb)*mean(pylcx));
            pdfyl2x(:,j)=pyl2x(:)/((ub-lb)*mean(pyl2x));
            pdfylsx(:,j)=pylsx(:)/((ub-lb)*mean(pylsx));
            

end

%matlabpool close


        pdfc=hist(reshape(pdftempc,NNx*NNz*NNy,1),lb:step:ub);
        pdf2=hist(reshape(pdftemp2,NNx*NNz*NNy,1),lb:step:ub);
        pdfs=hist(reshape(pdftemps,NNx*NNz*NNy,1),lb:step:ub);

pc=pdfc/((ub-lb)*mean(pdfc));
p2=pdf2/((ub-lb)*mean(pdf2));
ps=pdfs/((ub-lb)*mean(pdfs));

for j=1:NNy
        temp1=hist(reshape(pdftempc(:,:,j),NNx*NNz,1),lb:step:ub);
        temp2=hist(reshape(pdftemp2(:,:,j),NNx*NNz,1),lb:step:ub);
        temp3=hist(reshape(pdftemps(:,:,j),NNx*NNz,1),lb:step:ub);
           pcy(:,j)=temp1(:)/((ub-lb)*mean(temp1));
            p2y(:,j)=temp2(:)/((ub-lb)*mean(temp2));
            psy(:,j)=temp3(:)/((ub-lb)*mean(temp3));
	
end    



QWmean=sum(QWtemp)/(NNx*NNy*NNz);





