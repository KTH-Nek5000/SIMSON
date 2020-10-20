% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [uflocc,uflocs,ufloclam2,locx,locy,locz, Q,R, QWm,Qvf,Rvf,QWvfm]=comp_eig_ser(u,v,w,xF,yF,zF,kxvec,kzvec,ex,ez,locgrid,locbox,boxsize,skipstep,nproc,xfour,xfd,zfd);
%
% Function that computes the eigenvalues and eigenvectors of the strain
% rate tensor. Local principal coordinates are determined.
% Resample the fluctuating velocity field in the local coordinates.
%
% Resample the velocity field in the local coordinates and compute Q R to
% get Q R scatter plot in the local coordinates

%
% Transform velocity to physical space
%

%[phys,NNx,NNy,NNz]=fou2phys(fou,0,0);
% sp=size(phys);
% NNx=sp(1);
% NNz=sp(2);
% NNy=sp(3)/3;

% now phys is the physical-space velocity
% NNx : grid points in x (one point shorter than xF)
% ...


%  u=phys(:,:,1+0*NNy:1*NNy);
%  v=phys(:,:,1+1*NNy:2*NNy);
%  w=phys(:,:,1+2*NNy:3*NNy);
% clearvars phys
sp=size(u);
NNx=sp(1);
NNz=sp(2);
NNy=sp(3);
fprintf('The cut global Volume size: x=%d y=%d z=%d',NNx,NNy,NNz);
if xfour==1
%%%%%%%% for homogeneous x direction %%%%%%%%%%%%%%%%%%%
for j=1:NNy
    temp=[0 0 0];
   for k=1:NNz
    for i=1:NNx
        temp(1)=temp(1)+u(i,k,j);
        temp(2)=temp(2)+v(i,k,j);
        temp(3)=temp(3)+w(i,k,j);  
    end
   end
   for ll=1:3
     umean(j,ll)=temp(ll)/(NNx*NNz);
   end
end

    
for j=1:NNy
        ufx(:,:,j)=u(:,:,j)-umean(j,1);
        ufy(:,:,j)=v(:,:,j)-umean(j,2);
        ufz(:,:,j)=w(:,:,j)-umean(j,3);
 end

else
%%%%%%%%%% for nonhomogeneous x direction %%%%%%%%%%%%%%%%%
for j=1:NNy
    for i=1:NNx
      temp=[0 0 0];
  	 for k=1:NNz

        temp(1)=temp(1)+u(i,k,j);
        temp(2)=temp(2)+v(i,k,j);
        temp(3)=temp(3)+w(i,k,j);  
    	end
   	for ll=1:3
     	umean(i,j,ll)=temp(ll)/NNz;
   	end
   end
end

    
for j=1:NNy
    for i=1:NNx
        ufx(i,:,j)=u(i,:,j)-umean(i,j,1);
        ufy(i,:,j)=v(i,:,j)-umean(i,j,2);
        ufz(i,:,j)=w(i,:,j)-umean(i,j,3);

   end
 end
end



%%%%%%%% section for testing a fraction of y start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%% section for testing a fraction of y end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% if 1

%fprintf('The maximum relative difference between two methods: %f  \n',maxdiff);

% QWtemp=zeros(NNy,1);
% Q=zeros(NNx,NNz,NNy);
% R=zeros(NNx,NNz,NNy);
% parfor j=1:NNy
%    for k=1:NNz
%     for i=1:NNx
% %       sr =[dxfoup(i,k,j) .5*(dyfoup(i,k,j)+dxfoup(i,k,NNy+j)) .5*(dzfoup(i,k,j)+dxfoup(i,k,2*NNy+j));
% %           .5*(dxfoup(i,k,NNy+j)+dyfoup(i,k,j)) dyfoup(i,k,NNy+j) .5*(dzfoup(i,k,NNy+j)+dyfoup(i,k,2*NNy+j));
% %           .5*(dxfoup(i,k,2*NNy+j)+dzfoup(i,k,j)) .5*(dyfoup(i,k,2*NNy+j)+dzfoup(i,k,NNy+j)) dzfoup(i,k,2*NNy+j)];
% 
%       vgt =[dxfoup(i,k,j) dyfoup(i,k,j) dzfoup(i,k,j);
%             dxfoup(i,k,NNy+j) dyfoup(i,k,NNy+j) dzfoup(i,k,NNy+j);
%             dxfoup(i,k,2*NNy+j) dyfoup(i,k,2*NNy+j) dzfoup(i,k,2*NNy+j)];
%       or =[0. .5*(vgt(1,2)-vgt(2,1)) .5*(vgt(1,3)-vgt(3,1));
%           .5*(vgt(2,1)-vgt(1,2)) 0. .5*(vgt(2,3)-vgt(3,2));
%           .5*(vgt(3,1)-vgt(1,3)) .5*(vgt(3,2)-vgt(2,3)) 0.];
%       
% %        Q(i,k,j)=0;
% %        R(i,k,j)=0;
%         for ll=1:3
%             for mm=1:3
%                 Q(i,k,j) = Q(i,k,j) + (-0.5)*vgt(ll,mm)*vgt(mm,ll);
%                 QWtemp(j) = QWtemp(j) + (-0.5)*or(ll,mm)*or(mm,ll);
%                 for nn=1:3
%                     R(i,k,j) = R(i,k,j) + (-1./3.)*vgt(ll,mm)*vgt(mm,nn)*vgt(nn,ll);
%                 end
%             end
%         end
%         
%     
%     end
%    end
% end
% 
% % QWmean=0.;
% % for j=1:NNy
% %  QWmean=QWmean+QWtemp(j)/(NNx*NNy*NNz);
% % end
% QWmean=sum(QWtemp)/(NNx*NNy*NNz);


%retau = 180;
%locgrid = 0.01;
locnum = locbox+1;
locmid =(locnum+1)/2;


locv = -(locmid-1)*locgrid:locgrid:(locmid-1)*locgrid;

[locx, locy, locz]=meshgrid(locv); 


xs=boxsize(1);
xe=boxsize(2);
ys=boxsize(3);
ye=boxsize(4);
zs=boxsize(5);
ze=boxsize(6);

nyloc=floor((ye-ys)/skipstep(2))+1;
%num=zeros(1,nyloc);
ufnum=zeros(locnum*locnum*locnum,1);
unum=zeros(locnum*locnum*locnum,1);

uflocctemp=zeros(locnum*locnum*locnum,1);
uflocstemp=zeros(locnum*locnum*locnum,1);
ufloclam2temp=zeros(locnum*locnum*locnum,1);
ulocctemp=zeros(locnum*locnum*locnum,1);
ulocstemp=zeros(locnum*locnum*locnum,1);
uloclam2temp=zeros(locnum*locnum*locnum,1);
% qloctemp=zeros(locnum*locnum*locnum,nyloc);
% rloctemp=zeros(locnum*locnum*locnum,nyloc);


%sr=zeros(3,3);

%
% Running parallel on nproc processes
%
 try
     matlabpool(nproc)
 catch
     matlabpool close
     matlabpool(nproc)
 end
for jj=1:nyloc
    j=ys+(jj-1)*skipstep(2)
   for k=zs:skipstep(3):ze
	
    for i=xs:skipstep(1):xe
        


      
      sr =[dxfoup(i,k,j) .5*(dyfoup(i,k,j)+dxfoup(i,k,NNy+j)) .5*(dzfoup(i,k,j)+dxfoup(i,k,2*NNy+j));
          .5*(dxfoup(i,k,NNy+j)+dyfoup(i,k,j)) dyfoup(i,k,NNy+j) .5*(dzfoup(i,k,NNy+j)+dyfoup(i,k,2*NNy+j));
          .5*(dxfoup(i,k,2*NNy+j)+dzfoup(i,k,j)) .5*(dyfoup(i,k,2*NNy+j)+dzfoup(i,k,NNy+j)) dzfoup(i,k,2*NNy+j)];

      
      omr =[(dyfoup(i,k,2*NNy+j)-dzfoup(i,k,NNy+j)) (dzfoup(i,k,j)-dxfoup(i,k,2*NNy+j)) (dxfoup(i,k,NNy+j)-dyfoup(i,k,j))];

      [eigvec,eigval]=eig(sr);
%
%      lc: most compressing direction; 
%      l2: lambda2
%      ls: most stretching direction
%
      lc=eigvec(:,1);
      l2=eigvec(:,2);
      ls=eigvec(:,3);

 
%      if dot(l2,omr) < 0
%          l2= -l2;
%      end
%       if lc(1) < 0
%           lc = -lc;
%       end
%       if dot(ls,cross(lc,l2)) < 0
%           ls = -ls;
%       end
	l2=l2/norm(l2);
	lc=lc/norm(lc);
	l2=l2*sign(dot(l2,omr));
%	lc=lc*sign(lc(1));
	ls=cross(l2,lc);
%	lc=lc/norm(lc);
%      ls=cross(l2,lc);
      

      lv=inv([lc ls l2]);
      
      locmidx=xF(i);
      locmidy=yF(j); 
      locmidz=zF(k);
           
      loc2gx=locmidx+locx(:)*lc(1)+locy(:)*ls(1)+locz(:)*l2(1);
      loc2gy=locmidy+locx(:)*lc(2)+locy(:)*ls(2)+locz(:)*l2(2);
      loc2gz=locmidz+locx(:)*lc(3)+locy(:)*ls(3)+locz(:)*l2(3);
      if zfd==0
 	loc2gz=mod(loc2gz-zF(1),zF(NNz+1)-zF(1))+zF(1);
      end
%      fprintf('max loc2gx=%f, loc2gy=%f,loc2gz=%f \n',max(loc2gx),max(loc2gy),max(loc2gz));
%       for ll=1:locnum
%           for mm = 1:locnum
%               for nn = 1:locnum
%                   loc2gx(ll,mm,nn)=locmidx+locx(ll,mm,nn)*dot(lc,[1 0 0])+locy(ll,mm,nn)*dot(ls,[1 0 0])+locz(ll,mm,nn)*dot(l2,[1 0 0]);
%                   loc2gy(ll,mm,nn)=locmidy+locx(ll,mm,nn)*dot(lc,[0 1 0])+locy(ll,mm,nn)*dot(ls,[0 1 0])+locz(ll,mm,nn)*dot(l2,[0 1 0]);
%                   loc2gz(ll,mm,nn)=locmidz+locx(ll,mm,nn)*dot(lc,[0 0 1])+locy(ll,mm,nn)*dot(ls,[0 0 1])+locz(ll,mm,nn)*dot(l2,[0 0 1]);
%               end
%           end
%       end
%       
 

        [z,x, y] = meshgrid(zF(1:NNz),xF(1:NNx),yF(1:NNy));



%
% 3D interpretation of the fields
%
       ufloc2g1 =interp3(z,x,y,ufx,loc2gz,loc2gx,loc2gy,'linear',0);
       ufloc2g2 =interp3(z,x,y,ufy,loc2gz,loc2gx,loc2gy,'linear',0);
       ufloc2g3 =interp3(z,x,y,ufz,loc2gz,loc2gx,loc2gy,'linear',0);  
 %       fprintf('min ufx=%f  ufloc2g1=%f,  \n',min(ufx),min(ufloc2g1));     
       uloc2g1 =interp3(z,x,y,u(1:NNx,:,1:NNy),loc2gz,loc2gx,loc2gy,'linear',0);
       uloc2g2 =interp3(z,x,y,v(1:NNx,:,1:NNy),loc2gz,loc2gx,loc2gy,'linear',0);
       uloc2g3 =interp3(z,x,y,w(1:NNx,:,1:NNy),loc2gz,loc2gx,loc2gy,'linear',0);  
       
%        qtemp=interp3(z,x,y,Q,loc2gz,loc2gx,loc2gy,'linear',0);
%        rtemp=interp3(z,x,y,R,loc2gz,loc2gx,loc2gy,'linear',0);
%        
%        qloctemp(:,jj) =qloctemp(:,jj)+qtemp(:);
%        rloctemp(:,jj) =rloctemp(:,jj)+rtemp(:);

%       dloctemp(;,jj)=dloctemp(;,jj)+qtemp^3+(27./4.)*rtemp^2;
       
%       qloc=qloctemp/QWmean;
%       rloc=rloctemp/(QWmean^1.5);
 

%       for ll=1:locnum
%           for mm = 1:locnum
%               for nn = 1:locnum              
%                 uflc(ll,mm,nn)=lv(1,1)*ufloc2g1(ll,mm,nn)+lv(1,2)*ufloc2g2(ll,mm,nn)+lv(1,3)*ufloc2g3(ll,mm,nn);
%                 ufls(ll,mm,nn)=lv(2,1)*ufloc2g1(ll,mm,nn)+lv(2,2)*ufloc2g2(ll,mm,nn)+lv(2,3)*ufloc2g3(ll,mm,nn);
%                 uflam2(ll,mm,nn)=lv(3,1)*ufloc2g1(ll,mm,nn)+lv(3,2)*ufloc2g2(ll,mm,nn)+lv(3,3)*ufloc2g3(ll,mm,nn);         
%               
%               end
%           end
%       end   
%tc2=cputime-tc      

%
% get the three components of the (fluctuating) velocity vector in the local coordinates
%
                uflc=lv(1,1)*ufloc2g1(:) +lv(1,2)*ufloc2g2(:) +lv(1,3)*ufloc2g3(:); 
                ufls=lv(2,1)*ufloc2g1(:)+lv(2,2)*ufloc2g2(:)+lv(2,3)*ufloc2g3(:);
                uflam2=lv(3,1)*ufloc2g1(:)+lv(3,2)*ufloc2g2(:)+lv(3,3)*ufloc2g3(:) ; 
                
                ulc=lv(1,1)*uloc2g1(:) +lv(1,2)*uloc2g2(:) +lv(1,3)*uloc2g3(:); 
                uls=lv(2,1)*uloc2g1(:)+lv(2,2)*uloc2g2(:)+lv(2,3)*uloc2g3(:);
                ulam2=lv(3,1)*uloc2g1(:)+lv(3,2)*uloc2g2(:)+lv(3,3)*uloc2g3(:) ; 
                
                uflocctemp=uflocctemp+uflc;
                uflocstemp=uflocstemp+ufls;
                ufloclam2temp=ufloclam2temp+uflam2;
                
                ulocctemp=ulocctemp+ulc;
                ulocstemp=ulocstemp+uls;
                uloclam2temp=uloclam2temp+ulam2;
      

		ufnum=ufnum+double(sqrt(uflc.^2+ufls.^2+uflam2.^2)>0);
		unum=unum+double(sqrt(ulc.^2+uls.^2+ulam2.^2)>0);
%            if sqrt(uflc.^2+ufls.^2+uflam2.^2)>0
%                num(jj)=num(jj)+1;
%        fprintf(' num=%d    \n',num(jj));    
%            end

%totalnum=totalnum+1;
      
    end
  end
end
matlabpool close



%totalnum=sum(num);

if ufnum==0
ufnum=1;
end
if unum==0
unum=1;
end
%fprintf('ufnum=  unum=\n');
%reshape(ufnum,size(locx))
%reshape(unum,size(locx))
       uflocc = uflocctemp./ufnum;
       uflocs = uflocstemp./ufnum;
       ufloclam2 = ufloclam2temp./ufnum;
       ulocc = ulocctemp./unum;
       ulocs = ulocstemp./unum;
       uloclam2 = uloclam2temp./unum;
       

       
    uflocc=reshape(uflocc,size(locx));
   uflocs=reshape(uflocs,size(locx));
   ufloclam2=reshape(ufloclam2,size(locx));
   ulocc=reshape(ulocc,size(locx));
   ulocs=reshape(ulocs,size(locx));
   uloclam2=reshape(uloclam2,size(locx));
   
 %
 % Compute gradient of velocity field int the local coordinates
 %
 
%  if matgrad==1
%  [vgl(:,:,:,1,1), vgl(:,:,:,1,2), vgl(:,:,:,1,3)]=gradient(ulocc,locgrid); 
%  [vgl(:,:,:,2,1), vgl(:,:,:,2,2), vgl(:,:,:,2,3)]=gradient(ulocs,locgrid);
%  [vgl(:,:,:,3,1), vgl(:,:,:,3,2), vgl(:,:,:,3,3)]=gradient(uloclam2,locgrid);
%  else
     [mloc]= fds(1, locnum, 4, 0,locv);
     for k=1:locnum
  
                vgl(:,:,k,1,2)=mloc*ulocc(:,:,k);
                vgl(:,:,k,2,2)=mloc*ulocs(:,:,k);
                vgl(:,:,k,3,2)=mloc*uloclam2(:,:,k);
                
                vgl(:,:,k,1,1)=(mloc*ulocc(:,:,k)')';
                vgl(:,:,k,2,1)=(mloc*ulocs(:,:,k)')';
                vgl(:,:,k,3,1)=(mloc*uloclam2(:,:,k)')';

                vfgl(:,:,k,1,2)=mloc*uflocc(:,:,k);
                vfgl(:,:,k,2,2)=mloc*uflocs(:,:,k);
                vfgl(:,:,k,3,2)=mloc*ufloclam2(:,:,k);
                
                vfgl(:,:,k,1,1)=(mloc*uflocc(:,:,k)')';
                vfgl(:,:,k,2,1)=(mloc*uflocs(:,:,k)')';
                vfgl(:,:,k,3,1)=(mloc*ufloclam2(:,:,k)')';
 
 
     end
    
     for i=1:locnum
         for j=1:locnum
                vgl(i,j,:,1,3)=mloc*reshape(ulocc(i,j,:),locnum,1);
                vgl(i,j,:,2,3)=mloc*reshape(ulocs(i,j,:),locnum,1);
                vgl(i,j,:,3,3)=mloc*reshape(uloclam2(i,j,:),locnum,1);
               vfgl(i,j,:,1,3)=mloc*reshape(uflocc(i,j,:),locnum,1);
                vfgl(i,j,:,2,3)=mloc*reshape(uflocs(i,j,:),locnum,1);
                vfgl(i,j,:,3,3)=mloc*reshape(ufloclam2(i,j,:),locnum,1);
            
         end
     end
%  end
 
 
 
 Q=zeros(locnum,locnum,locnum);
 R=zeros(locnum,locnum,locnum);
Qvf=zeros(locnum,locnum,locnum);
 Rvf=zeros(locnum,locnum,locnum);
 
 QWtemp=zeros(locnum,locnum,locnum);
 QWvftemp=zeros(locnum,locnum,locnum);
 for l=1:3
     for m=1:3
        Q(:,:,:)=Q(:,:,:)-0.5*(vgl(:,:,:,l,m).*vgl(:,:,:,m,l));
       QWtemp(:,:,:) = QWtemp(:,:,:) - 0.125*(vgl(:,:,:,l,m)-vgl(:,:,:,m,l)).*(vgl(:,:,:,m,l)-vgl(:,:,:,l,m));
        Qvf(:,:,:)=Qvf(:,:,:)-0.5*(vfgl(:,:,:,l,m).*vfgl(:,:,:,m,l));
       QWvftemp(:,:,:) = QWvftemp(:,:,:) - 0.125*(vfgl(:,:,:,l,m)-vfgl(:,:,:,m,l)).*(vfgl(:,:,:,m,l)-vfgl(:,:,:,l,m));
  %       Q=Q-0.5*(vgl(:,l,m).*vgl(:,m,l));
         for n=1:3
             R(:,:,:)=R(:,:,:)-(1./3.)*vgl(:,:,:,l,m).*vgl(:,:,:,m,n).*vgl(:,:,:,n,l);
             Rvf(:,:,:)=Rvf(:,:,:)-(1./3.)*vfgl(:,:,:,l,m).*vfgl(:,:,:,m,n).*vfgl(:,:,:,n,l);
 %           R=R-(1./3.)*vgl(:,l,m).*vgl(:,m,n).*vgl(:,n,l);
         end
     end
 end
 
 Q=reshape(Q,locnum*locnum*locnum,1);
 R=reshape(R,locnum*locnum*locnum,1);
 QWm=sum(QWtemp(:))/(locnum*locnum*locnum);
 Qvf=reshape(Qvf,locnum*locnum*locnum,1);
 Rvf=reshape(Rvf,locnum*locnum*locnum,1);
 QWvfm=sum(QWvftemp(:))/(locnum*locnum*locnum);
 
 %end   
       
