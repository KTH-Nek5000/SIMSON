
% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
clear all
close all
filename='t1000.u';
scalar=3
xfd=1;
xfour=1;
comp_eig=0;
comp_pdf=1;

fprintf('compact finite difference scheme in the x direction xfd= %d \n',xfd );

tic;
[vel,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv,kxvec,kzvec]= readdns(filename,scalar);

t2=toc;
fprintf('read dns time: %f s \n',t2 );


% the variables are:
% vel : velocity in Fourier space, contains all velocity components
% xF  : coordinates in x (note, there is one more point at the end)
% yF  : coordinates in y (note, top is index 1 and wall last index)
% zF  : coordinates in z (note, there is one more point at the end)
%
% t   : time
% Re  : Reynolds number
% rest not important


% find the available number of processes/workers for parallel computation 
out=findResource();
nproc=get(out,'ClusterSize')
%nproc=6 


[phys,NNx,NNy,NNz]=fou2phys(vel,0,0);
Retau=retauf(phys(:,:,1:NNy),yF,kxvec,kzvec,Re);

fprintf('friction Reynolds number: Retau = %f \n', Retau);

if 1
step=0.01 %bin size
lb=-1.0 %lower bound
ub=1.0  %upper bound
dim=3 % 1: streamwise direction, 2: wall normal, 3: spanwise
if dim==1
    dimstr='x';
elseif dim==2
    dimstr='y';
elseif dim==3
    dimstr='z';
end

%Retau=180

[pc,p2,ps,pdfylcx,pdfyl2x,pdfylsx,NNy,Q,R,QWmean]=pdf_par(phys,xF,yF,zF,kxvec,kzvec,0,0,lb,ub,step,nproc,xfour,xfd,dim);

x=lb:step:ub;
save('pdfang.mat','x','pc','p2','ps','pdfylcx','pdfyl2x','pdfylsx')



%fprintf('Retau= %f \n',Retau);
F8=nearest((yF(1)-yF(8))*Retau);
F16=nearest((yF(1)-yF(16))*Retau);
F32=nearest((yF(1)-yF(32))*Retau);
F64=nearest((yF(1)-yF(64))*Retau);
FL8=nearest((yF(NNy-7)-yF(NNy))*Retau);   
%Q=Q/QWmean;
%R=R/(QWmean^1.5);
% [xx, yy]=meshgrid(floor(min(R(:))):0.1:ceil(max(R(:))), floor(min(Q(:))):0.1:ceil(max(Q(:))));

Ry=[R(:,:,8) R(:,:,16) R(:,:,32) R(:,:,64)];
Qy=[Q(:,:,8) Q(:,:,16) Q(:,:,32) Q(:,:,64)];
[xx, yy]=meshgrid(floor(min(Ry(:))):0.1:ceil(max(Ry(:))), floor(min(Qy(:))):0.1:ceil(max(Qy(:))));
 DD=(27./4.)*(xx.^2)+yy.^3;



if comp_pdf==1
    
pdf=figure();
pfig1=plot(x,ps,'-');
hold on
pfig2=plot(x,p2,'--');
hold on
pfig3=plot(x,pc,'-.');
xlim([lb-step ub+step]);
xlabel('cos(\theta_{\omega,i})')
ylabel('p.d.f.')
legend([pfig1 pfig2 pfig3], {'\theta_{\omega,1}','\theta_{\omega,2}','\theta_{\omega,3}'})
saveas(pdf,strcat('pdfangle_FD',num2str(xfd),'.eps'),'eps');

pdfyl1=figure();

    pyfigx1=plot(x,pdfylsx(:,8),'-');
    hold on
    pyfigx2=plot(x,pdfylsx(:,16),'--');
    hold on
    pyfigx3=plot(x,pdfylsx(:,32),'-.');
    hold on
    pyfigx4=plot(x,pdfylsx(:,64),':');
    hold on
    pyfigx5=plot(x,pdfylsx(:,NNy-7),'+');

    xlabel(strcat('cos(\theta_{\lambda_1,',dimstr,'})'))
    ylabel('p.d.f.')

    legend([pyfigx1 pyfigx2 pyfigx3 pyfigx4 pyfigx5], {strcat('y^+=',num2str(F8)),strcat('y^+=',num2str(F16)),strcat('y^+=',num2str(F32)),strcat('y^+=',num2str(F64)),strcat('Lw y^+=',num2str(FL8))})
saveas(pdfyl1,strcat('pdfl1',dimstr,'_FD',num2str(xfd),'.eps'),'eps');
    

    pdfyl2=figure();
    pyfigx1=plot(x,pdfyl2x(:,8),'-');
    hold on
    pyfigx2=plot(x,pdfyl2x(:,16),'--');
    hold on
    pyfigx3=plot(x,pdfyl2x(:,32),'-.');
    hold on
    pyfigx4=plot(x,pdfyl2x(:,64),':');
        hold on
    pyfigx5=plot(x,pdfyl2x(:,NNy-7),'+');


    xlabel(strcat('cos(\theta_{\lambda_2,',dimstr,'})'))
    ylabel('p.d.f.')

       legend([pyfigx1 pyfigx2 pyfigx3 pyfigx4 pyfigx5], {strcat('y^+=',num2str(F8)),strcat('y^+=',num2str(F16)),strcat('y^+=',num2str(F32)),strcat('y^+=',num2str(F64)),strcat('Lw y^+=',num2str(FL8))})
    saveas(pdfyl2,strcat('pdfl2',dimstr,'_FD',num2str(xfd),'.eps'),'eps');

    pdfyl3=figure();
    pyfigx1=plot(x,pdfylcx(:,8),'-');
    hold on
    pyfigx2=plot(x,pdfylcx(:,16),'--');
    hold on
    pyfigx3=plot(x,pdfylcx(:,32),'-.');
    hold on
    pyfigx4=plot(x,pdfylcx(:,64),':');
        hold on
    pyfigx5=plot(x,pdfylcx(:,NNy-7),'+');


    xlabel(strcat('cos(\theta_{\lambda_1,',dimstr,'})'))
    ylabel('p.d.f.')

    legend([pyfigx1 pyfigx2 pyfigx3 pyfigx4 pyfigx5], {strcat('y^+=',num2str(F8)),strcat('y^+=',num2str(F16)),strcat('y^+=',num2str(F32)),strcat('y^+=',num2str(F64)),strcat('Lw y^+=',num2str(FL8))})
    saveas(pdfyl3,strcat('pdfl3',dimstr,'_FD',num2str(xfd),'.eps'),'eps');
end


% [xx, yy]=meshgrid(floor(min(min(R(:,:,8)))):0.1:ceil(max(max(R(:,:,8)))), floor(min(min(Q(:,:,8)))):0.1:ceil(max(max(Q(:,:,8)))));
%DD=(27./4.)*(xx.^2)+yy.^3;
    pdfqr1=figure();
plot(R(:,:,8),Q(:,:,8),'.');
hold on
contour(xx,yy,DD,[0. 0.])
    xlabel('R')
    ylabel('Q')
    title(strcat('Q R scatter plot at y^+=',num2str(F8)));

   saveas(pdfqr1,strcat('qr_yp',num2str(F8),'_FD',num2str(xfd),'.eps'),'psc2');
   hold off
   
   
%   [xx, yy]=meshgrid(floor(min(min(R(:,:,16)))):0.1:ceil(max(max(R(:,:,16)))), floor(min(min(Q(:,:,16)))):0.1:ceil(max(max(Q(:,:,16)))));
% DD=(27./4.)*(xx.^2)+yy.^3;
      pdfqr2=figure();
plot(R(:,:,16),Q(:,:,16),'.');
hold on
contour(xx,yy,DD,[0. 0.])

    xlabel('R')
    ylabel('Q')
    title(strcat('Q R scatter plot at y^+=',num2str(F16)));

   saveas(pdfqr2,strcat('qr_yp',num2str(F16),'_FD',num2str(xfd),'.eps'),'psc2');
    hold off 
   
% [xx, yy]=meshgrid(floor(min(min(R(:,:,32)))):0.1:ceil(max(max(R(:,:,32)))), floor(min(min(Q(:,:,32)))):0.1:ceil(max(max(Q(:,:,32)))));
% DD=(27./4.)*(xx.^2)+yy.^3;  
   pdfqr3=figure();
plot(R(:,:,32),Q(:,:,32),'.');
hold on
contour(xx,yy,DD,[0. 0.])
    xlabel('R')
    ylabel('Q')
    title(strcat('Q R scatter plot at y^+=',num2str(F32)));

   saveas(pdfqr3,strcat('qr_yp',num2str(F32),'_FD',num2str(xfd),'.eps'),'psc2');
  hold off
  
%  [xx, yy]=meshgrid(floor(min(min(R(:,:,64)))):0.1:ceil(max(max(R(:,:,64)))), floor(min(min(Q(:,:,64)))):0.1:ceil(max(max(Q(:,:,64)))));
% DD=(27./4.)*(xx.^2)+yy.^3; 

      pdfqr4=figure();
plot(R(:,:,64),Q(:,:,64),'.');
hold on
contour(xx,yy,DD,[0. 0.])

    xlabel('R')
    ylabel('Q')
%      set(gca, 'XTick', 0.,'YTick', 0.);
    title(strcat('Q R scatter plot at y^+=',num2str(F64)));

   saveas(pdfqr4,strcat('qr_yp',num2str(F64),'_FD',num2str(xfd),'.eps'),'psc2');
   hold off
   
%     
%    pdfqr1=figure();
% plot(R(:,:,8),Q(:,:,8),'.');
% % hold on
% % qrd=plot(rloc,qloc,);
%     xlabel('R/<Qw>^{1.5}')
%     ylabel('Q/<Qw>')
%     title(strcat('Q R scatter plot at y^+=',num2str(F8)));
% 
%    saveas(pdfqr1,'qr1.eps','eps');
%   
%       pdfqr2=figure();
% plot(R(:,:,16),Q(:,:,16),'.');
% % hold on
% % qrd=plot(rloc,qloc,);
%     xlabel('R/<Qw>^{1.5}')
%     ylabel('Q/<Qw>')
%     title(strcat('Q R scatter plot at y^+=',num2str(F16)));
% 
%    saveas(pdfqr2,'qr2.eps','eps');
%    
%    
%   
%    pdfqr3=figure();
% plot(R(:,:,32),Q(:,:,32),'.');
% % hold on
% % qrd=plot(rloc,qloc,);
%     xlabel('R/<Qw>^{1.5}')
%     ylabel('Q/<Qw>')
%     title(strcat('Q R scatter plot at y^+=',num2str(F32)));
% 
%    saveas(pdfqr3,'qr3.eps','eps');
%   
%       pdfqr4=figure();
% plot(R(:,:,64),Q(:,:,64),'.');
% % hold on
% % qrd=plot(rloc,qloc,);
%     xlabel('R/<Qw>^{1.5}')
%     ylabel('Q/<Qw>')
%     title(strcat('Q R scatter plot at y^+=',num2str(F64)));
% 
%    saveas(pdfqr4,'qr4.eps','eps');
%   
% 
%     
   

end
 t3=toc-t2;
 
fprintf('pdf computation time: %f s \n',t3);   

if comp_eig==1

%local uniform grid resolution and box size
locgrid=0.02
locboxsize=0.8
locbox=locboxsize/locgrid
%global measurement volume start and end location (x y z respectively)
nl=50;
nh=78;
boxsize=[nl nh nl nh nl nh]
%the number of grid points to be skipped for the averaging 
skipstep=4
locmid=locbox/2+1;

[uflocc,uflocs,ufloclam2,locx,locy,locz,locv,totalnum, qloc,rloc, QWm]=comp_eig_par(phys,xF,yF,zF,kxvec,kzvec,0,0,locgrid,locbox,boxsize,skipstep,nproc,xfour,xfd);
%loc1
%uflocc 
if 1
qloc=qloc/QWm;
rloc=rloc/(QWm^1.5);
[xr, yq]=meshgrid(floor(min(rloc(:))):0.1:ceil(max(rloc(:))), floor(min(qloc(:))):0.1:ceil(max(qloc(:))));
Dqr=(27./4.)*(xr.^2)+yq.^3;
%Dqr=((27./4.)*(xr.^2)+yq.^3)*QWm^3;

locx=locx*Retau;
locy=locy*Retau;
locz=locz*Retau;

t4=toc-t3;
fprintf('eigen computation time: %f s \n',t4);

save('ufchanbig.mat','locx','locy','locz','uflocc','uflocs','ufloclam2','totalnum')



 fig1=figure();
 quiver(locx(:,:,locmid),locy(:,:,locmid),uflocc(:,:,locmid),uflocs(:,:,locmid))
    xlabel('\lambda_3')
    ylabel('\lambda_1')
 title('\lambda_2=0')
 saveas(fig1,strcat('ufloc.l2.0_',num2str(nl),'_',num2str(nh),'_loc',num2str(locbox),'_FD',num2str(xfd),'.eps'),'eps');


 fig2=figure();
 quiver(squeeze(locy(:,locmid,:)),squeeze(locz(:,locmid,:)),squeeze(uflocs(:,locmid,:)),squeeze(ufloclam2(:,locmid,:)))
    xlabel('\lambda_1')
    ylabel('\lambda_2')
 title('\lambda_3=0')
  saveas(fig2,strcat('ufloc.l3.0_',num2str(nl),'_',num2str(nh),'_loc',num2str(locbox),'_FD',num2str(xfd),'.eps'),'eps');


 fig3=figure();
 quiver(squeeze(locx(locmid,:,:)),squeeze(locz(locmid,:,:)),squeeze(uflocc(locmid,:,:)),squeeze(ufloclam2(locmid,:,:)))
    xlabel('\lambda_3')
    ylabel('\lambda_2')
 title('\lambda_1=0')
  saveas(fig3,strcat('ufloc.l1.0_',num2str(nl),'_',num2str(nh),'_loc',num2str(locbox),'_FD',num2str(xfd),'.eps'),'eps');

  
  
fig4=figure();
qr=plot(rloc,qloc,'.');
hold on
dqrfig=contour(xr,yq,Dqr,[0. 0.]);
    xlabel('R/<Qw>^{1.5}')
    ylabel('Q/<Qw>')
    title('Q R scatter plot in the strain rate tensor eigenframe');

   saveas(fig4,strcat('qr_',num2str(nl),'_',num2str(nh),'_loc',num2str(locbox),'_FD',num2str(xfd),'.eps'),'psc2');
  

fprintf('visualization time: %f s\n',toc-t4);

end
end
%whos

fprintf('Total time: %f s\n',toc);





