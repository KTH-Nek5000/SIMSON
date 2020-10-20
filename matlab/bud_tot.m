clear all
%close all


% Scalings:
% lstar=nu/utau
% Budget: inner: utau^3/lstar = utau^4/nu
%         
%         outer: utau^3/delta
%                multiply with y/delta
%                ---> f = y/delta * delta/utau^3 = y/utau^3



%v=load('vel3293.data.25');
%b=load('budtot3293.data.25');   
%v=load('vel1342.data.10');
%b=load('budtot1342.data.10');   
v=load('vel0688.data');
b=load('budtot0688.data');   
%v=load('vel0171.data');
%b=load('budtot0171.data');   


% Re_tau=1000: 3293.25
% Re_tau=550:  1342.10
% Re_tau=360:   688
% Re_tau=180:   171

re = 450;

y=v(:,1);
ny=length(v);
yl = v(ny,1);

[re_theta,re_deltastar,utau,lstar,H12] = comp_re(yl,v(:,2),re);
[res,deltas,utau,lstar,H12] = comp_re1(yl,v(:,2),re);
disp(sprintf('re_theta = %f re_tau=%f',re_theta,deltas(1)/lstar));


% inner scaling
fact  = utau^3/lstar;
factl = lstar;
% outer scaling
%factl = re_deltastar/re;
%fact  = 1/factl;

yy = 1;%b(:,1)/factl;

b(:,2)=-b(:,2);

figure(1)
title('k budget')
hold on


plot(b(:,1)/factl,b(:,2)/fact.*yy,'r')
plot(b(:,1)/factl,b(:,3)/fact.*yy,'k')
plot(b(:,1)/factl,b(:,4)/fact.*yy,'c')
plot(b(:,1)/factl,b(:,5)/fact.*yy,'g')
plot(b(:,1)/factl,b(:,6)/fact.*yy,'m')
plot(b(:,1)/factl,b(:,7)/fact.*yy,'b')

sum1(:,2) =+b(:,2)+b(:,3)+b(:,4)+b(:,5)+b(:,6)+b(:,7);
plot(b(:,1)/factl,sum1(:,2)/fact.*yy,'y')



xlabel('\ity^+')
ylabel('loss   (inner scaling)   gain')

legend('convection','production','dissipation','turbdiff','velp','visdiff','sum')

%-------compute integrals-----------
N=length(y);
Ly=y(end);
[y,D] = chebdif(N,1);
D=D*(-2/Ly);

DD=D;DD(1,:)=0;DD(1,1)=1;iDD=inv(DD);

a1=iDD*b(:,2);s1=a1(end)-a1(1);
a1=iDD*max(b(:,2),0);s1p=a1(end)-a1(1);
a1=iDD*min(b(:,2),0);s1n=a1(end)-a1(1);
a1=iDD*b(:,3);s2=a1(end)-a1(1);
a1=iDD*b(:,4);s3=a1(end)-a1(1);
a1=iDD*b(:,5);s4=a1(end)-a1(1);
a1=iDD*b(:,6);s5=a1(end)-a1(1);
a1=iDD*b(:,7);s6=a1(end)-a1(1);

disp('Turbulent kinetic energy k')
disp('Relevant source terms; value ; value/prod')
disp(sprintf('Prod: %f15  %f15',s2,abs(s2/s2)))
disp(sprintf('Diss: %f15  %f15',s3,abs(s3/s2)))
disp(sprintf('Conv: %f15  %f15',s1,abs(s1/s2)))
disp(sprintf('   p: %f15  %f15',s1p,abs(s1p/s2)))
disp(sprintf('   n: %f15  %f15',s1n,abs(s1n/s2)))
disp(sprintf('Sum:  %f15  %f15',s1+s2+s3+s4+s5+s6,abs(s1+s2+s3+s4+s5+s6)/s2))
disp('---------------------------------')
%--------------------------




figure(2)
title('uu budget')
hold on

b(:,8) = -b(:,8);

plot(b(:,1)/factl,b(:,8)/fact.*yy,'r')
plot(b(:,1)/factl,b(:,9)/fact.*yy,'k')
plot(b(:,1)/factl,b(:,10)/fact.*yy,'c')
plot(b(:,1)/factl,b(:,11)/fact.*yy,'g')
plot(b(:,1)/factl,b(:,12)/fact.*yy,'m')
plot(b(:,1)/factl,b(:,13)/fact.*yy,'m--')
plot(b(:,1)/factl,b(:,14)/fact.*yy,'b')

sum2(:,2) =+b(:,8)+b(:,9)+b(:,10)+b(:,11)+b(:,12)+b(:,13)+b(:,14);
plot(b(:,1)/factl,sum2(:,2)/fact.*yy,'y')



xlabel('\ity^+')
ylabel('loss   (inner scaling)   gain')

legend('convection','production','dissipation','turbdiff','pdiff','pstrain','visdiff','sum')


%-------compute integrals-----------
a1=iDD*b(:,8);s1=a1(end)-a1(1);
a1=iDD*max(b(:,8),0);s1p=a1(end)-a1(1);
a1=iDD*min(b(:,8),0);s1n=a1(end)-a1(1);
a1=iDD*b(:,9);s2=a1(end)-a1(1);
a1=iDD*b(:,10);s3=a1(end)-a1(1);
a1=iDD*b(:,11);s4=a1(end)-a1(1);
a1=iDD*b(:,12);s5=a1(end)-a1(1);
a1=iDD*b(:,13);s6=a1(end)-a1(1);
a1=iDD*b(:,14);s7=a1(end)-a1(1);

disp('Reynolds stress uu')
disp('Relevant source terms; value ; value/prod')
disp(sprintf('Prod: %f15  %f15',s2,abs(s2/s2)))
disp(sprintf('Diss: %f15  %f15',s3,abs(s3/s2)))
disp(sprintf('Conv: %f15  %f15',s1,abs(s1/s2)))
disp(sprintf('   p: %f15  %f15',s1p,abs(s1p/s2)))
disp(sprintf('   n: %f15  %f15',s1n,abs(s1n/s2)))
disp(sprintf('Pstr: %f15  %f15',s6,abs(s6/s2)))
disp(sprintf('Sum:  %f15  %f15',s1+s2+s3+s4+s5+s6+s7,abs(s1+s2+s3+s4+s5+s6+s7)/s2))
disp('---------------------------------')
%--------------------------

figure(3)
title('vv budget')
hold on

b(:,15) = -b(:,15);

plot(b(:,1)/factl,b(:,15)/fact.*yy,'r')
plot(b(:,1)/factl,b(:,16)/fact.*yy,'k')
plot(b(:,1)/factl,b(:,17)/fact.*yy,'c')
plot(b(:,1)/factl,b(:,18)/fact.*yy,'g')
plot(b(:,1)/factl,b(:,19)/fact.*yy,'m')
plot(b(:,1)/factl,b(:,20)/fact.*yy,'m--')
plot(b(:,1)/factl,b(:,21)/fact.*yy,'b')

sum3(:,2) =+b(:,15)+b(:,16)+b(:,17)+b(:,18)+b(:,19)+b(:,20)+b(:,21);
plot(b(:,1)/factl,sum3(:,2)/fact.*yy,'y')



xlabel('\ity^+')
ylabel('loss   (inner scaling)   gain')

legend('convection','production','dissipation','turbdiff','pdiff','pstrain','visdiff','sum')


%-------compute integrals-----------
a1=iDD*b(:,15);s1=a1(end)-a1(1);
a1=iDD*max(b(:,15),0);s1p=a1(end)-a1(1);
a1=iDD*min(b(:,15),0);s1n=a1(end)-a1(1);
a1=iDD*b(:,16);s2=a1(end)-a1(1);
a1=iDD*b(:,17);s3=a1(end)-a1(1);
a1=iDD*b(:,18);s4=a1(end)-a1(1);
a1=iDD*b(:,19);s5=a1(end)-a1(1);
a1=iDD*b(:,20);s6=a1(end)-a1(1);
a1=iDD*b(:,21);s7=a1(end)-a1(1);

disp('Reynolds stress vv')
disp('Relevant source terms; value ; value/pstrain')
disp(sprintf('Prod: %f15  %f15',s2,abs(s2/s6)))
disp(sprintf('Diss: %f15  %f15',s3,abs(s3/s6)))
disp(sprintf('Conv: %f15  %f15',s1,abs(s1/s6)))
disp(sprintf('   p: %f15  %f15',s1p,abs(s1p/s6)))
disp(sprintf('   n: %f15  %f15',s1n,abs(s1n/s6)))
disp(sprintf('Pstr: %f15  %f15',s6,abs(s6/s6)))
disp(sprintf('Sum:  %f15  %f15',s1+s2+s3+s4+s5+s6+s7,abs(s1+s2+s3+s4+s5+s6+s7)/s6))
disp('---------------------------------')
%--------------------------

figure(4)
title('ww budget')
hold on

b(:,22) = -b(:,22);
b(:,26) = 0;
 
plot(b(:,1)/factl,b(:,22)/fact.*yy,'r')
plot(b(:,1)/factl,b(:,23)/fact.*yy,'k')
plot(b(:,1)/factl,b(:,24)/fact.*yy,'c')
plot(b(:,1)/factl,b(:,25)/fact.*yy,'g')
plot(b(:,1)/factl,b(:,26)/fact.*yy,'m')
plot(b(:,1)/factl,b(:,27)/fact.*yy,'m--')
plot(b(:,1)/factl,b(:,28)/fact.*yy,'b')

sum4(:,2) =+b(:,22)+b(:,23)+b(:,24)+b(:,25)+b(:,26)+b(:,27)+b(:,28);
plot(b(:,1)/factl,sum4(:,2)/fact.*yy,'y')



xlabel('\ity^+')
ylabel('loss   (inner scaling)   gain')

legend('convection','production','dissipation','turbdiff','pdiff','pstrain','visdiff','sum')

%-------compute integrals-----------
a1=iDD*b(:,22);s1=a1(end)-a1(1);
a1=iDD*max(b(:,22),0);s1p=a1(end)-a1(1);
a1=iDD*min(b(:,22),0);s1n=a1(end)-a1(1);
a1=iDD*b(:,23);s2=a1(end)-a1(1);
a1=iDD*b(:,24);s3=a1(end)-a1(1);
a1=iDD*b(:,25);s4=a1(end)-a1(1);
a1=iDD*b(:,26);s5=a1(end)-a1(1);
a1=iDD*b(:,27);s6=a1(end)-a1(1);
a1=iDD*b(:,28);s7=a1(end)-a1(1);

disp('Reynolds stress ww')
disp('Relevant source terms; value ; value/pstrain')
disp(sprintf('Prod: %f15  %f15',s2,abs(s2/s6)))
disp(sprintf('Diss: %f15  %f15',s3,abs(s3/s6)))
disp(sprintf('Conv: %f15  %f15',s1,abs(s1/s6)))
disp(sprintf('   p: %f15  %f15',s1p,abs(s1p/s6)))
disp(sprintf('   n: %f15  %f15',s1n,abs(s1n/s6)))
disp(sprintf('Pstr: %f15  %f15',s6,abs(s6/s6)))
disp(sprintf('Sum:  %f15  %f15',s1+s2+s3+s4+s5+s6+s7,abs(s1+s2+s3+s4+s5+s6+s7)/s6))
disp('---------------------------------')
%--------------------------


figure(5)
title('uv budget')
hold on

b(:,29) = -b(:,29);

plot(b(:,1)/factl,b(:,29)/fact.*yy,'r')
plot(b(:,1)/factl,b(:,30)/fact.*yy,'k')
plot(b(:,1)/factl,b(:,31)/fact.*yy,'c')
plot(b(:,1)/factl,b(:,32)/fact.*yy,'g')
plot(b(:,1)/factl,b(:,33)/fact.*yy,'m')
plot(b(:,1)/factl,b(:,34)/fact.*yy,'m--')
plot(b(:,1)/factl,b(:,35)/fact.*yy,'b')

sum5(:,2) =+b(:,29)+b(:,30)+b(:,31)+b(:,32)+b(:,33)+b(:,34)+b(:,35);
plot(b(:,1)/factl,sum5(:,2)/fact.*yy,'y')



xlabel('\ity^+')
ylabel('loss   (inner scaling)   gain')

legend('convection','production','dissipation','turbdiff','pdiff','pstrain','visdiff','sum')



figure(6)
title('uu+vv+ww budget')
hold on

plot(b(:,1)/factl,.5*(b(:,8)+b(:,15)+b(:,22))/fact.*yy,'r')
plot(b(:,1)/factl,.5*(b(:,9)+b(:,16)+b(:,23))/fact.*yy,'k')
plot(b(:,1)/factl,.5*(b(:,10)+b(:,17)+b(:,24))/fact.*yy,'c')
plot(b(:,1)/factl,.5*(b(:,11)+b(:,18)+b(:,25))/fact.*yy,'g')
plot(b(:,1)/factl,.5*(b(:,12)+b(:,19)+b(:,26))/fact.*yy,'m')
plot(b(:,1)/factl,.5*(b(:,13)+b(:,20)+b(:,27))/fact.*yy,'m--')
plot(b(:,1)/factl,.5*(b(:,14)+b(:,21)+b(:,28))/fact.*yy,'b')




xlabel('\ity^+')
ylabel('loss   (inner scaling)   gain')

legend('convection','production','dissipation','turbdiff','pdiff','pstrain','visdiff')


%-------compute integrals-----------
a1=iDD*b(:,29);s1=a1(end)-a1(1);
a1=iDD*max(b(:,29),0);s1p=a1(end)-a1(1);
a1=iDD*min(b(:,29),0);s1n=a1(end)-a1(1);
a1=iDD*b(:,30);s2=a1(end)-a1(1);
a1=iDD*b(:,31);s3=a1(end)-a1(1);
a1=iDD*b(:,32);s4=a1(end)-a1(1);
a1=iDD*b(:,33);s5=a1(end)-a1(1);
a1=iDD*b(:,34);s6=a1(end)-a1(1);
a1=iDD*b(:,35);s7=a1(end)-a1(1);

disp('Reynolds stress uv')
disp('Relevant source terms; value ; value/prod')
disp(sprintf('Prod: %f15  %f15',s2,abs(s2/s2)))
disp(sprintf('Diss: %f15  %f15',s3,abs(s3/s2)))
disp(sprintf('Conv: %f15  %f15',s1,abs(s1/s2)))
disp(sprintf('   p: %f15  %f15',s1p,abs(s1p/s2)))
disp(sprintf('   n: %f15  %f15',s1n,abs(s1n/s2)))
disp(sprintf('Pstr: %f15  %f15',s6,abs(s6/s2)))
disp(sprintf('Sum:  %f15  %f15',s1+s2+s3+s4+s5+s6+s7,abs(s1+s2+s3+s4+s5+s6+s7)/s2))
disp('---------------------------------')
%--------------------------







%
% FILE OUTPUT
%


fid=fopen('bud-1000-k.out','wt');
fprintf(fid,'DNS of a turbulent zero-pressure gradient boundary layer\n');
fprintf(fid,'References:\n');
fprintf(fid,'Schlatter and Orlu, J. Fluid Mech., 659 (2010)\n');
fprintf(fid,'Schlatter et al., 2009, Bulletin APS, 54:19, page 59\n');
fprintf(fid,'Li and Schlatter, 2011, private communciation\n');
fprintf(fid,'\n');
fprintf(fid,'Integral quantities:\n');
fprintf(fid,'Re_{\\theta}   = %14.3f\n',res(3));
fprintf(fid,'Re_{\\delta^*} = %14.3f\n',res(2));
fprintf(fid,'Re_{\\tau}     = %14.4f\n',deltas(1)/lstar);
fprintf(fid,'H_{12}        = %14.6f\n',H12);
fprintf(fid,'c_f           = %14.9f\n',2*(utau/v(end,2))^2);
fprintf(fid,'\n');
fprintf(fid,'Wall-normal profiles:\n');
fprintf(fid,['y/\\delta_{99}       y+          conv+          prod+' ...
	     '          diss+          t-diff+        velp+' ...
	      '          vis-diff+      residual+\n']);


for i=1:ny
  fprintf(fid,'%13.7f  %13.7f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  \n',v(i,1)/deltas(1),v(i,1)/lstar,b(i,2)/fact,b(i,3)/fact,b(i,4)/fact,b(i,5)/fact,b(i,6)/fact,b(i,7)/fact,sum1(i,2)/fact);
end
fclose(fid);



fid=fopen('bud-1000-uu.out','wt');
fprintf(fid,'DNS of a turbulent zero-pressure gradient boundary layer\n');
fprintf(fid,'References:\n');
fprintf(fid,'Schlatter and Orlu, J. Fluid Mech., 659 (2010)\n');
fprintf(fid,'Schlatter et al., 2009, Bulletin APS, 54:19, page 59\n');
fprintf(fid,'Li and Schlatter, 2011, private communciation\n');
fprintf(fid,'\n');
fprintf(fid,'Integral quantities:\n');
fprintf(fid,'Re_{\\theta}   = %14.3f\n',res(3));
fprintf(fid,'Re_{\\delta^*} = %14.3f\n',res(2));
fprintf(fid,'Re_{\\tau}     = %14.4f\n',deltas(1)/lstar);
fprintf(fid,'H_{12}        = %14.6f\n',H12);
fprintf(fid,'c_f           = %14.9f\n',2*(utau/v(end,2))^2);
fprintf(fid,'\n');
fprintf(fid,'Wall-normal profiles:\n');
fprintf(fid,['y/\\delta_{99}       y+          conv+          prod+' ...
	     '          diss+          t-diff+        pdiff+' ...
	      '         pstrain+       vis-diff+      residual+\n']);


for i=1:ny
  fprintf(fid,'%13.7f  %13.7f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  \n',v(i,1)/deltas(1),v(i,1)/lstar,b(i,8)/fact,b(i,9)/fact,b(i,10)/fact,b(i,11)/fact,b(i,12)/fact,b(i,13)/fact,b(i,14)/fact,sum2(i,2)/fact);
end
fclose(fid);

fid=fopen('bud-1000-vv.out','wt');
fprintf(fid,'DNS of a turbulent zero-pressure gradient boundary layer\n');
fprintf(fid,'References:\n');
fprintf(fid,'Schlatter and Orlu, J. Fluid Mech., 659 (2010)\n');
fprintf(fid,'Schlatter et al., 2009, Bulletin APS, 54:19, page 59\n');
fprintf(fid,'Li and Schlatter, 2011, private communciation\n');
fprintf(fid,'\n');
fprintf(fid,'Integral quantities:\n');
fprintf(fid,'Re_{\\theta}   = %14.3f\n',res(3));
fprintf(fid,'Re_{\\delta^*} = %14.3f\n',res(2));
fprintf(fid,'Re_{\\tau}     = %14.4f\n',deltas(1)/lstar);
fprintf(fid,'H_{12}        = %14.6f\n',H12);
fprintf(fid,'c_f           = %14.9f\n',2*(utau/v(end,2))^2);
fprintf(fid,'\n');
fprintf(fid,'Wall-normal profiles:\n');
fprintf(fid,['y/\\delta_{99}       y+          conv+          prod+' ...
	     '          diss+          t-diff+        pdiff+' ...
	      '         pstrain+       vis-diff+      residual+\n']);


for i=1:ny
  fprintf(fid,'%13.7f  %13.7f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  \n',v(i,1)/deltas(1),v(i,1)/lstar,b(i,15)/fact,b(i,16)/fact,b(i,17)/fact,b(i,18)/fact,b(i,19)/fact,b(i,20)/fact,b(i,21)/fact,sum3(i,2)/fact);
end
fclose(fid);

fid=fopen('bud-1000-ww.out','wt');
fprintf(fid,'DNS of a turbulent zero-pressure gradient boundary layer\n');
fprintf(fid,'References:\n');
fprintf(fid,'Schlatter and Orlu, J. Fluid Mech., 659 (2010)\n');
fprintf(fid,'Schlatter et al., 2009, Bulletin APS, 54:19, page 59\n');
fprintf(fid,'Li and Schlatter, 2011, private communciation\n');
fprintf(fid,'\n');
fprintf(fid,'Integral quantities:\n');
fprintf(fid,'Re_{\\theta}   = %14.3f\n',res(3));
fprintf(fid,'Re_{\\delta^*} = %14.3f\n',res(2));
fprintf(fid,'Re_{\\tau}     = %14.4f\n',deltas(1)/lstar);
fprintf(fid,'H_{12}        = %14.6f\n',H12);
fprintf(fid,'c_f           = %14.9f\n',2*(utau/v(end,2))^2);
fprintf(fid,'\n');
fprintf(fid,'Wall-normal profiles:\n');
fprintf(fid,['y/\\delta_{99}       y+          conv+          prod+' ...
	     '          diss+          t-diff+        pdiff+' ...
	      '         pstrain+       vis-diff+      residual+\n']);


for i=1:ny
  fprintf(fid,'%13.7f  %13.7f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  \n',v(i,1)/deltas(1),v(i,1)/lstar,b(i,22)/fact,b(i,23)/fact,b(i,24)/fact,b(i,25)/fact,b(i,26)/fact,b(i,27)/fact,b(i,28)/fact,sum4(i,2)/fact);
end
fclose(fid);

fid=fopen('bud-1000-uv.out','wt');
fprintf(fid,'DNS of a turbulent zero-pressure gradient boundary layer\n');
fprintf(fid,'References:\n');
fprintf(fid,'Schlatter and Orlu, J. Fluid Mech., 659 (2010)\n');
fprintf(fid,'Schlatter et al., 2009, Bulletin APS, 54:19, page 59\n');
fprintf(fid,'Li and Schlatter, 2011, private communciation\n');
fprintf(fid,'\n');
fprintf(fid,'Integral quantities:\n');
fprintf(fid,'Re_{\\theta}   = %14.3f\n',res(3));
fprintf(fid,'Re_{\\delta^*} = %14.3f\n',res(2));
fprintf(fid,'Re_{\\tau}     = %14.4f\n',deltas(1)/lstar);
fprintf(fid,'H_{12}        = %14.6f\n',H12);
fprintf(fid,'c_f           = %14.9f\n',2*(utau/v(end,2))^2);
fprintf(fid,'\n');
fprintf(fid,'Wall-normal profiles:\n');
fprintf(fid,['y/\\delta_{99}       y+          conv+          prod+' ...
	     '          diss+          t-diff+        pdiff+' ...
	      '         pstrain+       vis-diff+      residual+\n']);


for i=1:ny
  fprintf(fid,'%13.7f  %13.7f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  \n',v(i,1)/deltas(1),v(i,1)/lstar,b(i,29)/fact,b(i,30)/fact,b(i,31)/fact,b(i,32)/fact,b(i,33)/fact,b(i,34)/fact,b(i,35)/fact,sum5(i,2)/fact);
end
fclose(fid);




break




