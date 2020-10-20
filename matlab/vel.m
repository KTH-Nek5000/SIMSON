%close all
%clear all
set(0,'DefaultLineLineWidth',1)


% 2500
xl=6000;zl=240;nx=8192;nz=768;
% 2500_fine
% xl=3000;zl=120;nx=4096;nz=480;
% 4000
% xl=3000;zl=120;nx=3072;nz=256;
% Skote
% xl=600;zl=34;nx=640;nz=128;

% high
%xl=750,zl=60,yl=50;nx=2048;nz=384;

re = 450;


p_ferrante = 0;
p_jimenez = 0;
p_spalart = 0;
p_ramis = 0;
p_osterlund = 0;
p_wu = 0;

cc='b';

% m fine
% r long
% b short


%u_fs=load('u_fs.dat');
%w_fs=load('w_fs.dat');
%beta = atan(w_fs(:,2)./u_fs(:,2));


%figure(1)
%hold on
%plot(u_fs(:,1),u_fs(:,2),'r')
%plot(w_fs(:,1),w_fs(:,2),'b')
%plot(u_fs(:,1),beta,'g')
%legend('u fs','w fs','beta')


%v=load('vel1130.data');
%v=load('vel2000.data');
%v=load('vel1130.data');
%v=load('vel0686.data');
%v=load('vel0364.data');
%v=load('vel0310.data');

%v=load('../2500_fine/vel0686.data');
%v=load('../4000/vel2000.data');
%v=load('../2500_fine/vel0500.data');


% 4000
%v=load('vel1500.data');
%v=load('vel5000.data');
%v=load('vel3500.data');
%v=load('vel4500.data');
%v=load('vel4000.data');
%v=load('vel1130.data');
%v=load('../2500/vel1500.data');cc='r';
%v=load('../4000/vel4000.data');cc='g';

%v=load('vel4000.data.21-99');cc='g';
%v=load('vel4000.data.100-199');cc='b';
%v=load('vel4000.data.100-481');cc='r';
%v=load('../4000/vel1000.data');cc='g';
%v=load('../2500_fine/vel1000.data');cc='k';


%v=load('vel1776.data');cc='b';
%v=load('vel1278.data');cc='b';
%v=load('vel0783.data');cc='b';
%v=load('vel1000.data');cc='b';

%v=load('../4000/vel4000.data.50');cc='k';




%
%COMPARISON OF TRANSITION LOCATION
%

% RE = 670
%v=load('vel0500.data');cc='k';  % probably OK (405)
%v=load('~/vel0405.data');cc='y';  % probably OK (405)
%v=load('vel0364.data');cc='k'; % PROBABLY JUST NOT OK
%v=load('~/scratch/TRIP/BL_strong3/stats/vel0364.data');cc='b'; % OK
%v=load('~/scratch/TRIP/BL_strong6/vel0512.data');cc='g';  % NOT OK
%v=load('~/scratch/TRIP/BL_strong7/vel1017.data');cc='r';  % NOT OK
%v=load('~/scratch/TRIP/BL_strong8/vel0870.data');cc='m';  % NOT OK
%v=load('~/scratch/TRIP/BL_strong9/vel0414.data');cc='y';  % OK
%v=load('~/scratch/TRIP/bl_low1/vel0898.data');cc='y';re=225; % OK

% RE = 1100
%v=load('vel0783.data');cc='k';
%v=load('~/scratch/TRIP/BL_strong3/stats/vel0790.data');cc='b';
%v=load('~/scratch/TRIP/BL_strong6/vel0942.data');cc='g';
%v=load('~/scratch/TRIP/BL_strong7/vel1489.data');cc='r';  % GOOD AGREEMENT
%v=load('~/scratch/TRIP/BL_strong8/vel1341.data');cc='m';
%v=load('~/scratch/TRIP/BL_strong9/vel0838.data');cc='y';
%v=load('~/scratch/TRIP/bl_low1/vel1777.data');cc='y';re=225; % OK

% RE = 1550
%v=load('vel1278.data');cc='k';
%v=load('~/scratch/TRIP/BL_strong3/stats/vel1278.data');cc='b';
%v=load('~/scratch/TRIP/BL_strong6/vel1438.data');cc='g';
%v=load('~/scratch/TRIP/BL_strong7/vel1981.data');cc='r';  % GOOD AGREEMENT
%v=load('~/scratch/TRIP/BL_strong8/vel1883.data');cc='m';
%v=load('~/scratch/TRIP/BL_strong9/vel1337.data');cc='y';
%v=load('~/scratch/TRIP/bl_low1/vel2756.data');cc='y';re=225; % OK
%v=load('~/scratch/TRIP/bl_high/vel0116.data');cc='k';re=1800; %  not OK

% RE = 2000
%v=load('vel1816.data');cc='k';
%v=load('~/scratch/TRIP/BL_strong3/stats/vel1824.data');cc='b';
%v=load('~/scratch/TRIP/BL_strong6/vel1970.data');cc='g';
%v=load('~/scratch/TRIP/BL_strong7/vel2500.data');cc='r'; % bit far
%v=load('~/scratch/TRIP/BL_strong8/vel2398.data');cc='m';
%v=load('~/scratch/TRIP/BL_strong9/vel1874.data');cc='y';
%v=load('~/scratch/TRIP/bl_low1/vel4000.data');cc='y';re=225; % a bit far
%v=load('~/scratch/TRIP/bl_high/vel0258.data');cc='k';re=1800; % OK

%------------------------------



%v=load('vel4500.data');cc='k';


%-----------------------------------------------

%v=load('vel3165.data');cc='k';
%v=load('vel3165.data.10');cc='b';
%v=load('vel3165.data.25');cc='g';  % OK
%v=load('vel3165.data.50');cc='r';


%v=load('vel1278.data');cc='k';
%v=load('vel1278.data.10');cc='b';   % OK
%v=load('vel1278.data.25');cc='g';
%v=load('vel1278.data.50');cc='r';

%v=load('vel4639.data.25');cc='k';
%v=load('vel4639.data.10');cc='b';   
%v=load('vel4639.data.25');cc='g';   % OK
%v=load('vel4639.data.50');cc='r';

%v=load('vel1000.data.10');cc='r';


% Re_tau = 550: 1278
%          360: 0686
%         1000: 3165
%v=load('vel0686.data');cc='r';
%v=load('vel1278.data.10');cc='r';
%v=load('vel3165.data.10');cc='r';
%-------------------
% DATA FOR WEBSITE
%v=load('vel0364.data');cc='r';
%v=load('vel0686.data');cc='r';
%v=load('vel1130.data.10');cc='r';
%v=load('vel1816.data.10');cc='r';
%v=load('vel2500.data.10');cc='r';
%v=load('vel3165.data.25');cc='r';
%v=load('vel3500.data.25');cc='r';
%v=load('vel4000.data.25');cc='r';
%v=load('vel4500.data.25');cc='r';
%v=load('vel4639.data.25');cc='r';
%-------------------

v=load('vel0171.data');cc='b';


% Ramis      2532,       3640, 4080
% Osterlund: 2532, 3058, 3651,      4613

% File   Re_theta   Ramis   Osterlund
% 5000   4270
% 4639   4061         B5  
% 4500   3971
% 4000   3628         B4         X
% 3500   3273
% 3165   3032                   X
% 3000   2911
% 2500   2536         B2         X
% 2000   2148
% 1816   2001 *
% 1776   1968
% 1500   1740
% 1278   1551 *
% 1130   1420
% 1000   1303
%  783   1100 *
%  686   1000 
%  500    819
%  364    677
 
% structure: 
% u v w urms vrms wrms uv uw vw omx omy omz omxrms
% omyrms omzrms p prms .... 
% 25-26  S(u) F(u)
ny=length(v);
yl=v(end,1);



[re_theta_u,re_deltastar_u,utau,lstar_u,H12_u] = comp_re(yl,v(:,2),re);
[res,deltas,utau,lstar,H12] = comp_re1(yl,v(:,2),re);


disp(sprintf('domain: %f x %f x %f',xl,yl,zl));
disp(sprintf('resolution: %i x %i x %i',nx,ny,nz));
disp(sprintf('re_theta = %f',re_theta_u));
disp(sprintf('re_dstar = %f',re_deltastar_u));
disp(sprintf('re_tau   = %f',res(1)*utau));
disp(sprintf('tau_rms  = %f',v(1,16)/utau*lstar));
disp(sprintf('H12      = %f',H12));
disp(sprintf('H12_99   = %f',deltas(5)/deltas(6)));


disp(sprintf('    deltax+=%f',xl/nx/lstar_u));
disp(sprintf('    deltaz+=%f',zl/nz/lstar_u));
disp(sprintf('    deltay+=%f',v(2,1)/lstar_u));
disp(sprintf('    deltay+=%f',(v((ny+1)/2,1)-v((ny-1)/2,1))/lstar_u));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                             %%
%% Mean Velocity Profile                                       %%
%%                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on

plot(v(:,1)/lstar_u,v(:,2)/utau,cc,'linewidth',2)

yp=0.1:0.1:15;
if (find(strcmp(get(get(gca,'children'),'displayname'),'log1')))
else
  plot(yp,yp,'g','displayname','log1')
end
yp=[5 2000];
if (find(strcmp(get(get(gca,'children'),'displayname'),'log2')))
else
  plot(yp,1/0.41*log(yp)+5.2,'g','displayname','log2')
end





set(gca,'xscale','log')

load ~/DATA/SPALART/spalart.mat;

if p_spalart==1
  plot(spalart300.u(:,3),spalart300.u(:,4),'m--')
  plot(spalart670.u(:,3),spalart670.u(:,4),'m--')
  plot(spalart1410.u(:,3),spalart1410.u(:,4),'m--')
end

jim = load('~/DATA/chan2003.dat');
if (p_jimenez==1)
  if (find(strcmp(get(get(gca,'children'),'displayname'),'jim_chan2003')))
  else
    plot(jim(:,2),jim(:,3),'b--','displayname','jim_chan2003');
  end
end




load('~/DATA/wu_moin_2009.mat');
if (p_wu==1)
  plot(wu_u(:,1),wu_u(:,2),'k--');
end



axis([1 3e3 0 27])
xlabel('\ity^+')
ylabel('\itU^+')




% Osterlund data

load('~/DATA/OSTERLUND/zpg_v1.mat');

if (p_osterlund==1)
  for i=1:3
    plot(sw1(i).yp,sw1(i).up,'k--')
  end

  %for i=11:12
  %  plot(sw1(i).yp,sw1(i).up,'k--')
  %end
end  

% Ramis
load /home/x_phisc/DATA/Ramis/20110714/ZPGTBL_JFM.mat


R70=load('~/DATA/Ramis/20090615/finalG_Exp70.mat');
R78=load('~/DATA/Ramis/20090615/finalG_Exp78.mat');
R79=load('~/DATA/Ramis/20090615/finalG_Exp79.mat');



if (p_ramis==1)
  plot(R70.H50(:,6),R70.H50(:,7),'k.-')
  plot(R78.H50(6:end,6),R78.H50(6:end,7),'k.-')
  plot(R79.H50(3:end,6),R79.H50(3:end,7),'k.-')

  plot(ZPGTBL_T.B5.data(:,6),ZPGTBL_T.B5.data(:,7),'ko-')
  plot(ZPGTBL_T.B2.data(:,6),ZPGTBL_T.B2.data(:,7),'ko-')
  plot(ZPGTBL_T.B4.data(:,6),ZPGTBL_T.B4.data(:,7),'ko-')


end




ferrante=load('~/DATA/ferrante.dat');
if (p_ferrante==1)
  plot(ferrante(:,1),ferrante(:,2),'m');
end

jim1100=load('~/DATA/JIMENEZ/TBL2000/profiles/jim.1100.prof');
jim1551=load('~/DATA/JIMENEZ/TBL2000/profiles/jim.1551.prof');
jim1968=load('~/DATA/JIMENEZ/TBL2000/profiles/jim.1968.prof');
if (p_jimenez==1)
  plot(jim1100(:,2),jim1100(:,7),'r-.');
  plot(jim1551(:,2),jim1551(:,7),'r-.');
  plot(jim1968(:,2),jim1968(:,7),'r-.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                             %%
%% Velocity Defect Profile                                     %%
%%                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(50)
hold on

% Rotta-Clauser thickness
%lfact = lstar*res(2)
lfact = deltas(1);

plot(v(:,1)/lfact,v(end,2)/utau-v(:,2)/utau,cc,'linewidth',2)


yd=[0.01 1];
if (find(strcmp(get(get(gca,'children'),'displayname'),'log2')))
else
  plot(yd,-1/0.41*log(yd)+3,'g','displayname','log2')
end


if (p_jimenez==1)
  i=find(jim1100(:,7)/jim1100(end,7)>0.99);d99=jim1100(i(1),2);
  plot(jim1100(:,2)/d99,jim1100(end,7)-jim1100(:,7),'r-.');
  i=find(jim1551(:,7)/jim1551(end,7)>0.99);d99=jim1551(i(1),2);
  plot(jim1551(:,2)/d99,jim1551(end,7)-jim1551(:,7),'r-.');
  i=find(jim1968(:,7)/jim1968(end,7)>0.99);d99=jim1968(i(1),2);
  plot(jim1968(:,2)/d99,jim1968(end,7)-jim1968(:,7),'r-.');
end



axis([0.001 2 0 20])
xlabel('\ity/\delta_{99}')
ylabel('\itU^+')

set(gca,'xscale','log')
%break


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                             %%
%% vorticity fluctuations                                      %%
%%                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(250)
hold on
title('Vorticity fluctuations - inner scaling')

ff = utau/lstar;
plot(v(:,1)/lstar,v(:,14)/(ff),cc,'linewidth',1)
plot(v(:,1)/lstar,v(:,15)/(ff),cc,'linewidth',1)
plot(v(:,1)/lstar,v(:,16)/(ff),cc,'linewidth',1)

if (p_jimenez==1)
  plot(jim1100(:,2),jim1100(:,20),'r-.');
  plot(jim1100(:,2),jim1100(:,21),'r-.');
  plot(jim1100(:,2),jim1100(:,22),'r-.');
end

figure(251)
hold on
title('Vorticity fluctuations - mixed scaling')

ff = sqrt(utau/lstar * 1/deltas(2));
plot(v(:,1)/deltas(1),v(:,14)/(ff),cc,'linewidth',1)
plot(v(:,1)/deltas(1),v(:,15)/(ff),cc,'linewidth',1)
plot(v(:,1)/deltas(1),v(:,16)/(ff),cc,'linewidth',1)

figure(252)
hold on
title('Vorticity fluctuations - outer scaling')

ff = 1/deltas(2);
plot(v(:,1)/deltas(1),v(:,14)/(ff),cc,'linewidth',1)
plot(v(:,1)/deltas(1),v(:,15)/(ff),cc,'linewidth',1)
plot(v(:,1)/deltas(1),v(:,16)/(ff),cc,'linewidth',1)


%break




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                 %%
%% Indicator Quantity y+*du+/dy+                                   %%
%%                                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
hold on

[y,D] = chebdif(ny,1);
D=D*(-2/yl);

du=D*v(:,2);


%figure(100)
%plot(v(:,1)./lstar_u, du./(utau/lstar_u))
%break


%figure(2000)
%hold on
%plot(v(:,1)./lstar_u , 1./(v(:,1)./lstar_u.*du./(utau/lstar_u)),cc,'linewidth',2)


plot(v(:,1)./lstar_u , v(:,1)./lstar_u.*du./(utau/lstar_u),cc,'linewidth',2)




%plot(v(:,1)./lstar_u , (v(:,1)./lstar_u+5).*du./(utau/lstar_u),'g','linewidth',2)


set(gca,'xscale','log')
plot([10 1000],[1/0.384 1/0.384],'g--')
plot([10 1000],[1/0.4 1/0.4],'g')
plot([10 1000],[1/0.42 1/0.42],'g--')

% Nagib profile
yp=logspace(-1,3,100);

kappa = 0.384;
b0=1e-2*kappa^-1;
b1=1.1e-2;
b2=1.1e-4;
P23 = b0*(1+b1*yp+b2*yp.^2)./(1+b1*yp+b2*yp.^2+kappa*b0*b2*yp.^3);
h1=-1e-2;
h2=6e-3;
h3=9.977e-4;
h4=2.2e-5;
h5=1e-6;
P45 = (1-b0)*(1+h1*yp+h2*yp.^2)./(1+h1*yp+h2*yp.^2+h3*yp.^3+h4*yp.^4+h5*yp.^5);

plot(yp,yp.*(P23+P45),'k');



if (p_jimenez==1)
  % channel
  plot(jim(:,2),jim(:,2).*jim(:,7),'b--');
  
  % Boundary layer
   plot(jim1100(:,2),jim1100(:,2).*jim1100(:,16),'r-.')
   plot(jim1551(:,2),jim1551(:,2).*jim1551(:,16),'r-.')
   plot(jim1968(:,2),jim1968(:,2).*jim1968(:,16),'r-.')
end



axis([1 3000 0 6])
xlabel('\ity^+')
ylabel('\Xi')

%abort('333')

% du+/dy+ as fct of y+
%
figure(20)
hold on
plot(v(:,1)./lstar_u ,1-du./(utau/lstar_u),cc,'linewidth',2)


%break


%
% Total Shear Stress
%

if (1==1)

  figure(11)
  hold on
  lfact = v(ny,2)/utau*re_deltastar_u/re;
  dd=-v(:,8)./utau^2+du./(utau/lstar_u);
  ddd=D*dd*lfact;
  plot(v(:,1)/lfact,ddd,cc)
  set(gca,'xscale','log')
  axis([1e-4 4e-1 -6 0])
  title('total shear stress derivative')


  figure(12)
  title('total shear stress')
  hold on
  plot(v(:,1)/lstar_u,-v(:,8)./utau^2,cc)
  plot(v(:,1)/lstar_u,+du./(utau/lstar_u),cc)
  plot(v(:,1)/lstar_u,dd,cc)


  figure(13)
  title('1-dudy and 1-tau')
  hold on
  plot(v(:,1)/lstar_u,1-du./(utau/lstar_u),cc)
  plot(v(:,1)/lstar_u,1-dd,'r--')
  grid on

%  break
  
  %figure(13)
  %hold on
  %du=D*D*v(:,2);
  %plot(v(:,1)/lstar_u,du)

end

%
% Reynolds stresses
%

v(:,5) = v(:,5).^2;
v(:,6) = v(:,6).^2;
v(:,7) = v(:,7).^2;

%
% Reynolds stresses inner scaling
%

figure(4)
hold on
plot(v(:,1)/lstar_u,sqrt(v(:,5)/utau.^2),cc)
plot(v(:,1)/lstar_u,sqrt(v(:,6)/utau.^2),cc)
plot(v(:,1)/lstar_u,sqrt(v(:,7)/utau.^2),cc)

plot(v(:,1)/lstar_u,v(:,8)/utau.^2,cc)
%plot(v(:,1)/lstar_u,v(:,9)/utau.^2,'r--')
%plot(v(:,1)/lstar_u,v(:,10)/utau.^2,'g--')


if (p_spalart==1)
plot(spalart1410.u(:,3), sqrt(spalart1410.u(:,5)),'m--')
plot(spalart1410.u(:,3), sqrt(spalart1410.u(:,6)),'m--')
plot(spalart1410.u(:,3), sqrt(spalart1410.u(:,7)),'m--')
plot(spalart1410.u(:,3),(spalart1410.u(:,8)),'m--')

plot(spalart670.u(:,3), sqrt(spalart670.u(:,5)),'m--')
plot(spalart670.u(:,3), sqrt(spalart670.u(:,6)),'m--')
plot(spalart670.u(:,3), sqrt(spalart670.u(:,7)),'m--')
plot(spalart670.u(:,3), (spalart670.u(:,8)),'m--')

plot(spalart300.u(:,3), sqrt(spalart300.u(:,5)),'m--')
plot(spalart300.u(:,3), sqrt(spalart300.u(:,6)),'m--')
plot(spalart300.u(:,3), sqrt(spalart300.u(:,7)),'m--')
plot(spalart300.u(:,3), (spalart300.u(:,8)),'m--')
end

if (p_ferrante==1)
  plot(ferrante(:,1),sqrt(ferrante(:,3)),'m');
  plot(ferrante(:,1),sqrt(ferrante(:,4)),'m');
  plot(ferrante(:,1),(ferrante(:,5)),'m');
end

if (p_jimenez==1)
  plot(jim1100(:,2),jim1100(:,3),'r-.')
  plot(jim1100(:,2),jim1100(:,4),'r-.')
  plot(jim1100(:,2),jim1100(:,5),'r-.')
  plot(jim1100(:,2),jim1100(:,6),'r-.')
  
  plot(jim1551(:,2),jim1551(:,3),'r-.')
  plot(jim1551(:,2),jim1551(:,4),'r-.')
  plot(jim1551(:,2),jim1551(:,5),'r-.')
  plot(jim1551(:,2),jim1551(:,6),'r-.')
  
  plot(jim1968(:,2),jim1968(:,3),'r-.')
  plot(jim1968(:,2),jim1968(:,4),'r-.')
  plot(jim1968(:,2),jim1968(:,5),'r-.')
  plot(jim1968(:,2),jim1968(:,6),'r-.')
end


if (p_wu==1)
  plot(wu_urms(:,1),wu_urms(:,2),'k--');
  plot(wu_vrms(:,1),wu_urms(:,2),'k--');
  plot(wu_wrms(:,1),wu_urms(:,2),'k--');
  plot(wu_uv(:,1),wu_uv(:,2),'k--');
end

if (p_ramis==1)
  plot(R70.H50(:,6),R70.H50(:,8),'k.-')
  plot(R78.H50(6:end,6),R78.H50(6:end,8),'k.-')
  plot(R79.H50(3:end,6),R79.H50(3:end,8),'k.-')

  plot(ZPGTBL_T.B5.data(:,6),ZPGTBL_T.B5.data(:,8),'ko-')
  plot(ZPGTBL_T.B2.data(:,6),ZPGTBL_T.B2.data(:,8),'ko-')
  plot(ZPGTBL_T.B4.data(:,6),ZPGTBL_T.B4.data(:,8),'ko-')
end

%marusic=load('marusic_1107.data');
%plot(marusic(:,1),sqrt(marusic(:,2)),'g')


xlabel('\ity^+')
ylabel('\itu_{i,\rmrms}^+,  (u\primev\prime)^+')

axis([0 700 -1 3])

%
% Reynolds stresses outer scaling
%


figure(5)
hold on
lfact = re_deltastar_u/re;
plot(v(:,1)/lfact,(v(:,5)),cc)
plot(v(:,1)/lfact,(v(:,6)),cc)
plot(v(:,1)/lfact,(v(:,7)),cc)
plot(v(:,1)/lfact,(v(:,8)),cc)

xlabel('\ity/\delta^*')
ylabel('\itu_i\primeu_j\prime')

axis([0 7 -3e-3 17e-3])


if (p_spalart==1)
  plot(spalart1410.u(:,3)/90, (spalart1410.u(:,5))/spalart1410.u(end,4).^2,'m--')
  plot(spalart1410.u(:,3)/90, (spalart1410.u(:,6))/spalart1410.u(end,4).^2,'m--')
  plot(spalart1410.u(:,3)/90, (spalart1410.u(:,7))/spalart1410.u(end,4).^2,'m--')
  plot(spalart1410.u(:,3)/90, (spalart1410.u(:,8))/spalart1410.u(end,4).^2,'m--')

  plot(spalart670.u(:,3)/52, (spalart670.u(:,5))/spalart670.u(end,4).^2,'m--')
  plot(spalart670.u(:,3)/52, (spalart670.u(:,6))/spalart670.u(end,4).^2,'m--')
  plot(spalart670.u(:,3)/52, (spalart670.u(:,7))/spalart670.u(end,4).^2,'m--')
  plot(spalart670.u(:,3)/52, (spalart670.u(:,8))/spalart670.u(end,4).^2,'m--')
end

%
% Wind tunnel quantity
%

figure(6)
hold on
plot(v(:,2),sqrt(v(:,5)),cc)

if (p_spalart==1)
  plot(spalart1410.u(:,4)./spalart1410.u(end,4),sqrt(spalart1410.u(:,5))./spalart1410.u(end,4),'m--')
  plot(spalart670.u(:,4)./spalart670.u(end,4),sqrt(spalart670.u(:, 5))./spalart670.u(end,4),'m--')
end
  
  
xlabel('\itU')
ylabel('\itu_{\rmrms}')
axis([0 1 0 0.14])


if (p_ramis==1)
  plot(R70.H50(:,2)/R70.H50(end,2),R70.H50(:,3)/R70.H50(end,2),'k.-')
  plot(R78.H50(6:end,2)/R78.H50(end,2),R78.H50(6:end,3)/R78.H50(end,2),'k.-')
  plot(R79.H50(3:end,2)/R79.H50(end,2),R79.H50(3:end,3)/R79.H50(end,2),'k.-')

  plot(ZPGTBL_T.B5.data(:,7)/ZPGTBL_T.B5.data(end,7),ZPGTBL_T.B5.data(:,8)/ZPGTBL_T.B5.data(end,7),'ro-')
  plot(ZPGTBL_T.B2.data(:,7)/ZPGTBL_T.B2.data(end,7),ZPGTBL_T.B2.data(:,8)/ZPGTBL_T.B2.data(end,7),'ro-')
  plot(ZPGTBL_T.B4.data(:,7)/ZPGTBL_T.B4.data(end,7),ZPGTBL_T.B4.data(:,8)/ZPGTBL_T.B4.data(end,7),'ro-')
end

if (p_ferrante==1)
  ff = ferrante(end,2);
  plot(ferrante(:,2)/ff,sqrt(ferrante(:,3))/ff,'m')
end
  
if (p_jimenez==1)
  ff=jim1100(end,7);
  plot(jim1100(:,7)/ff,jim1100(:,3)/ff,'r-.')
  ff=jim1551(end,7);
  plot(jim1551(:,7)/ff,jim1551(:,3)/ff,'r-.')
  ff=jim1968(end,7);
  plot(jim1968(:,7)/ff,jim1968(:,3)/ff,'r-.')
end
  
%
% Pressure Fluctuations
%

figure(7)
subplot(2,1,1)
hold on
lfact = lstar;
pfact = utau;

plot(v(:,1)/lfact,(v(:,18)/pfact).^1,cc)
xlabel('\ity^+')
ylabel('\itp_{\rmrms}/\rho u_\tau U_0')

subplot(2,1,2)
hold on
lfact = deltas(1);
pfact = utau.^2;
plot(v(:,1)/lfact,(v(:,18)/pfact).^1,cc)
xlabel('\ity/\delta_{99}')
ylabel('\it(p_{\rmrms}^{+})^2')
set(gca,'xscale','log')


% STREAMWISE MOMENTUM EQUATION

figure(51)
hold on

t1 = (v(:,2)/utau).*(v(:,21)/(utau/lstar));
t2 = (v(:,3)/utau).* (D*v(:,2)/(utau/lstar));
t3 = -D*D*v(:,2)/(utau/lstar/lstar);
t4 =  D*v(:,8)/(utau*utau/lstar);
t5 =  v(:,23)/(utau*utau/lstar);   % pressure
t6 = -v(:,22)/(utau/lstar/lstar);  %
t7 = v(:,24)/(utau*utau/lstar);    %

%t5=t5*0;
%t6=t6*0;

plot(v(:,1)/lstar, t1,'b' )
plot(v(:,1)/lstar, t2,'r' )
plot(v(:,1)/lstar, t3,'g' )
plot(v(:,1)/lstar, t4,'m' )
plot(v(:,1)/lstar, t5,'k' )
plot(v(:,1)/lstar, t6,'c' )
plot(v(:,1)/lstar, t7,'r--' )

plot(v(:,1)/lstar, t1+t2+t3+t4+t5+t6+t7,'y')
%plot(v(:,1)/lstar, -t6,'r')
%break

figure(52)
hold on 

DD=D;
DD(1,:)=0;
DD(1,1)=1;


%t8 = -(t1+t2+t3+t5+t6+t7);
%t9 = inv(DD)*t8;
%plot(v(:,1)/lstar,t9/lstar,'r');
%plot(v(:,1)/lstar,v(:,8)/utau/utau,'b');


plot(v(:,1)/lstar,t1+t2,cc);
%plot(v(:,1)/lstar,t1,'b--');
set(gca,'xscale','log');


% Skewness and flatness



v(find(v(:,1)>2*deltas(1),1):end,25) = 0;
v(find(v(:,1)>2*deltas(1),1):end,26) = 3;

figure(300)
hold on
plot(v(:,1)/lstar,v(:,25),cc)
xlabel('\ity^+');
ylabel('S(u)');
ylim([-3 3 ])
figure(301)
hold on
plot(v(:,1)/lstar,v(:,26),cc)
xlabel('\ity^+');
ylabel('F(u)');
ylim([0 50 ])


% v and w
figure(302)
hold on
plot(v(:,1)/lstar,v(:,3)/utau,cc)
xlabel('\ity^+');
ylabel('v');




%
% FILE OUTPUT
%
[y,D] = chebdif(ny,1);
D=D*(-2/yl);
du=D*v(:,2);

fid=fopen('vel.out','wt');
fprintf(fid,'%% DNS of a turbulent zero-pressure gradient boundary layer\n');
%fprintf(fid,'Reference: Schlatter et al., 2010, private communication\n ');
%fprintf(fid,'Reference: Schlatter et al., Phys. Fluids 21, 051702 (2009)\n ');
fprintf(fid,'%% References:\n');
fprintf(fid,'%% Schlatter and Orlu, 2010, J. Fluid Mech., 659 (2010)\n');
fprintf(fid,'%% Schlatter et al., 2009, Bulletin APS, 54:19, page 59\n');

fprintf(fid,'%%\n');
fprintf(fid,'%% Integral quantities:\n');
fprintf(fid,'%% Re_{\\theta}   = %14.3f\n',res(3));
fprintf(fid,'%% Re_{\\delta^*} = %14.3f\n',res(2));
fprintf(fid,'%% Re_{\\tau}     = %14.4f\n',deltas(1)/lstar);
fprintf(fid,'%% H_{12}        = %14.6f\n',H12);
fprintf(fid,'%% c_f           = %14.9f\n',2*(utau/v(end,2))^2);
fprintf(fid,'%%\n');
fprintf(fid,'%% Wall-normal profiles:\n');
fprintf(fid,'%% y/\\delta_{99}       y+          U+          urms+       vrms+       wrms+       uv+         p+            prms+       pu+         pv+         S(u)        F(u)        dU+/dy+     V+        omxrms^+    omyrms^+    omzrms^+ \n');


for i=1:ny
  fprintf(fid,'%13.7f  %13.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f \n',v(i,1)/deltas(1),v(i,1)/lstar,v(i,2)/utau,sqrt(v(i,5))/utau,sqrt(v(i,6))/utau,sqrt(v(i,7))/utau,v(i,8)/utau.^2,(v(i,17)-v(1,17))/utau.^2,v(i,18)/utau^2,v(i,19)/utau^3,v(i,20)/utau^3,v(i,25),v(i,26),du(i)/utau*lstar,v(i,3)/utau,v(i,14)/utau*lstar,v(i,15)/utau*lstar,v(i,16)/utau*lstar);
end
fclose(fid);
break




