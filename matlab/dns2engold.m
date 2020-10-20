% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
% Program to output DNS data on EnSight Gold format
% Arguments as follows
% Flag Filename [Mean flow file] Casename [start stop step]
%
% Flag:
%  -d to subtract mean flow file
%  -u for more than one file
%  -px for zero padding in streamwise Fourier space
%  -pz for zero padding in spanwise Fourier space
%

%
% Clear all variabels
%
clear all;

%
% Default Settings
%
nsteps=1;
nstart=1;
ninc=1;
casename='test';
titlestring='Test';
ex=0;                  % Extra number of points in x-direction in
                       % physical space
ez=0;                  % Extra number of points in z-direction in
                       % physcial space
l2=true;               % Compute lambda2
l2d=0.;                % Initiate lambda2 data
grad=false;             % Compute velocity gradient
gradp=0.;
vort=false;            % Compute vorticity
vortf=0.;
%
% Read mean flow and other field
%
[mf,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv]=readdns('bl6000_mf_streak_new_fringe_t0.u');
[vel,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv]=readdns('ebl6000.u');

%[mf,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv]=readdns('ch001-av000.u');
%[vel,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv]=readdns('ch001-av275.u');%
%
% Make Fourier transform
%
%[cnf,NNx,NNy,NNz]=fou2phys(vel-mf,ex,ez);
[cnf,NNx,NNy,NNz]=fou2phys(vel,ex,ez);

%
% Compute wavenumber vectors
%
kxvec=linspace(0,2*pi/Lx*(NNx/2-1),NNx/2);
kzvec=linspace(0,2*pi/Lz*(NNz/2-1),NNz/2);
kzvec=[kzvec -fliplr(kzvec(2:end))];

%
% Compute physical grid
%
x=linspace(0,Lx,NNx);
z=linspace(0,Lz,NNz);

xd=[x(2:end)-x(1:end-1) 0];
yFd=[yF(1:end-1)'-yF(2:end)' 0];
zd=[z(2:end)-z(1:end-1) 0];

[X,Y,Z]=meshgrid(x,yF,z);

%
% Compute velocity gradient tensor
%
if grad
  [dxfou,dyfou,dzfou]=gradfield(vel,yF,kxvec,kzvec);
  dxfoup=fou2phys(dxfou,ex,ez);
  dyfoup=fou2phys(dyfou,ex,ez);
  dzfoup=fou2phys(dzfou,ex,ez);
  gradp=cat(3,dxfoup,dyfoup,dzfoup);
end
%
% Compute vorticity
%
if vort
%  [vortf] = curlfield(X,Y,Z,zzu,zzv,zzw);
end

%
% Compute mean flow
%
%  for j=1:nyp
%    for i=1:nx
%      umeanr = 0.
%      for k=1:nz
%	umeanr = umeanr + vel(i,j,k,1)
%      end
%      umeanr = umeanr/nz
%      for k=1:nz
%	ar(i,j,k,1) = vel(i,j,k,1) - umeanr
%      end
%    end
%  end

%
% Compute lambda2
%
if l2
  [l2d,dxfou,dyfou,dzfou]=lambda2(vel,yF,kxvec,kzvec,ex,ez);
end

%
% Write case file
%
status = writecase(casename,vel,NNx,NNy,NNz,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv,nsteps,nstart,ninc,vort,grad,l2);
status = writeengold(casename,titlestring,NNx,NNy,NNz,xd,yFd,zd,Lx,Ly,Lz,cnf,X,Y,Z,vort,vortf,grad,gradp,l2,l2d);
