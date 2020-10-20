clear all
close all
filename='t2000.u';
[vel,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv]= readdns(filename);

% the variables are:
% vel : velocity in Fourier space, contains all velocity components
% xF  : coordinates in x (note, there is one more point at the end)
% yF  : coordinates in y (note, top is index 1 and wall last index)
% zF  : coordinates in z (note, there is one more point at the end)
%
% t   : time
% Re  : Reynolds number
% rest not important

[phys,NNx,NNy,NNz]=fou2phys(vel,0,0);

% now phys is the physical-space velocity
% NNx : grid points in x (one point shorter than xF)
% ...

u=phys(:,:,1+0*NNy:1*NNy);
v=phys(:,:,1+1*NNy:2*NNy);
w=phys(:,:,1+2*NNy:3*NNy);

% now we have u,v,w as individual arrays

surf(xF(1:end-1),zF(1:end-1),u(:,:,60)')
title(sprintf('This is the x/z plane at y=%f',yF(60)))

% just to show an application...