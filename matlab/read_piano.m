
clear all
close all

nx = 768+1;
nz = 1024+1;

p = load('/home/guest/pschlatt/pianoxz.dat');
x = reshape(p(:,1),nx,nz);
z = reshape(p(:,2),nx,nz);
f = reshape(p(:,3),nx,nz);

pcolor(x,z,f);
shading flat

