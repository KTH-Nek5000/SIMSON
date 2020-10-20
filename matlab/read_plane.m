clear all
close all

% read binary file
fid=fopen('~/plane_xz_0.dat','r','ieee-le');
eor = fread(fid,1,'int32');
nx = fread(fid,1,'int32');
nz = fread(fid,1,'int32');
eor = fread(fid,2,'int32');
x = fread(fid,nx,'float64');
eor = fread(fid,2,'int32');
z = fread(fid,nz,'float64');
eor = fread(fid,2,'int32');
p = fread(fid,nx*nz,'float64');
p=reshape(p,nx,nz);
fclose(fid);

% plot plane
pcolor(x,z,p')
shading interp
axis tight
