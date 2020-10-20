clear all
close all


name = '/scratch/pschlatt/pod/plane.dat';

npl = 1250;

fid=fopen(name,'r','ieee-be.l64');

% read header
eor=fread(fid,1,'int32');
re=fread(fid,1,'float64');
log=fread(fid,1,'int32');
xl=fread(fid,1,'float64');
zl=fread(fid,1,'float64');
t=fread(fid,1,'float64');
z0=fread(fid,1,'float64');
eor=fread(fid,2,'int32');

nx=fread(fid,1,'int32'); 
ny=fread(fid,1,'int32'); 
nz=fread(fid,1,'int32'); 
nsymm=fread(fid,1,'int32'); 
eol=fread(fid,2,'int32');  

tpl=fread(fid,1,'int32');
ivar=fread(fid,1,'int32');
cpl=fread(fid,1,'float64');
fltype=fread(fid,1,'int32');
dstar=fread(fid,1,'float64');
eol=fread(fid,1,'int32');


% do the scaling
xl=xl/dstar;
yl=2/dstar;
zl=zl/dstar;

%grid
x=xl/nx*[0:nx-1];
for i=1:ny
%    y(i)=(cos(pi*(i-1)/(ny-1))+1)/2*yl;
    y(i)=(cos(pi*(i-1)/(ny-1)));
end
z=zl/nz*[0:nz]-zl/2;
[zz,xx] = meshgrid(z,x);



uxy=zeros(nx,ny,npl);

% read planes
for i=1:npl

   
  eor=fread(fid,1,'int32');
  t(i)=fread(fid,1,'float64')/dstar;
  disp(sprintf('Plane %0.5i t=%4.2f',i,t(i)));
  
  eor=fread(fid,1,'float64');
  eor=fread(fid,2,'int32');
  d=fread(fid,[nx ny],'float64');


  uxy(:,:,i) = d;
  eor=fread(fid,1,'int32');

  % plot plane
  if (1==0) 
    figure(10)
    clf
    hold on
    pcolor(x,y,uxy(:,:,i)')
    contour(x,y,uxy(:,:,i)',[0.1:0.1:1],'k');
    hold off
    axis([0 xl -1 1])
    caxis([0 1])
    shading interp
    title(sprintf('t=%4.1f',t(i)))
    drawnow
  end
  
end

fclose(fid);

