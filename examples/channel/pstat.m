clear all
close all
system("./read_stat.sh xy.stat stat_vel.in  stat_vel.data");

stat=load('stat_vel.data');
y=stat(:,1);
ny=size(y,1);
[yy D]=chebdif(ny,1);
re=4200;

dudy = D*stat(:,2);
utau = sqrt((dudy(end)-dudy(1))/2/re);
lstar = 1/utau/re;


figure(1)
plot(y,stat(:,2))
xlabel('y/h')
ylabel('u/u_{cl}')

figure(2)
hold on
yp=linspace(0,10,20);
plot(yp,yp,'g')

yp = linspace(5,300,20);
plot(yp,1/0.41*log(yp)+5.2,'g')

yy = y(1:(end+1)/2)+1;
uu = ( stat(1:(ny+1)/2,2) + stat(ny:-1:(ny+1)/2,2) )/2;
semilogx(yy/lstar,uu/utau);
xlabel('y^+')
ylabel('U^+')



  
