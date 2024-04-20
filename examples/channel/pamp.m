clear all
close all

amp=load('amp.stat');
re=4200;
reb = amp(1,13);
t=amp(:,1);
  
figure (1)
utau1 = sqrt(amp(:,14)/re);
utau2 = sqrt(-amp(:,15)/re);
retau1 = utau1*re;
retau2 = utau2*re;
retau = sqrt((mean(amp(floor(end/2):end,14))-mean(amp(floor(end/2):end,15)))/2/re)*re;


hold on
plot(t,retau1,'r');
plot(t,retau2,'b'); 
plot([t(floor(end/2)) t(end)],[retau retau],'g')
legend('upper wall','lower wall')
xlabel('time')
ylabel('Re_\tau')
title(sprintf('Channel flow with Re_b=%5.0f, Re_t=%5.1f',reb,retau))


  


