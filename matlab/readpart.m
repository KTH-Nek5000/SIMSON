%close all
clear all


for i=1:2000
  name = 'part-00069.stat';
  
  name=sprintf('part-%05i.stat',i);
  
  fid=fopen(name,'r','ieee-be.l64');
  % read header
  eor=fread(fid,1,'int32');
  t=fread(fid,1,'float64');
  npart=fread(fid,1,'int32');
  eor=fread(fid,1,'int32');
  
  % read data
  for i=1:npart
    eor=fread(fid,1,'int32');
    part(i,:)=fread(fid,10,'float64');
    eor=fread(fid,1,'int32');
  end
  
  fclose(fid);
  
  figure(1)
  plot(mod(part(:,1),4*pi),part(:,2),'.')
  ylim([-1 1])
  xlim([0 4*pi])
  figure(2)
  hist(part(:,2),20)
  xlim([-1 1])
  ylim([0 200])
  
  
end