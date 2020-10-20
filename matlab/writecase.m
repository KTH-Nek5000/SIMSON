% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [status] = writecase(casename,vel,NNx,NNy,NNz,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv,nsteps,nstart,ninc,vort,grad,l2)
%
% This function creates a EnSight Gold case file
%
  timeset=1;
  flag = 0;
  fid = fopen([casename '.case'],'w');
  %
  % Write format section
  %
  fprintf(fid,'FORMAT\n');
  fprintf(fid,'type: ensight gold\n');
  fprintf(fid,'\n');
  %
  % Write geometry section
  %
  fprintf(fid,'GEOMETRY\n');
  fprintf(fid,['model: 1  ' casename '.geo\n']);
  fprintf(fid,'\n');
  %
  % Write variable section
  %
  fprintf(fid,'VARIABLE\n');
  fprintf(fid,'constant per case:  1 X-length   %1.8f\n',Lx);
  fprintf(fid,'constant per case:  1 Y-length   %1.8f\n',Ly);
  fprintf(fid,'constant per case:  1 Z-length   %1.8f\n',Lz);
  fprintf(fid,'constant per case:  1 Re         %1.8f\n',Re);
  fprintf(fid,'constant per case:  1 Nx         %1.0f\n',NNx);
  fprintf(fid,'constant per case:  1 Ny         %1.0f\n',NNy);
  fprintf(fid,'constant per case:  1 Nz         %1.0f\n',NNz);
  fprintf(fid,'constant per case:  1 Flowtype   %1.0f\n',flowtype);

  fprintf(fid,'constant per case:  1 Rlam       %1.8f\n',rlam);
  fprintf(fid,'constant per case:  1 Spanv      %1.8f\n',spanv);
  fprintf(fid,['vector per node:    1 Velocity   ' casename '.vel\n']);
  if grad
      fprintf(fid,['tensor asym per node:    1 Gradient   ' casename '.grad\n']);
  end
  if l2
      fprintf(fid,['scalar per node:    1 Lambda2   ' casename '.l2\n']);
  end
  if vort
      fprintf(fid,['vector per node:    1 Vorticity   ' casename '.vort\n']);
  end
  %
  % Write time section
  %
  fprintf(fid,'\n');
  fprintf(fid,'TIME\n');
  fprintf(fid,'time set:                     %1.0f\n',timeset);
  fprintf(fid,'number of steps:              %1.0f\n',nsteps);
  fprintf(fid,'filename start number:        %1.0f\n',nstart);
  fprintf(fid,'filename increment:           %1.0f\n',ninc);
  fprintf(fid,'time values:                  %1.0f\n',t(1));
  for j=2:nsteps
    fprintf(fid,'%1.8f',t(j));
    if ( mod(j,6) == 0 )
      fprintf(fid,'\n');
    end
  end
  % 
  % Close the file
  %
  status = fclose(fid);
