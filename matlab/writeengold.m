% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [status] = writeengold(casename,titlestring,NNx,NNy,NNz,xd,yFd,zd,Lx,Ly,Lz,cnf,X,Y,Z,vort,vortf,grad,gradp,l2,l2p)
%
% Write geometry and data files
%
%
% Write geometry file
%
% Note that big endian is not supported by the EnSight Gold format
  fwriteid = fopen([casename '.geo'],'w','ieee-le.l64');
%  fwriteid = fopen([casename '.geo'],'w');
  
  ts = textline('C Binary');
  fwrite(fwriteid,ts,'uchar');

  ts = textline('Structured DNS grid');
  fwrite(fwriteid,ts,'uchar');
  
  ts = textline(titlestring);
  fwrite(fwriteid,ts,'uchar');
  
  ts = textline('node id assign');
  fwrite(fwriteid,ts,'uchar');
  
  ts = textline('element id assign');
  fwrite(fwriteid,ts,'uchar');
  
% Optional
%  ts = textline('extents');
%  fwrite(fwriteid,ts,'uchar');
%  fwrite(fwriteid,[0 Lx 0 Ly 0 Lz],'int32');

  ts = textline('part');
  fwrite(fwriteid,ts,'uchar');
  
  fwrite(fwriteid,1,'int32');

  ts = textline('Main region');
  fwrite(fwriteid,ts,'uchar');
  %ts = textline('Fringe region');

  ts = textline('block rectilinear');
  fwrite(fwriteid,ts,'uchar');
  %
  % Write coordinate array lengths
  %
  fwrite(fwriteid,[NNx NNy NNz],'int32');
  %
  % Write coordinate arrays
  %
  fwrite(fwriteid,X(1,:,1),'float32');
  fwrite(fwriteid,Y(:,1,1),'float32');
  fwrite(fwriteid,Z(1,1,:),'float32');
      
  status = fclose(fwriteid);
  %
  % Write data file
  %
  fwriteid = fopen([casename '.vel'],'w');
  ts = textline('Velocity');
  fwrite(fwriteid,ts,'uchar');
  
  ts = textline('part');
  fwrite(fwriteid,ts,'uchar');
  
  fwrite(fwriteid,1,'int32');

  ts = textline('block');
  fwrite(fwriteid,ts,'uchar');
  %
  % Data are written in the order x, y, and z direction
  %
  %%surf(x,z,ecmf(:,:,iy)-mff(:,:,iy))
  for i=1:3
    for zz=1:NNz
      for yy=1:NNy
          fwrite(fwriteid,cnf(:,zz,yy+(i-1)*NNy),'float32');
      end
    end
  end
  status = fclose(fwriteid);
  
  %
  % Write velocity gradient data
  %
  if grad
    fwriteid = fopen([casename '.grad'],'w','ieee-le.l64');
    ts = textline('Gradient');
    fwrite(fwriteid,ts,'uchar');

    ts = textline('part');
    fwrite(fwriteid,ts,'uchar');

    fwrite(fwriteid,1,'int32');

    ts = textline('block');
    fwrite(fwriteid,ts,'uchar');
    %
    % Data are written in the order x, y, and z direction
    %
    for i=1:9
      for zz=1:NNz
	for yy=1:NNy
          fwrite(fwriteid,gradp(:,zz,yy+(i-1)*NNy),'float32');
	end
      end
    end
    status = fclose(fwriteid);
  end

  %
  % Write vorticity file
  %
  if vort
    fwriteid = fopen([casename '.vort'],'w','ieee-le.l64');
    ts = textline('Vorticity');
    fwrite(fwriteid,ts,'uchar');

    ts = textline('part');
    fwrite(fwriteid,ts,'uchar');

    fwrite(fwriteid,1,'int32');

    ts = textline('block');
    fwrite(fwriteid,ts,'uchar');
    %
    % Data are written in the order x, y, and z direction
    %
    for i=1:3
      for zz=1:NNz
	for yy=1:NNy
          fwrite(fwriteid,vortf(:,zz,yy+(i-1)*NNy),'float32');
	end
      end
    end
    status = fclose(fwriteid);
  end

  %
  % Write lambda2 file
  %
  if l2
    fwriteid = fopen([casename '.l2'],'w','ieee-le.l64');
    ts = textline('Lambda2');
    fwrite(fwriteid,ts,'uchar');

    ts = textline('part');
    fwrite(fwriteid,ts,'uchar');

    fwrite(fwriteid,1,'int32');

    ts = textline('block');
    fwrite(fwriteid,ts,'uchar');
    %
    % Data are written in the order x, y, and z direction
    %
    for zz=1:NNz
      for yy=1:NNy
	fwrite(fwriteid,l2p(:,zz,yy),'float32');
      end
    end
    status = fclose(fwriteid);
  end
