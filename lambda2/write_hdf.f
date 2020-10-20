c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine write_hdf(ur,ui,ar,ai,dur,dui,nx,nyp,nz,eta,xcoord,
     &     zcoord,t,fileout,scalar)
c
c     Creates HDF file
c
#ifdef HDF
      include 'hdf.inc'
#endif

      integer i,nx,nyp,nz, j, k,l
      real*4 data2(nx,nyp,nz)
      integer scalar
      real t
      real ur(nx/2+1,nyp,nz,3+scalar),ui(nx/2+1,nyp,nz,3+scalar)
      real ar(nx/2+1,nyp,nz,2),ai(nx/2+1,nyp,nz,2)
      real dur(nx/2+1,nyp,nz,3,3),dui(nx/2+1,nyp,nz,3,3)
      real eta(nyp), xcoord(nx), zcoord(nz)
      character*80 fileout

#ifdef HDF
      integer rank
      parameter(rank=3)
      integer sd_id, sds_id, status,dim_id
      integer sfstart, sfwdata, sfcreate, sfendacc, sfend, sfsdmname
      integer sfdimid, start(rank), sfsdscale, sfsnatt
      integer dimsizes(rank), stride(rank)
      character*10 vel_name(1:6), dim_name(rank)
c
c     Store names, coordinates, meta-data, etc
c
      dim_name(1) = 'x'
      dim_name(2) = 'y'
      dim_name(3) = 'z'
      dimsizes(1) = nx
      dimsizes(2) = nyp
      dimsizes(3) = nz
      vel_name(1) = 'u'
      vel_name(2) = 'v'
      vel_name(3) = 'w'
      vel_name(4) = 'umean'
      vel_name(5) = 'lambda2'
      vel_name(6) = 'dudy'
      stride(1) = 1
      stride(2) = 1
      stride(3) = 1
      start(1) = 0
      start(2) = 0
      start(3) = 0

      sd_id = sfstart(fileout,DFACC_CREATE)
      status = sfsnatt(sd_id, 'Time', DFNT_FLOAT, 1, real(t,4))
c
c     Write data, dimensions etc to datasets. The hdf-stride function
c     is not used because it is slow. Also, one cannot pass arguments by
c     reference and convert to single precision at the same time, hence
c     one has to rewrite the data to single precision separately.
c
      do i=1,6
         sds_id=sfcreate(sd_id,vel_name(i),DFNT_FLOAT, rank, dimsizes)
c
c     Assign name and scale to dimensions
c
         dim_id = sfdimid(sds_id,0)
         status = sfsdmname(dim_id,dim_name(1))
         status = sfsdscale(dim_id,dimsizes(1),
     &        DFNT_FLOAT, real(xcoord,4))
         dim_id = sfdimid(sds_id,1)
         status = sfsdmname(dim_id,dim_name(2))
         status = sfsdscale(dim_id,dimsizes(2),
     &        DFNT_FLOAT, real(eta,4))
         dim_id = sfdimid(sds_id,2)
         status = sfsdmname(dim_id,dim_name(3))
         status = sfsdscale(dim_id,dimsizes(3),
     &        DFNT_FLOAT, real(zcoord,4))

         if (i.le.3) then
            do l = 1,nz
               do k = 1,nyp
                  do j=1,nx/2
                     data2(2*j-1,k,l) = real(ur(j,k,l,i),4)
                     data2(2*j,k,l) = real(ui(j,k,l,i),4)
                  end do
               end do
            end do

            status = sfwdata(sds_id,start,stride,
     &           dimsizes,data2(1,1,1))
            if (status.eq.-1) write(*,*) 'WARNING: Writing failed'

         elseif (i.le.5) then
            do l = 1,nz
               do k = 1,nyp
                  do j=1,nx/2
                     data2(2*j-1,k,l) = real(ar(j,k,l,i-3),4)
                     data2(2*j,k,l) = real(ai(j,k,l,i-3),4)
                  end do
               end do
            end do

            status =sfwdata(sds_id,start,stride,
     &           dimsizes,data2(1,1,1))

         else
            do l = 1,nz
               do k = 1,nyp
                  do j=1,nx/2
                     data2(2*j-1,k,l) = real(dur(j,k,l,1,2),4)
                     data2(2*j,k,l) = real(dui(j,k,l,1,2),4)
                  end do
               end do
            end do

            status=sfwdata(sds_id,start,stride,
     &           dimsizes,data2(1,1,1))

         end if
         status = sfendacc(sds_id)
      end do

      status = sfend(sd_id)
#endif
      end subroutine write_hdf
