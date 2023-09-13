      program SET_MAP
! ===============================================
! to set various maps
! ===============================================
      implicit none
! index TRIP
      integer            ::  iXX, iYY, jXX, jYY, dXX, dYY
      integer            ::  nXX, nYY                        !! x-y matrix GLOBAL
      real               ::  gsize                           !! grid size [degree]
      real               ::  west, north, east, south
      integer            ::  ios
! input
      real,allocatable   ::  uparea(:,:)
      real,allocatable   ::  rivwth(:,:)
      real,allocatable   ::  elevtn(:,:)

      real,allocatable   ::  glon(:), glat(:)


      integer*8          ::  id
      real               ::  lat0, lon0, egm96
      real               ::  lat,  lon,  area,  width, elv

      real               ::  upa_thrs, upa_min, upa_max, wth_thrs
      integer            ::  isFound, istep, nn

      integer            ::  kXX, kYY

! file
      character*128      ::  regmap
      character*128      ::  rfile1, rlist
      character*128      ::  wfile1
! ===============================================
      west=-180
      east=180
      south=-90
      north=90

      nXX=21600
      nYY=10800
      gsize=1./60.

      regmap='./1min_river/'

      allocate(uparea(nXX,nYY),rivwth(nXX,nYY),elevtn(nXX,nYY))

      rfile1='./'//trim(regmap)//'/uparea.bin'
      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nXX*nYY,status='old',iostat=ios)
      read(11,rec=1) uparea
      close(11)

      rfile1='./'//trim(regmap)//'/width.bin'
      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nXX*nYY,status='old',iostat=ios)
      read(11,rec=1) rivwth
      close(11)

      rfile1='./'//trim(regmap)//'/elevtn.bin'
      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nXX*nYY,status='old',iostat=ios)
      read(11,rec=1) elevtn
      close(11)

      do iYY=1, nYY
        do iXX=1, nXX
          if( uparea(iXX,iYY)>0 ) uparea(iXX,iYY)=uparea(iXX,iYY)*1.e-6
        end do
      end do

      allocate(glon(nXX),glat(nYY))
      do iYY=1, nYY
        glat(iYY)=  90.0-gsize*(iYY-0.5)
      end do
      do iXX=1, nXX
        glon(iXX)=-180.0+gsize*(iXX-0.5)
      end do


! ===============================================
      rlist='./HydroWeb_input.txt'
      open(11, file=rlist, form='formatted')
      read(11,*)

      wfile1='./HydroWeb_alloc.txt'
      open(21, file=wfile1, form='formatted')
      write(21,'(a)') '  ID  lat_ori  lon_ori  lat_MERIT  lon_MERIT  area_MERIT  width_MERIT  elv_MERIT'

      upa_thrs=10000.
      upa_min =1000.

 1000 continue
      read(11,*,end=1090) id, lat0, lon0, egm96
!!      print *, id, lat0, lon0, egm96

      iXX=int( (lon0 -west)/gsize )+1
      iYY=int( (north-lat0)/gsize )+1

      isFound=0
      upa_max=0
      wth_thrs=1

! ===================================
! large river with width data

      if( uparea(iXX,iYY)>upa_min .and. rivwth(iXX,iYY)>=wth_thrs )then
        upa_max=uparea(iXX,iYY)
        kXX=iXX
        kYY=iYY
        isFound=1
      endif

      if( upa_max<upa_thrs )then  !! seach larger stream nearby
        nn=1
        do dYY=-nn, nn
          do dXX=-nn, nn
            jXX=iXX+dXX
            jYY=iYY+dYY
            if( jXX<=0 ) jXX=jXX+nXX
            if( jXX>nXX) jXX=jXX-nXX
            if( jYY>0 .and. jYY<=nYY )then
              if( uparea(jXX,jYY)>upa_max .and. uparea(jXX,jYY)>upa_min .and. rivwth(jXX,jYY)>=wth_thrs )then  !! larger than 2x local river
                upa_max=uparea(jXX,jYY)
                kXX=jXX
                kYY=jYY
                isFound=1
              endif
            endif
          end do
        end do
      endif

      if( isFound==0 .or. upa_max<upa_thrs )then  !! still not found
        nn=2
        do dYY=-nn, nn
          do dXX=-nn, nn
            jXX=iXX+dXX
            jYY=iYY+dYY
            if( jXX<=0 ) jXX=jXX+nXX
            if( jXX>nXX) jXX=jXX-nXX
            if( jYY>0 .and. jYY<=nYY )then
              if( uparea(jXX,jYY)>upa_max .and. uparea(jXX,jYY)>upa_min .and. rivwth(jXX,jYY)>=wth_thrs )then  !! larger than 2x local river
                upa_max=uparea(jXX,jYY)
                kXX=jXX
                kYY=jYY
                isFound=1
              endif
            endif
          end do
        end do
      endif

! =========================================
!  allow smaller rivers

      if( isFound==0 )then
        istep=1
        upa_min=500
        wth_thrs=1

        nn=0
 2100   continue
        do dYY=-nn, nn
          do dXX=-nn, nn
            jXX=iXX+dXX
            jYY=iYY+dYY
            if( jXX<=0 ) jXX=jXX+nXX
            if( jXX>nXX) jXX=jXX-nXX
            if( jYY>0 .and. jYY<=nYY )then
              if( uparea(jXX,jYY)>upa_min .and. rivwth(jXX,jYY)>=wth_thrs )then  !! larger than 2x local river
                if( uparea(jXX,jYY)>upa_max )then
                  upa_max=uparea(jXX,jYY)
                  kXX=jXX
                  kYY=jYY
                  isFound=1
                endif
              endif
            endif
          end do
        end do
        if( isFound==0 )then
          if( nn<5 )then
            nn=nn+1
            goto 2100
          else
            if( istep==1 )then
              istep=istep+1
              upa_min=500
              wth_thrs=0
              nn=0
              goto 2100
            elseif( istep==2 )then
              istep=istep+1
              upa_min=100
              wth_thrs=1
              nn=0
              goto 2100
            elseif( istep==3 )then
              istep=istep+1
              upa_min=1000
              wth_thrs=0
              nn=0
              goto 2100
            elseif( istep==4 )then
              istep=istep+1
              upa_min=500
              wth_thrs=0
              nn=0
              goto 2100
            elseif( istep==5 )then
              istep=istep+1
              upa_min=100
              wth_thrs=0
              nn=0
              goto 2100
            endif
          endif
        endif
      endif

      if( isFound==0 )then
        lon=-999
        lat=-999
        elv=-999
        area=-999
        width=-999
        write(6,'(a12,i16,6f10.2)') 'Not Found: ', id, lat0, lon0
        write(21,'(i16,7f15.4)') id, lat0, lon0, lat, lon, area, width, elv
        goto 1000
      else
        lon=glon(kXX)
        lat=glat(kYY)
        elv=elevtn(kXX,kYY)
        area=uparea(kXX,kYY)
        width=rivwth(kXX,kYY)
        write(21,'(i16,7f15.4)') id, lat0, lon0, lat, lon, area, width, elv
        goto 1000
      endif
 1090 continue

      close(11)
      close(21)

! ====================
      end program SET_MAP






       subroutine nextGRID(iXX, iYY, jXX, jYY)
! ===============================================
       implicit none
       integer            ::  iXX, iYY, jXX, jYY
! ===============================================
       iXX=jXX
       iYY=jYY
       return
       end


