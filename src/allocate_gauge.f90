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
      real,allocatable   ::  uparea(:,:)                 !! drainage area (GRID base)
      integer,allocatable::  basin(:,:)                  !! next grid X
      integer,allocatable::  biftag(:,:)                  !! next grid X
      integer,allocatable::  class(:,:)                  !! flow direction conbined
      integer,allocatable::  nextXX(:,:), nextYY(:,:)
      integer,allocatable::  mask(:,:)

      real,allocatable   ::  glon(:), glat(:)
      real               ::  cnt_lon, cnt_lat, lon1, lat1, lon2, num
      integer            ::  imask

      integer*8          ::  id, bsn
      real               ::  lat0, lon0, lat, lon, area0, area
      real               ::  err, err0, err1, err2, dd, diff
      real               ::  rate, rate0
      integer            ::  kXX, kYY, tag, nn
      integer            ::  iXX0, iYY0

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

      allocate(uparea(nXX,nYY), biftag(nXX,nYY), class(nXX,nYY), basin(nXX,nYY))
      allocate(nextXX(nXX,nYY),nextYY(nXX,nYY),  mask(nXX,nYY))

      rfile1='./'//trim(regmap)//'/uparea.bin'
      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nXX*nYY,status='old',iostat=ios)
      read(11,rec=1) uparea
      close(11)

      rfile1='./'//trim(regmap)//'/basin.bin'
      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nXX*nYY,status='old',iostat=ios)
      read(11,rec=1) basin
      close(11)

      rfile1='./'//trim(regmap)//'/nextxy.bin'
      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nXX*nYY,status='old',iostat=ios)
      read(11,rec=1) nextXX
      read(11,rec=2) nextYY
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
      call getarg(1,rlist)
      open(11, file=rlist, form='formatted')
      read(11,*)

      wfile1='./gauge_alloc.txt'
      open(21, file=wfile1, form='formatted')
      write(21,'(a)') '  ID  basin  lat_ori  lon_ori  area_ori  lat_MERIT  lon_MERIT  area_MERIT  diff  error  ix  iy  bsn_lat  bsn_lon'

 1000 continue
      read(11,*,end=1090) id, lat0, lon0, area0
      print *, id, lat0, lon0, area0

      iXX=int( (lon0 -west)/gsize )+1
      iYY=int( (north-lat0)/gsize )+1

      err0=1.e20
      err1=1.e20
      rate0=1.e20
      tag=0

      nn=5       !!  search domain
      if( area0>1000    ) nn=6
      if( area0>10000   ) nn=7
      if( area0>100000  ) nn=8
      if( area0>300000  ) nn=9
      if( area0>1000000 ) nn=10

      do dYY=-nn, nn
        do dXX=-nn, nn
          jXX=iXX+dXX
          jYY=iYY+dYY
          if( jXX<=0 ) jXX=jXX+nXX
          if( jXX>nXX) jXX=jXX-nXX

          if( jYY>0 .and. jYY<=nYY )then
            if( uparea(jXX,jYY)>area0*0.05 )then
              err=(uparea(jXX,jYY)-area0)/area0
              err2=err
              dd=(  (jYY-iYY)**2.+(jXX-iXX)**2. )**0.5

!              if( err>0 ) err2=err+0.02*dd
!              if( err<0 ) err2=err-0.02*dd
              if( err>0 ) err2=err+0.2*dd
              if( err<0 ) err2=err-0.2*dd


              if( err2>=0 )then
                rate=(1+err2)
              elseif( err2>-1 .and. err2<0 )then
                rate=1./(1+err2)
                rate=min(rate,1000.)
              else
                rate=1000
              endif

              if( rate<rate0 )then
                err0=err2
                err1=err
                kXX=jXX
                kYY=jYY
                area=uparea(kXX,kYY)
                lon=west +gsize*(kXX-0.5)
                lat=north-gsize*(kYY-0.5)

                if( err0>=0 )then
                  rate0=(1+err0)
                elseif( err0>-1 .and. err0<0 )then
                  rate0=1/(1+err0)
                  rate0=min(rate0,1000.)
                else
                  rate0=1000
                endif

              endif
            endif
          endif
        end do
      end do
      if( err0<1.e20 .and. area0>0 )then
        iXX0=kXX
        iYY0=kYY
        diff=area-area0
        bsn=basin(iXX0,iYY0)

        mask(:,:)=0
        do iYY=1, nYY
          do iXX=1, nXX
            if( basin(iXX,iYY)/=bsn )then
              mask(iXX,iYY)=-9
            else
              mask(iXX,iYY)=0
              if( nextXX(iXX,iYY)<0 ) mask(iXX,iYY)=2
            endif
          end do
        end do
        mask(iXX0,iYY0)=1

        do iYY=1, nYY
          do iXX=1, nXX
            if( mask(iXX,iYY)==0 )then
              jXX=iXX
              jYY=iYY
              do while( mask(jXX,jYY)==0 )
                kXX=nextXX(jXX,jYY)
                kYY=nextYY(jXX,jYY)
                jXX=kXX
                jYY=kYY
              end do
              imask=mask(jXX,jYY)

              jXX=iXX
              jYY=iYY
              do while( mask(jXX,jYY)==0 )
                mask(jXX,jYY)=imask
                kXX=nextXX(jXX,jYY)
                kYY=nextYY(jXX,jYY)
                jXX=kXX
                jYY=kYY
              end do
            endif
          end do
        end do

        lon1=glon(iXX0)
        lat1=glat(iYY0)
        cnt_lat=0
        cnt_lon=0
        num=0
        do iYY=1, nYY
          do iXX=1, nXX
            if( mask(iXX,iYY)==1) then
              num=num+1

              cnt_lat=cnt_lat+glat(iYY)

              lon2=glon(iXX)
              if( lon1> 120 .and. lon2<-120 ) lon2=lon2+360
              if( lon1<-120 .and. lon2> 120 ) lon2=lon2-360
              cnt_lon=cnt_lon+lon2
            endif
          end do
        end do
        cnt_lat=cnt_lat/num
        cnt_lon=cnt_lon/num
        if( cnt_lon<-180 ) cnt_lon=cnt_lon+360
        if( cnt_lon> 180 ) cnt_lon=cnt_lon-360

        write(21,'(i16,i8, 7f15.4,f10.4,2i6,2f10.4)') id, bsn, lat0, lon0, area0, lat, lon, area, &
          diff, err1, kXX, kYY, cnt_lat, cnt_lon
      else
        kXX=-999
        kYY=-999
        lon=0
        lat=0
        area=0
        err0=-9
        err1=-9
        diff=0
        bsn=-9
        write(21,'(i16,i8, 7f15.4,f10.4,2i6)') id, bsn, lat0, lon0, area0, lat, lon, area, diff, err1, kXX, kYY
        write(6,'(i16,i8, 7f15.4,f10.4,2i6)') id, bsn, lat0, lon0, area0, lat, lon, area, diff, err1, kXX, kYY, -9, -9
      endif

      goto 1000
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


