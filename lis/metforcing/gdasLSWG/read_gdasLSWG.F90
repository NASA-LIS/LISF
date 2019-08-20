!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: retgdasLSWG
! \label{retgdasLSWG}
!
! !REVISION HISTORY:
!  20 Oct 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine retgdasLSWG(n, order, m, metdata)
! !USES: 
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only : LIS_calendar, LIS_doy2date
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_verify
  use gdasLSWG_forcingMod, only : gdasLSWG_struc

  implicit none

  
! !ARGUMENTS:
!
  integer,   intent(in) :: n
  integer,   intent(in) :: order
  integer,   intent(in) :: m
  real                  :: metdata(LIS_rc%met_nf(m), LIS_rc%ngrid(n))

! !DESCRIPTION:
!
!   This subroutine reads the GDAS (ASCII) data specied for the LSWG 
!   work. At the start of the simulation, the code reads the entire dataset 
!   into memory and then indexes to the right location at other times. The
!   code also performs spatial interpolation to the LIS grid using 
!   a neighor search approach. 
!  
!   The data has 26 levels - We assume that there are 25 layers of data. 
!   The rh data is only specified for 21 levels. The remaining levels are
!   assumed to be zero. 
!EOP
  
  integer               :: ftn
  logical               :: file_exists
  integer               :: ios, iret
  integer               :: c,r, kk
  integer               :: yr, doy, mo, da, hr
  real                  :: lat, lon
  integer               :: status
  type(ESMF_Time)       :: time1 
  integer               :: i, iv, t,v
  logical               :: readflag
  integer               :: cindex, rindex, tindex
  real                  :: lis1d(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1             :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1             :: lb(gdasLSWG_struc(n)%mi)
!Indexed by param identifier (prefix p), layer type (t), and then actual layer(z)
  integer, parameter    :: num2ignore = 14
  real                  :: junk(num2ignore)
  real                  :: temp_tmp(26)
  real                  :: temp_rh(21)
  integer, parameter    :: zPLEVEL(26)= (/1000,975,950,925,900,850,800,750,700,650,600, &
       550,500,450,400,350,300,250,200,150,100,70,50,30,20,10/)

  if(gdasLSWG_struc(n)%startRead) then
     gdasLSWG_struc(n)%startRead = .false. 
     inquire(file=trim(gdasLSWG_struc(n)%gdasLSWGfile),exist=file_exists)
     
     if(.not.file_exists) then 
        write(LIS_logunit,*) 'GDASLSWG file',trim(gdasLSWG_struc(n)%gdasLSWGfile), &
             ' does not exist'
        write(LIS_logunit,*) 'Program stopping ..'
        call LIS_endrun()
     endif
     
     write(LIS_logunit,*) 'reading GDASLSWG file ',trim(gdasLSWG_struc(n)%gdasLSWGfile)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=trim(gdasLSWG_struc(n)%gdasLSWGfile))
     
     readflag = .true. 
     do while(readflag) 
        read(ftn,*,iostat=ios) yr, doy, hr, lat, lon, junk, temp_tmp, temp_rh 
        call LIS_doy2date(yr, doy, mo, da)
        call ESMF_TimeSet(time1, yy=yr, &
             mm = mo, dd = da, h = hr, calendar=LIS_calendar, &
             rc=status)
        call LIS_verify(status, 'ESMF_TimeSet in readGDASLSWG')
        t = nint((time1 - gdasLSWG_struc(n)%startTime)/gdasLSWG_struc(n)%timestep)+1
        
        cindex = nint((lon-gdasLSWG_struc(n)%gridDesci(5))/&
             gdasLSWG_struc(n)%gridDesci(9))+1
        rindex = nint((lat-gdasLSWG_struc(n)%gridDesci(4))/&
             gdasLSWG_struc(n)%gridDesci(10))+1
        tindex = rindex+(cindex-1)*nint(gdasLSWG_struc(n)%gridDesci(2))

        gdasLSWG_struc(n)%tmp(t,tindex,:) = temp_tmp(:)
        gdasLSWG_struc(n)%rh(t,tindex,:) = temp_rh(:)

        if(ios.ne.0) then 
           readflag = .false. 
        endif        
     enddo

     call LIS_releaseUnitNumber(ftn)
  endif

  if(order.eq.1) then 

     t = nint((gdasLSWG_struc(n)%btime1 - gdasLSWG_struc(n)%startTime)/&
          gdasLSWG_struc(n)%timestep)+1

  else
     t = nint((gdasLSWG_struc(n)%btime2 - gdasLSWG_struc(n)%startTime)/&
          gdasLSWG_struc(n)%timestep)+1
  endif

  ! Assuming that we will have no missing data
  lb = .true. 
 
! interpolating temperature values
  iv = 0
  do v=1,26
     iv = iv + 1
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,gdasLSWG_struc(n)%tmp(t,:,v),&
          lo,lis1d,gdasLSWG_struc(n)%mi, LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdasLSWG_struc(n)%n113,LIS_rc%udef,iret)

!     if(v==3) then 
!        open(100,file='varfield.bin',form='unformatted')
!        write(100) lis1d
!        close(100)
!        stop
!     endif

     do r=1, LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              metdata(iv,LIS_domain(n)%gindex(c,r)) = lis1d(c+(r-1)*LIS_rc%lnc(n))
           endif
        enddo
     enddo
  enddo
! interpolating rh values
  do v=1,21
     iv = iv+1
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,gdasLSWG_struc(n)%rh(t,:,v),&
          lo,lis1d,gdasLSWG_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdasLSWG_struc(n)%n113,LIS_rc%udef,iret)
     do r=1, LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              metdata(iv,LIS_domain(n)%gindex(c,r)) = lis1d(c+(r-1)*LIS_rc%lnc(n))
           endif
        enddo
     enddo
  enddo
! remaining rh levels specified to be 0
  do v=22,26
     iv = iv+1
     do r=1, LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              metdata(iv,LIS_domain(n)%gindex(c,r)) = 0.00
           endif
        enddo
     enddo
  enddo

!pressure values
  do v=1,26
     iv = iv+1
     do r=1, LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              metdata(iv,LIS_domain(n)%gindex(c,r)) = zPlevel(v)*100.0
           endif
        enddo
     enddo
  enddo

end subroutine retgdasLSWG

