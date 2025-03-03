!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_gswp2
! \label{read_gswp2}
!  
! !REVISION HISTORY:
!
! 20Feb2004; Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_gswp2(order,n, findex, yr,mo,da,hr,mn,ss)
! !USES:
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_metforcingMod, only:LIS_forc
  use LIS_logMod, only : LIS_logunit
  use LIS_gswpMod, only : getgswp_timeindex
  use gswp2_forcingMod, only : gswp2_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  the GSWP2 files, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[yr]
!    current year
!  \item[mo]
!    current month
!  \item[da]
!    current day of the month
!  \item[hr]
!    current hour of day
!  \item[mn]
!    current minute
!  \item[ss]
!    current second
!   \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[getgswp\_timeindex](\ref{getgswp_timeindex}) \newline
!    computes the GSWP2 time index
!  \item[gswp2\_process\_forcing](\ref{gswp2_process_forcing}) \newline
!    reads the GSWP2 forcing from the native NETCDF file
!  \item[interp\_gswp2](\ref{interp_gswp2}) \newline
!    interpolates GSWP2 forcing to the LIS grid
!  \item[gswp2\_fillgaps](\ref{fillgaps_gswp2}) \newline
!    fills in the mismatches due to GSWP2 and LIS masks
!  \end{description}
!EOP
  integer :: yr,mo,da,hr,mn,ss
  integer :: tmp_yr, tmp_mo
  integer :: index1,c,r,order
  character*8 :: cyear,cmo
  character(len=LIS_CONST_PATH_LEN) :: ffile
  real :: fvar(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real, allocatable, dimension(:) :: fvar1
  real :: tempgswp2(LIS_rc%ngrid(n))
  integer :: tid,qid,pid,wid,swid,lwid,rid,sid
  integer :: ncid,status
  integer :: cindex,rindex
  integer :: land_id, land_len
  real,dimension(gswp2_struc(n)%ncold*gswp2_struc(n)%nrold) :: fvar2

  write(LIS_logunit,*)'read_gswp21',yr,mo,da,hr,mn,ss
  call getgswp_timeindex(yr,mo,da,hr,mn,ss,index1)
  write(LIS_logunit,*)'read_gswp22',yr,mo,da,hr,mn,ss,index1
#if ( defined GSWP2_OPENDAP )
  tmp = LIS_rc%gridDesc(4)/1000.0
  write(cslat, '(f8.2)') tmp
  tmp = LIS_rc%gridDesc(7)/1000.0
  write(cnlat, '(f8.2)') tmp
  tmp = LIS_rc%gridDesc(5)/1000.0
  write(cwlon, '(f8.2)') tmp
  tmp = LIS_rc%gridDesc(8)/1000.0
  write(celon, '(f8.2)') tmp
  write(c_index,'(i8)') index1
  
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading tair '
  call system("./gswp2_scripts/gettair.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/tair.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ! Note:  GSWP2 forcing data are written into monthly files with
  ! a frequency of 3 hours.  The first entry corresponds to 03:00:00
  ! of the first day of the given month.  The last entry corresponds
  ! to 00:00:00 of the first day of the next month.
  !
  ! E.g.; for Tair_cru198207.nc, the data run from 1982-07-01T03:00:00
  ! through 1982-08-01T00:00:00, inclusive.
  !
  ! So, when you are at hour 0 on the first day of the month,
  ! you need to open the file corresponding to the previous month.
  if ( da == 1 .and. hr == 0 .and. mn == 0 .and. ss == 0 ) then
     if ( mo == 1 ) then
        tmp_mo = 12
        tmp_yr = yr - 1
     else
        tmp_mo = mo - 1
        tmp_yr = yr
     endif
  else
     tmp_mo = mo
     tmp_yr = yr
  endif
  write(LIS_logunit,*) tmp_yr,tmp_mo,da,hr,mn,ss,index1
  write(cyear, '(I4)') tmp_yr
  write(cmo, '(I2.2)') tmp_mo
  ffile = trim(gswp2_struc(n)%tair)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading tair ',trim(ffile)

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_dimid(ncid, "land",land_id)
  status = nf90_inquire_dimension(ncid, land_id, len=land_len)
  status = nf90_close(ncid)
#endif
  gswp2_struc(n)%vector_len = land_len
  allocate(fvar1(land_len))

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Tair",tid)
  status = nf90_get_var(ncid, tid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex, gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,&
       fvar2,tempgswp2)

  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
!           if(LIS_domain(n)%gindex(cindex,rindex).lt.10) print*, 'forc ',cindex, rindex, fvar(cindex,rindex), lisdom(n)%gindex(cindex,rindex)
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

!  open(100,file='tair.bin',form='unformatted')
!  write(100) fvar
!  close(100)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(1,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif              
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(1,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(1,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo

#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading qair '
  call system("./gswp2_scripts/getqair.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/qair.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%qair)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading qair ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Qair",qid)
  status = nf90_get_var(ncid, qid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex1 = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(2,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(2,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(2,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(2,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(2,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo

#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading swdown '
  call system("./gswp2_scripts/getswdown.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/swdown.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%swdown)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading swdown ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "SWdown",swid)
  status = nf90_get_var(ncid, swid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(3,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(3,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(3,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(3,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(3,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo

#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading lwdown '
  call system("./gswp2_scripts/getlwdown.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/lwdown.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%lwdown)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading LWdown ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "LWdown",lwid)
  status = nf90_get_var(ncid, lwid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)


#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(4,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(4,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(4,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(4,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(4,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo
#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Wind '
  call system("./gswp2_scripts/getwind.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/wind.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%wind)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Wind ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Wind",wid)
  status = nf90_get_var(ncid, wid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)


#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(5,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(5,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(5,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(5,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(5,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo

#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Pressure '
  call system("./gswp2_scripts/getpsurf.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/psurf.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%psurf)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading psurf ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "PSurf",pid)
  status = nf90_get_var(ncid, pid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(7,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(7,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(7,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(7,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(7,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo

#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Rainf '
  call system("./gswp2_scripts/getrainf.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/rainf.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%rainf)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Rainf ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Rainf",rid)
  status = nf90_get_var(ncid, rid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(8,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(8,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(8,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(8,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(8,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo
#if ( defined GSWP2_OPENDAP )
#else
  ffile = trim(gswp2_struc(n)%snowf)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Snowf ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Snowf",sid)
  status = nf90_get_var(ncid, sid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(9,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(9,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(9,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(9,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(9,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo
        
#if ( defined GSWP2_OPENDAP )
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Rainf_C '
  call system("./gswp2_scripts/getrainf_c.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp2_scripts/rainf_c.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswp2_struc(n)%rainf_c)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  write(LIS_logunit,*)'MSG: GSWP2 forcing -- Reading Rainf_C ',trim(ffile)

#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Rainf_C",rid)
  status = nf90_get_var(ncid, rid,fvar1,&
       start = (/1,index1/),&
       count = (/land_len,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  tempgswp2 = -9999.0
  call gswp2_process_forcing(n, fvar1, fvar2)
  call interp_gswp2(n,findex,gswp2_struc(n)%ncold*gswp2_struc(n)%nrold,fvar2,tempgswp2)
  do c=1,LIS_rc%lnc(n)
     do r=1,LIS_rc%lnr(n)
     !do r=LIS_rc%lnr(n),1,-1
     !   rindex = 150-r+1
        rindex = r
        cindex = c
        if(LIS_domain(n)%gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = tempgswp2(LIS_domain(n)%gindex(cindex,rindex))
        endif
     enddo
  enddo

  call fillgaps_gswp2(n,1,fvar)

#endif
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(trim(LIS_rc%met_tinterp(findex)).eq."linear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(10,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif(order==2) then 
                 gswp2_struc(n)%metdata2(10,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           elseif(trim(LIS_rc%met_tinterp(findex)).eq."trilinear") then 
              if ( order == 1 ) then 
                 gswp2_struc(n)%metdata1(10,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              elseif ( order == 2 ) then 
                 gswp2_struc(n)%metdata2(10,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              else 
                 gswp2_struc(n)%metdata3(10,LIS_domain(n)%gindex(c,r)) = fvar(c,r)
              endif
           endif
        endif
     enddo
  enddo
  deallocate(fvar1)

end subroutine read_gswp2
!BOP
! 
! !ROUTINE: gswp2_process_forcing
! \label{gswp2_process_forcing}
! 
! !INTERFACE: 
subroutine gswp2_process_forcing(n, gswp21d, gswp2_filled)
! !USES: 
   use gswp2_forcingMod, only : gswp2_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                                        :: n
   real, dimension(gswp2_struc(n)%vector_len), intent(in)      :: gswp21d
   real, dimension(gswp2_struc(n)%ncold*gswp2_struc(n)%nrold), intent(out) :: gswp2_filled

   real, dimension(gswp2_struc(n)%ncold, gswp2_struc(n)%nrold) :: gswp22d

!
! !DESCRIPTION: 
!  
!  Converts the native GSWP2 vector data to a gridded data, reordered
!  to the LIS format. 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[gswp21d]
!    native GSWP2 data as a 1-d vector.
!  \item[gswp2\_filled]
!    Filled, 1-d vector GSWP2 field
!  \item[gswp22d]
!    reordered gridded GSWP2 field
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[gswp2\_expand\_to\_grid](\ref{gswp2_expand_to_grid}) \newline
!   expand the gswp2 native data, change to the LIS ordering style. 
!  \item[gspw\_2d\_to\_1d](\ref{gswp2_2d_to_1d}) \newline
!   convert the 2d field to a 1d vector
!  \end{description}
!EOP
   call gswp2_expand_to_grid(n, gswp21d, gswp22d)
   call gswp2_2d_to_1d(n, gswp22d, gswp2_filled)
   
end subroutine gswp2_process_forcing

!BOP
! 
! !ROUTINE: interp_gswp2
! \label{interp_gswp2}
! 
! !INTERFACE: 
subroutine interp_gswp2(n,findex, land_len,f,tempgswp2)
! !USES: 
   use LIS_coreMod,     only : LIS_rc, LIS_domain
   use gswp2_forcingMod, only : gswp2_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                        :: n 
   integer, intent(in)                        :: findex
   integer, intent(in)                        :: land_len 
   real, dimension(land_len), intent(in)      :: f
   real, dimension(LIS_rc%ngrid(n)), intent(out) :: tempgswp2
! 
! !DESCRIPTION: 
!  Interpolates the GSWP2 data onto the LIS grid. 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[land\_len]
!    size of the GSWP2 data
!  \item[f]
!    native GSWP2 data 
!  \item[tempgswp2]
!    interpolated GSWP2 data
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP
   real                      :: gridDesco(50)    ! Input,output grid info arrays
   integer                   :: count1,mo,iret,i,j
   logical*1                 :: lb(gswp2_struc(n)%ncold*gswp2_struc(n)%nrold)
   logical*1                 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   real                      :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))

!------------------------------------------------------------------------     
! Initializing input and output grid arrays
!------------------------------------------------------------------------
    gridDesco = 0
    gridDesco = LIS_rc%gridDesc(n,:)
    mo = LIS_rc%lnc(n) * LIS_rc%lnr(n)
!------------------------------------------------------------------------
! Defining input data bitmap
!------------------------------------------------------------------------
    lb = .true.
    do i = 1, land_len
      if ( f(i) == LIS_rc%udef ) then
          lb(i) = .false.
      endif
    enddo
!------------------------------------------------------------------------
! Defining output data bitmap
!------------------------------------------------------------------------
    lo = .true.
!------------------------------------------------------------------------
! Interpolate data from GSWP2 grid to LIS grid
!------------------------------------------------------------------------
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
       call bilinear_interp(gridDesco,lb,f,lo,go,gswp2_struc(n)%mi,mo, & 
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            gswp2_struc(n)%w111,gswp2_struc(n)%w121,&
            gswp2_struc(n)%w211,gswp2_struc(n)%w221,&
            gswp2_struc(n)%n111,gswp2_struc(n)%n121,&
            gswp2_struc(n)%n211,gswp2_struc(n)%n221,LIS_rc%udef, iret)
       
    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
       call conserv_interp(gridDesco,lb,f,lo,go,gswp2_struc(n)%mi,mo, & 
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            gswp2_struc(n)%w112,gswp2_struc(n)%w122,&
            gswp2_struc(n)%w212,gswp2_struc(n)%w222,&
            gswp2_struc(n)%n112,gswp2_struc(n)%n122,&
            gswp2_struc(n)%n212,gswp2_struc(n)%n222,LIS_rc%udef, iret)
    endif
    
    count1 = 0
    do j = 1, LIS_rc%lnr(n)
       do i = 1, LIS_rc%lnc(n)
          if(LIS_domain(n)%gindex(i,j) .ne. -1) then
             tempgswp2(LIS_domain(n)%gindex(i,j)) = go(i+count1)
          endif
       enddo
       count1 = count1 + LIS_rc%lnc(n)
    enddo
    
  end subroutine interp_gswp2

!BOP
! 
! !ROUTINE: gswp2_expand_to_grid
!  \label{gswp2_expand_to_grid}
! 
! !INTERFACE: 
subroutine gswp2_expand_to_grid(n, array1, array2)
! !USES: 
   use gswp2_forcingMod, only : gswp2_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                                         :: n
   real, dimension(gswp2_struc(n)%vector_len), intent(in)       :: array1
   real, dimension(gswp2_struc(n)%ncold, gswp2_struc(n)%nrold), &
         intent(out) :: array2
! 
! !DESCRIPTION: 
!  changes the 1d native GSWP2 data to a 2-d gridded data, reordered to 
!  N-S, E-W style.
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[array1]
!    native 1-d GSWP2 data
!  \item[array2]
!    reordered 2-d GSWP2 data. 
!  \end{description}
!EOP
   integer :: c, r, cindex, rindex

   array2 = -9999.0

   do c = 1, gswp2_struc(n)%ncold
      do r = 1, gswp2_struc(n)%nrold
         rindex = gswp2_struc(n)%nrold - r + 1
         cindex = c
         if ( gswp2_struc(n)%gindex(cindex,rindex) /= -1 ) then
            array2(cindex,rindex) = array1(gswp2_struc(n)%gindex(cindex,rindex))
         endif
      enddo
   enddo

end subroutine gswp2_expand_to_grid

!BOP
!
! !ROUTINE: gswp2_2d_to_1d
! \label{gswp2_2d_to_1d}
! 
! !INTERFACE: 
subroutine gswp2_2d_to_1d(n, array2, array1)
! !USES: 
   use gswp2_forcingMod, only : gswp2_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                                        :: n
   real, dimension(gswp2_struc(n)%ncold, gswp2_struc(n)%nrold),          &
         intent(in)                                           :: array2
   real, dimension(gswp2_struc(n)%ncold*gswp2_struc(n)%nrold),           &
         intent(out)                                          :: array1
! 
! !DESCRIPTION: 
!  reshapes the 2d-array to 1d array. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[array1]
!    input 2-d array
!  \item[array2]
!    input 1-d array
!  \end{description}
! 
!EOP
   integer :: c,r
!   array1 = reshape(array2, (/1/))
   do r=1,gswp2_struc(n)%nrold
      do c=1,gswp2_struc(n)%ncold
         array1(c+(r-1)*gswp2_struc(n)%ncold) = array2(c,r)
      enddo
   enddo
end subroutine gswp2_2d_to_1d

!BOP
! 
! !ROUTINE: gswp2_mask
! \label{gswp2_mask}
! 
! !INTERFACE: 
subroutine gswp2_mask()
#if ( defined USE_NETCDF3 ||  defined USE_NETCDF4 )
! !USES: 
   use LIS_coreMod,     only : LIS_rc
   use gswp2_forcingMod, only : gswp2_struc
   use netcdf
! 
! !DESCRIPTION: 
!   Reads the GSWP2 mask and computes the masking index.  
!
!EOP
   implicit none

   integer :: ncid, mvarid, status, ierr
   integer :: n, c, r, cindex, rindex, count1
   real, allocatable, dimension(:,:) :: mask

   do n = 1, LIS_rc%nnest

      allocate(mask(gswp2_struc(n)%ncold, gswp2_struc(n)%nrold), stat=ierr)
      allocate(gswp2_struc(n)%gindex(gswp2_struc(n)%ncold, gswp2_struc(n)%nrold), &
               stat=ierr)

      status = nf90_open(path=gswp2_struc(n)%mfile,mode=nf90_nowrite,ncid=ncid)
      status = nf90_inq_varid(ncid, "landmask",mvarid)
      status = nf90_get_var(ncid, mvarid,mask)
      status = nf90_close(ncid)

      gswp2_struc(n)%gindex = -1

      count1 = 1
      do c = 1, gswp2_struc(n)%ncold
         do r = 1, gswp2_struc(n)%nrold
            rindex = gswp2_struc(n)%nrold - r + 1
            cindex = c
            if ( mask(cindex,rindex) == 1.0 ) then
               gswp2_struc(n)%gindex(cindex,rindex) = count1
               count1 = count1 + 1
            endif
         enddo
      enddo

      deallocate(mask)
   enddo
#else
   use LIS_logMod, only : LIS_logunit, LIS_endrun
    write(LIS_logunit,*)'ERR: gswp2_mask requires NetCDF support'
   call LIS_endrun()
#endif
end subroutine gswp2_mask

