!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: get_gdasbc
! \label{get_gdasbc}
!
!
! !REVISION HISTORY:

! 
! !INTERFACE:
subroutine get_gdasbc(n,findex)
! !USES:
  use LIS_coreMod, only        : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_tick, LIS_get_nstep
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use gdasbc_forcingMod,  only : gdasbc_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!
!EOP

  integer :: c,f,ferror,try
  integer :: order
  real*8  :: time1,time2,timenow
  real*8  :: dtime1, dtime2
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  character(len=LIS_CONST_PATH_LEN) :: name
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime     ! 1=move time 2 data into time 1  
  integer :: kk           ! Forecast member index
  integer :: nstep
  
!=== End Variable Definition =============================================
  try=-999
  nstep = LIS_get_nstep(LIS_rc,n)  

!====Assumption will be not to find or move any data
  gdasbc_struc(n)%findtime1=0
  gdasbc_struc(n)%findtime2=0
  movetime=0

!=== Determine Required NLDAS-2 Data Times (The previous hour and the future hour)
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 

  yr1 = LIS_rc%yr  !previous assimilation/forecast hour
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = 6*(int(real(LIS_rc%hr)/6.0))
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call LIS_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr2 = LIS_rc%yr  !next assimilation/forecast hour
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 6*(int(real(LIS_rc%hr)/6.0))
  mn2 = 0
  ss2 = 0
  ts2 = 6*60*60
  call LIS_tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
  
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     gdasbc_struc(n)%findtime1=1
     gdasbc_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  
  if ( nstep == 0 .or. nstep == 1 ) then
     gdasbc_struc(n)%gdasbctime1 = time1
     gdasbc_struc(n)%gdasbctime2 = time2
  endif
  
  if ( timenow > gdasbc_struc(n)%gdasbctime2 ) then
     movetime  = 1
     gdasbc_struc(n)%findtime2 = 1
  end if
  
  if(movetime.eq.1) then
     gdasbc_struc(n)%gdasbctime1=gdasbc_struc(n)%gdasbctime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           gdasbc_struc(n)%metdata1(:,f,c)=gdasbc_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(gdasbc_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        do kk= gdasbc_struc(n)%st_iterid, gdasbc_struc(n)%en_iterid
           call get_gdasbc_filename(n,kk,findex,&
                name,gdasbc_struc(n)%gdasbcdir,&
                yr1,mo1,da1,doy1,hr1)
           write(unit=LIS_logunit,fmt=*)'[INFO] getting file1a.. ',trim(name)
           order = 1
           call read_gdasbc(n,kk,findex,order,mo1,name,ferror)
        enddo

        if(ferror.ge.1) gdasbc_struc(n)%gdasbctime1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-2 data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1	   	


  if(gdasbc_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

        do kk= gdasbc_struc(n)%st_iterid, gdasbc_struc(n)%en_iterid
           call get_gdasbc_filename(n,kk,findex, &
                name,gdasbc_struc(n)%gdasbcdir,&
                yr2,mo2,da2,doy2,hr2)
           write(unit=LIS_logunit,fmt=*)'[INFO] getting file2a.. ',trim(name)
           order = 2
           call read_gdasbc(n,kk,findex,order,mo2,name,ferror)
        end do
       
        if(ferror.ge.1) then
           gdasbc_struc(n)%gdasbctime2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-2 data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1

end subroutine get_gdasbc


 !BOP
! !ROUTINE: get_gdasbc_filename
! \label{get_gdasbc_filename}
!
! !REVISION HISTORY:
!
! !INTERFACE:
 subroutine get_gdasbc_filename(n,kk,findex, &
              filename,gdasbcdir,yr,mo,da,doy,hr)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_forecastMod

   implicit none
! !ARGUMENTS:
   integer                    :: n
   integer                    :: kk
   integer                    :: findex
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)  :: gdasbcdir
   integer, intent(in)        :: yr,mo,da,doy,hr
!
! !DESCRIPTION:
!   This subroutine puts together GES DISC NLDAS-2 "A" file name for
!   1 hour file intervals
!
!  The arguments are:
!  \begin{description}
!  \item[gdasbcdir]
!    Name of the NLDAS-2 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[doy]
!   Julian day of year (needed for subdirectory structure)
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the timestamped GES DISC NLDAS-2 file
!  \end{description}
!
!EOP
   character*4  :: fyr
   character*3  :: fdoy
   character*2  :: fmo, fda, fhr
   integer      :: doy2

  !=== end variable definition =============================================

   write(unit=fyr, fmt='(i4.4)')  yr
   write(unit=fdoy,fmt='(i3.3)')  doy
   write(unit=fmo, fmt='(i2.2)')  mo
   write(unit=fda, fmt='(i2.2)')  da
   write(unit=fhr, fmt='(i2.2)')  hr
   
    !=== Assemble GES DISC NLDAS-2 filename:

    ! NLDAS_FOR[A]0125_H.A19850601.0000.002.grb

   filename = trim(gdasbcdir)//"/"//fyr//fmo//"/GDAS_BC_"//&
        fyr//fmo//fda//fhr//"00.d01.nc"
   
 end subroutine get_gdasbc_filename

