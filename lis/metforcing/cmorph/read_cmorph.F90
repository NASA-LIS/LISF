!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_cmorph
! \label{read_cmorph}
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with LIS_domain 3 (2x2.5)
!  29 Dec 2003: Luis Goncalves; Added code to use CMORPH precip data
!  06 Jan 2006: Yudong Tian; modified for LISv4.2
!
! !INTERFACE:
subroutine read_cmorph (n, kk, name_cmorph, findex, order, ferror_cmorph, iflg )
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod,  only : LIS_logunit, LIS_getNextUnitNumber, &
                          LIS_releaseUnitNumber
  use LIS_metforcingMod, only : LIS_forc
  use cmorph_forcingMod, only : cmorph_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n 
  integer, intent(in) :: kk
  character(len=*)    :: name_cmorph
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_cmorph
  integer             :: iflg           
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  CMORPH data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[name\_cmorph]
!    name of the CMORPH file
!  \item[ferror\_cmorph]
!    flag to indicate success of the call (=0 indicates success)
!  \item[iflg]
!    flag indicating which 1/2 hour to read
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_cmorph](\ref{interp_cmorph}) \newline
!    spatially interpolates the CMORPH data
!  \end{description}
!EOP

  integer :: index1

!==== Local Variables=======================

  integer  :: r,ios
  integer  :: i,j,xd,yd,ibad
  parameter(xd=4948,yd=1649)                 ! Dimension of original CMORPH 8Km data
  integer, parameter :: ncmorph=xd*yd
  parameter(ibad=-9999.0)                    ! Bad (missing data) value
  character*1  precip(xd,yd),timestamp(xd,yd)
  character*1  staid(xd,yd)
  character*1  testout1(xd,yd)               ! Reconfigured original precip array
  real     :: realprecip(xd,yd)
  real     :: testout(xd,yd)
  real, allocatable :: precip_regrid(:,:)    ! Interpolated precip array
  character(len=LIS_CONST_PATH_LEN) :: fname, zname          ! Filename variables
  logical           :: file_exists
  integer           :: ftn

!=== End Variable Definition =======================

  allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  precip_regrid = -1.0
! J.Case (2/10/2015) -- I initialized realprecip to 0 in LISv6.1.  Is that necessary?
! I also initialized precip forcing to "LIS_rc%udef". Again, is that necessary?
  realprecip = -1.0
  fname = name_cmorph
  zname = trim(name_cmorph)//".Z"
  if(order.eq.1) then 
     cmorph_struc(n)%metdata1 = -1.0
  elseif(order.eq.2) then 
     cmorph_struc(n)%metdata2 = -1.0
  endif

!----------------------------------------------------------------------
! First look for and read compressed file (.Z). If not found, read binary file  
!----------------------------------------------------------------------
 ferror_cmorph = 1
 ! compressed file 
 inquire(file=zname, EXIST=file_exists)
 if (file_exists) then
    if(LIS_masterproc) write(LIS_logunit,*) &
       "[INFO] Reading Z compressed CMORPH precipitation data:: ",zname
    call rdCmorZ(zname, precip, xd, yd, iflg)
 else
   ! J.Case (2/10/2015) -- Added code to handle ".gz" CMORPH file extensions.
   zname = trim(name_cmorph)//".gz"
   inquire(file=zname, EXIST=file_exists)
   if (file_exists) then
     if(LIS_masterproc) write(LIS_logunit,*) &
       "[INFO] Reading gz compressed CMORPH precipitation data:: ",zname
     call rdCmorZ(zname, precip, xd, yd, iflg)
   else
     ! binary file 
     inquire(file=fname, EXIST=file_exists)
     if (file_exists) then
       if(LIS_masterproc) write(LIS_logunit,*) &
        "[INFO] Reading binary CMORPH precipitation data:: ",trim(fname)
       ftn = LIS_getNextUnitNumber()
       open(unit=ftn,file=fname, status='old',access='direct', &
            form='unformatted',recl=xd*yd*3,iostat=ios)
       read (ftn,rec=iflg) precip, timestamp, staid
       call LIS_releaseUnitNumber(ftn)
     else ! either file does not exist 
       if(LIS_masterproc) write(LIS_logunit,*) "[WARN] Missing CMORPH precipitation data ",trim(fname)
       ferror_cmorph = 0
     endif
   endif ! J.Case added  for .gz support (2/10/2015)
 end if

!    write(*,*) fname,ios,iflg
!debug 
!   open(25,file='./read_1',&
!   form='unformatted',access='direct',recl=xd*yd)
!   write(25,rec=1)precip
!   close(25)
!enddebug   
   
!    
! flipping data from 60N->60S to 60S->60N LIS standard    
!
    r=yd
    do i = 1,yd
       do j = 1,xd
         testout1(j,r) = precip(j,i)
       enddo
       r = r - 1
    enddo
    do i = 1,yd
       do j = 1,xd
         if(ichar(testout1(j,i)).eq.255)then
           realprecip(j,i) = -1
         else
           realprecip(j,i) = ichar(testout1(j,i))*0.2
         endif
       enddo
    enddo
!debug   
!   open(25,file='./read_2',&
!   access='direct',recl=xd*yd*4)
!   write(25,rec=1)realprecip
!   CLOSE(25)
!enddebug   
     
!    
! shifting data from 0->360 to -180->+180 LIS standard    
!
    do i = 1,yd
       r=1
       do j = (xd/2)+1,xd
         testout(r,i) = realprecip(j,i)
         r = r + 1
       enddo
       do j = 1,xd/2
         testout(r,i) = realprecip(j,i)
         r = r + 1
       enddo
    enddo 
    do i = 1,yd
       do j = 1,xd
         realprecip(j,i) = testout(j,i)
       enddo
    enddo
!debug   
!   open(25,file='./read_3',&
!   form='unformatted',access='direct',recl=xd*yd*4)
!   write(25,rec=1)realprecip
!   CLOSE(25)
!enddebug   
!
!=== End of data reconfiguration

    if(ferror_cmorph.eq.1) then 
       call interp_cmorph(n, xd, yd, realprecip, LIS_rc%gridDesc(n,:), &
            LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid)
       do j = 1,LIS_rc%lnr(n)
          do i = 1,LIS_rc%lnc(n)
             if (precip_regrid(i,j) .ge. 0.0) then
                index1 = LIS_domain(n)%gindex(i,j)
                if(index1 .ne. -1) then
                   if(order.eq.1) then 
                      cmorph_struc(n)%metdata1(kk,1,index1) = &
                           precip_regrid(i,j)   !here is mm/h
                   elseif(order.eq.2) then 
                      cmorph_struc(n)%metdata2(kk,1,index1) = &
                           precip_regrid(i,j)   !here is mm/h
                   endif
                endif
             endif
          enddo
       enddo
       if (LIS_masterproc) then 
         write(LIS_logunit,*)"[INFO] Obtained CMORPH precipitation data:: ",trim(fname)
       endif
    else
      if (LIS_masterproc) then
        write(LIS_logunit,*)"[WARN] Missing CMORPH precipitation data ",trim(fname)
      endif
    endif

    deallocate (precip_regrid)
    
 end subroutine read_cmorph

 subroutine rdCmorZ(zname, precip, xd, yd, iflg)

  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber

  integer :: xd, yd, iflg
  character*1 ::  precip(xd,yd)
  character(len=*) :: zname
  character*1, allocatable :: buff(:, :, :)
  integer :: readzipf, dlen, rdlen 
  
  allocate(buff(xd, yd, 6))   ! precip, time, satid, precp, time, satid
  dlen = xd*yd*6
! J.Case (2/10/2015) -- likely bug -- should be readzipf, not cm_readzipf
!  rdlen=cm_readzipf(trim(zname)//char(0), buff, dlen) 
  rdlen=readzipf(trim(zname)//char(0), buff, dlen) 
  if(LIS_masterproc) then
     write(LIS_logunit, *)"rdlen=", rdlen, " CMORPH-8KM Zipfile reading success"
  endif
  if (rdlen .NE. dlen) then
    if(LIS_masterproc) then
       write(LIS_logunit, *)"rdlen/dlen=",rdlen,dlen,&
      "CMORPH-8KM Zipfile reading error ... assign undef now"
    endif
    buff = char(255)
  end if
  precip = buff(:, :, (iflg-1)*3+1 )   ! precip, time, satid, precp, time, satid
  deallocate(buff)
  return

end subroutine rdCmorZ
