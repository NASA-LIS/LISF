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
!               observations with LDT_domain 3 (2x2.5)
!  29 Dec 2003: Luis Goncalves; Added code to use CMORPH precip data
!  06 Jan 2006: Yudong Tian; modified for LDTv4.2
!
! !INTERFACE:
subroutine read_cmorph (n, name_cmorph, findex, order, ferror_cmorph, iflg )
! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod, only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_metforcingMod, only : LDT_forc
  use cmorph_forcingMod, only : cmorph_struc

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n 
  character(len=*)    :: name_cmorph
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_cmorph
  integer             :: iflg           
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  CMORPH data and interpolates to the LDT domain.
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

  integer   :: r,ios
  integer   :: i,j,xd,yd,ibad
  parameter(xd=4948,yd=1649)    ! Dimension of original CMORPH 8Km data
  integer, parameter :: ncmorph=xd*yd
  parameter(ibad=-9999.0)                        ! Bad (missing data) value
  character*1  precip(xd,yd),timestamp(xd,yd)
  character*1  staid(xd,yd)
  character*1  testout1(xd,yd)                   ! Reconfigured original precip array
  real      :: realprecip(xd,yd)
  real      :: testout(xd,yd)
  real, allocatable :: precip_regrid(:,:)        ! Interpolated precip array
  character(len=LDT_CONST_PATH_LEN) :: fname, zname              ! Filename variables
  logical           :: file_exists
  integer           :: ftn
!=== End Variable Definition =======================

  allocate (precip_regrid(LDT_rc%lnc(n),LDT_rc%lnr(n)))
  precip_regrid = -1.0
  realprecip = -1.0
  fname = name_cmorph
  zname = trim(name_cmorph)//".Z"
  if(order.eq.1) then 
     LDT_forc(n,findex)%metdata1 = -1.0
  elseif(order.eq.2) then 
     LDT_forc(n,findex)%metdata2 = -1.0
  endif
!----------------------------------------------------------------------
! First look for and read compressed file (.Z). If not found, read binary file  
!----------------------------------------------------------------------
 ferror_cmorph = 1
 ! compressed file 
 inquire(file=zname, EXIST=file_exists)
 if (file_exists) then
    if(LDT_masterproc) write(LDT_logunit,*) "Reading CMORPH precipitation data:: ", trim(zname)
    call rdCmorZ(zname, precip, xd, yd, iflg)
 else
    ! binary file 
    inquire(file=fname, EXIST=file_exists)
    if (file_exists) then
       if(LDT_masterproc) write(LDT_logunit,*) "Reading CMORPH precipitation data:: ", trim(fname)
       ftn = LDT_getNextUnitNumber()
       open(unit=ftn,file=fname, status='old',access='direct', &
            form='unformatted',recl=xd*yd*3,iostat=ios)
       read (ftn,rec=iflg) precip, timestamp, staid
       call LDT_releaseUnitNumber(ftn)
    else ! either file does not exist 
       if(LDT_masterproc) write(LDT_logunit,*) "Missing CMORPH precipitation data ", trim(fname) 
       ferror_cmorph = 0
    endif
  end if

!    write(*,*) fname,ios,iflg
!debug 
!   open(25,file='./read_1',&
!   form='unformatted',access='direct',recl=xd*yd)
!   write(25,rec=1)precip
!   close(25)
!enddebug   
   
!    
! flipping data from 60N->60S to 60S->60N LDT standard    
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
! shifting data from 0->360 to -180->+180 LDT standard    
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
       call interp_cmorph(n, xd, yd, realprecip, LDT_rc%gridDesc(n,:), &
            LDT_rc%lnc(n),LDT_rc%lnr(n),precip_regrid)
       do j = 1,LDT_rc%lnr(n)
          do i = 1,LDT_rc%lnc(n)
             if (precip_regrid(i,j) .ge. 0.0) then
                index1 = LDT_domain(n)%gindex(i,j)
                if(index1 .ne. -1) then
                   if(order.eq.1) then 
                      LDT_forc(n,findex)%metdata1(1,index1) = &
                           precip_regrid(i,j)   !here is mm/h
                   elseif(order.eq.2) then 
                      LDT_forc(n,findex)%metdata2(1,index1) = &
                           precip_regrid(i,j)   !here is mm/h
                   endif
                endif
             endif
          enddo
       enddo
       write(LDT_logunit,*) "Obtained CMORPH precipitation data:: ", trim(fname)
    else
       write(LDT_logunit,*) "Missing CMORPH precipitation data ", trim(fname)
    
    endif

    deallocate (precip_regrid)
    
 end subroutine read_cmorph


 subroutine rdCmorZ(zname, precip, xd, yd, iflg)

  use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod, only : LDT_logunit, LDT_getNextUnitNumber

  integer      :: xd, yd, iflg
  character*1  :: precip(xd,yd)
  character(len=*) :: zname
  character*1, allocatable :: buff(:, :, :)
  integer      :: cm_readzipf, dlen, rdlen 
  
  allocate(buff(xd, yd, 6))   ! precip, time, satid, precp, time, satid
  dlen = xd*yd*6
  rdlen=cm_readzipf(trim(zname)//char(0), buff, dlen) 
  if(LDT_masterproc) &  
    write(LDT_logunit, *)"rdlen=", rdlen, " CMORPH-8KM Zipfile reading success"
  if (rdlen .NE. dlen) then
    if(LDT_masterproc) &  
      write(LDT_logunit, *)"rdlen=", rdlen, " CMORPH-8KM Zipfile reading error ... assign undef now"
      buff = char(255)
  end if
  precip = buff(:, :, (iflg-1)*3+1 )   ! precip, time, satid, precp, time, satid
  deallocate(buff)
  return

 end subroutine rdCmorZ
