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
!
! !ROUTINE: read_TRMM3B42RT
! \label{read_TRMM3B42RT}
!
! !REVISION HISTORY: 
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with LIS_domain 3 (2x2.5)
!
!
! !INTERFACE:
subroutine read_TRMM3B42RT (n, name_TRMM3B42RT, findex, order, ferror_TRMM3B42RT)
! !USES:
 use LIS_coreMod, only           : LIS_rc, LIS_domain
 use LIS_logMod, only            : LIS_logunit, LIS_getNextUnitNumber, &
                                   LIS_releaseUnitNumber
 use LIS_metforcingMod, only     : LIS_forc
 use TRMM3B42RT_forcingMod, only : TRMM3B42RT_struc
 use LIS_constantsMod,      only : LIS_CONST_PATH_LEN
 
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)    :: name_TRMM3B42RT
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_TRMM3B42RT

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  TRMM 3B42RT data and interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[name\_TRMM3B42RT]
!    name of the 3 hour TRMM 3B42RT forecast file
!  \item[ferror\_TRMM3B42RT]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_TRMM3B42RT](\ref{interp_TRMM3B42RT}) \newline
!    spatially interpolates the TRMM 3B42RT data
!  \end{description}
!EOP

  integer :: index1

!==== Local Variables=======================
       
 integer :: ios
 integer :: i,j,xd,yd
 parameter(xd=1440,yd=400)                  ! Dimension of TRMM 3B42RT data

 real :: precip(xd,yd), tmp(xd, yd)   
 real, allocatable :: precip_regrid(:,:)        ! Interpolated precipitation array
 character(len=LIS_CONST_PATH_LEN) :: fname                ! Filename variables
 character(len=LIS_CONST_PATH_LEN) :: dfile 
 integer             ::  ftn
 logical :: file_exists

!=== End Variable Definition =======================

 fname = name_TRMM3B42RT
!------------------------------------------------------------------------
! Fill necessary arrays to assure not using old TRMM 3B42RT data
!------------------------------------------------------------------------
! J.Case (4/22/2013) -- Make consistent with Stg4/NMQ routines
 if(order.eq.1) then 
     TRMM3B42RT_struc(n)%metdata1 = LIS_rc%udef ! J.Case
 elseif(order.eq.2) then 
     TRMM3B42RT_struc(n)%metdata2 = LIS_rc%udef ! J.Case
 endif
 allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
 precip_regrid = -1.0 ! J.Case

!------------------------------------------------------------------------
! Find TRMM 3B42RT precip data, read it in and assign to forcing precip array.
! Must reverse grid in latitude dimension to be consistent with LDAS grid
!------------------------------------------------------------------------
 dfile = trim(fname)//".bin.gz"
 inquire(file=dfile, EXIST=file_exists) 
 if (file_exists) then 
    call rd3B42RTgz(dfile, precip, xd, yd) 
 else 
    dfile = trim(fname)//".bin"
    inquire(file=dfile, EXIST=file_exists) 
    if (file_exists) then 
       call rd3B42RTbin(dfile, precip, xd, yd) 
    else 
       dfile = trim(fname)//".1gd4r" 
       inquire(file=dfile, EXIST=file_exists) 
       if (file_exists) then 
         call rd3B42RT1gd4r(dfile, precip, xd, yd) 
       else 
         dfile = trim(fname)//" [.bin.gz|.bin.|1gd4r] missing"
         precip = -1.0
       end if
   end if
 end if
    
 write(LIS_logunit,*) "3B42RT data: ", dfile 

 ! because raw RT data goes from 0.125 to 359.875, need to swap 
 ! western/eastern hemisphere to make it go -179.875 ~ 179.875

 Do j=1, yd
   Do i=1, 720
     tmp(i, j) = precip(i+720, j)
     precip(i+720, j) = precip(i, j)
     precip(i, j) = tmp(i, j)
  End Do 
End Do

! J.Case (4/19/2013) -- Test print out of raw precip array
! write (99,*) precip

!------------------------------------------------------------------------
! Interpolating to desired LIS_domain and resolution
! Global precip datasets not used currently to force NLDAS
!------------------------------------------------------------------------
   !write(LIS_logunit,*) "Writing un-interpolated TRMM 3B42RT precipitation out "
   !open(71, file="huff-ungrid.1gd4r", access="direct", &
   !    recl=xd*yd*4, form="unformatted")
   ! write(71, rec=1) precip
   !close(71)

    call interp_TRMM3B42RT(n, xd, yd, precip, LIS_rc%gridDesc(n,:), &
         LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid, findex)
    
   !write(LIS_logunit,*) "Writing interpolated TRMM 3B42RT precipitation out "
   !open(73, file="huff-regrid.1gd4r", access="direct", &
   !    recl=LIS_rc%d%lnr*LIS_rc%d%lnc*4, form="unformatted")
   ! write(73, rec=1) precip_regrid
   !close(73)
   !write(LIS_logunit,*) "Writing interpolated TRMM 3B42RT precipitation out finished"

! J.Case (4/19/2013) -- Test print out of the regridded precip (on LIS grid).
! write (98,*) precip_regrid
   
    do j = 1,LIS_rc%lnr(n)
       do i = 1,LIS_rc%lnc(n)
          if (precip_regrid(i,j) .ne. -1.0) then
             index1 = LIS_domain(n)%gindex(i,j)
             if(index1 .ne. -1) then
                if(order.eq.1) then 
                   TRMM3B42RT_struc(n)%metdata1(1,index1) = precip_regrid(i,j)   !here is mm/h
                elseif(order.eq.2) then 
                   TRMM3B42RT_struc(n)%metdata2(1,index1) = precip_regrid(i,j)   !here is mm/h
                endif
             endif
          endif
       enddo
    enddo

! J.Case (4/19/2013) -- Test print out of the suppdata precip (on LIS grid).
! write (97,*) TRMM3B42RT_struc(n)%metdata1(1,:)

   ferror_TRMM3B42RT = 1

 deallocate (precip_regrid)

end subroutine read_TRMM3B42RT


subroutine rd3B42RT1gd4r(dfile, precip, xd, yd)
 use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
      LIS_releaseUnitNumber

 integer :: i,j,xd,yd, ftn
 real :: precip(xd,yd)
 character(len=*) :: dfile 

 ftn = LIS_getNextUnitNumber()
 open(unit=ftn,file=dfile, status='old', &
      &          access='direct',recl=xd*yd*4, form='unformatted') 
   read (ftn,rec=1) precip
 close(ftn)
 call LIS_releaseUnitNumber(ftn)

 return 
end subroutine rd3B42RT1gd4r

!===============================
subroutine rd3B42RTbin(dfile, precip, xd, yd)
 use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
      LIS_releaseUnitNumber

 integer, parameter :: nc=1440, nr=480
 integer :: i,j,xd,yd, ftn
 real :: precip(xd,yd), output(nc, nr)
 character(len=*) :: dfile
 integer*2 :: rr(nc, nr)

 ftn = LIS_getNextUnitNumber()
 open(unit=ftn,file=dfile, status='old', &
      &          access='direct',recl=nc*2, form='unformatted') 
   Do j=1, nr
    read (ftn,rec=(j+1)) rr(:, j)  !skip 2880-byte header
   End Do 
 close(ftn)
 call LIS_releaseUnitNumber(ftn)


!   Convert integer to real, flip N-S, and set undef values
        Do j=1, nr
          Do i = 1, nc 
           if( rr(i, j).GE.0 ) then
             output(i, nr-j+1) = rr(i, j)*0.01
           else
             output(i, nr-j+1) = -9999.0
           end if
          End Do
        End Do

       Do j=41, 440
        precip(:, j-40) = output(:, j)     ! subset to 50 N/S
       End Do 

    return 

end subroutine rd3B42RTbin 

!============ Read .bin.gz file =================
subroutine rd3B42RTgz(zipfile, output, xd, yd)
        integer, parameter :: nc=1440, nr=480
        integer :: xd, yd
        real output(xd, yd)
        character(len=*) zipfile
        integer*2 input(nc, nr), itmp(nc, nr+1), rtmp(nc)  ! tmp includes header
        character*1 array(nc*(nr+1)*2), ct    ! buffer space
        equivalence (itmp, array)
        integer readzipf, dlen, rdlen, i, j, l

        dlen=nc*(nr+1)*2

        rdlen = readzipf(trim(zipfile)//char(0), array, dlen)
            ! swap the bytes for big-endian representation
            Do l=1, nc*(nr+1)*2-1, 2
             ct = array(l)
             array(l) = array(l+1)
             array(l+1) = ct
            End Do
            if(rdlen .ne. dlen) then
              write(*, *) "rdlen=", rdlen, " File reading error ..."
              write(*, *) "Fill array with undef"
              itmp = -999
            end if

        Do j=1, nr
          input(:, j) = itmp(:, j+1)
        End Do

!   Flip N-S
        Do j=41, 440
          Do i = 1, nc
           if( input(i, j).GE.0 ) then
             output(i, nr-j+1-40) = input(i, j)*0.01
           else
             output(i, nr-j+1-40) = -9999.0
           end if
          End Do
        End Do

        return
end subroutine rd3B42RTgz

