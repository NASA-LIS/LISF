!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! !ROUTINE: rdhm356_temper_file
! \label{rdhm356_temper_file}
!
! !INTERFACE:
subroutine rdhm356_temper_file( name, rdhm356_temper_dir, yr, mo, da, hr)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE2 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[rdhm356\_temper\_dir]
!    Name of the STAGE II directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name]
!   name of the time-stamped STAGE II file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, mo, da, hr

  character(len=*) :: name
  character(40) :: rdhm356_temper_dir
  character(4) :: cyear
  character(2) :: cmon, cday, chour

!- Build the filename to be opened

   write ( cyear, '(i4)' ) yr
   if ( mo < 10 ) write ( cmon, '(a1,i1)' ) "0", mo
   if ( mo >= 10) write ( cmon, '(i2)' ) mo
   if ( da < 10 ) write ( cday, '(a1,i1)' ) "0", da
   if ( da >= 10) write ( cday, '(i2)' ) da
   if ( hr < 10 ) write ( chour, '(a1,i1)' ) "0", hr
   if ( hr >= 10) write ( chour, '(i2)' ) hr

!  tair0101198701z.gz
   name = trim(rdhm356_temper_dir)//'/tair'//cmon//cday//cyear//chour//'z'//'.gz'
end subroutine rdhm356_temper_file

!BOP
! !ROUTINE: rdhm356_precip_file
! \label{rdhm356_precip_file}
!
! !INTERFACE:
subroutine rdhm356_precip_file( name, rdhm356_precip_dir, yr, mo, da, hr)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE2 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[rdhm356\_precip\_dir]
!    Name of the STAGE II directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name]
!   name of the time-stamped STAGE II file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, mo, da, hr

  character(len=*) :: name
  character(40) :: rdhm356_precip_dir
  character(4) :: cyear
  character(2) :: cmon, cday, chour

!- Build the filename to be opened

   write ( cyear, '(i4)' ) yr
   if ( mo < 10 ) write ( cmon, '(a1,i1)' ) "0", mo
   if ( mo >= 10) write ( cmon, '(i2)' ) mo
   if ( da < 10 ) write ( cday, '(a1,i1)' ) "0", da
   if ( da >= 10) write ( cday, '(i2)' ) da
   if ( hr < 10 ) write ( chour, '(a1,i1)' ) "0", hr
   if ( hr >= 10) write ( chour, '(i2)' ) hr

!  xmrg1109200913z.gz
   name = trim(rdhm356_precip_dir)//'/xmrg'//cmon//cday//cyear//chour//'z'//'.gz'
end subroutine rdhm356_precip_file

