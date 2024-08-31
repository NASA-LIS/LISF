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
! !ROUTINE: LIS_process_cmd_args
! \label{LIS_process_cmd_args}
!
! !INTERFACE: 
subroutine LIS_process_cmd_args
! !USES: 
  use LIS_coreMod,   only : LIS_masterproc, LIS_rc
  use LIS_logMod,    only : LIS_endrun
#if ( defined AIX )
   use xlfutility, only : iargc
#endif
! !ARGUMENTS: 
!   None
!
! !DESCRIPTION: 
!    Processes command line arguments supplied to LIS at runtime.
!
! \begin{tabular}{ll}
! Option                 & Description                          \cr
! \texttt{-V} or 
! \texttt{--version}     & cause LIS to dump information
!                          about itself and then end.           \cr
! \texttt{-f <name>} or 
! \texttt{--file <name>} & specify the name of the LIS runtime
!                          configuration file.                  \cr
!                        & If this option is not specified,
!                          then LIS will default to reading
!                          ``lis.config''.
! \end{tabular}
!
!EOP

   implicit none

#if ( !(defined AIX) && !(defined GFORTRAN) )
   integer, external :: iargc
#endif

   integer :: i, numargs
   character(len=100) :: argi

   numargs=iargc()

   i = 1
   do
      if ( i > numargs ) then
         exit
      endif

      call getarg(i, argi) ; i = i + 1

      if ( trim(argi) == "--version" .or. trim(argi) == "-V" ) then
         if ( LIS_masterproc ) then
!            call LIS_version
            call LIS_endrun
         endif
      endif

      if ( trim(argi) == "--file" .or. trim(argi) == "-f" ) then
         call getarg(i, argi) ; i = i + 1
         LIS_rc%lis_config_file = trim(argi)
      endif
   enddo

end subroutine LIS_process_cmd_args
