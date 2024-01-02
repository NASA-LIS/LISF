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
! !ROUTINE: RAPID_routing_readrst
! \label{RAPID_routing_readrst}
!
! !REVISION HISTORY:
! 19 Jul 2021: Yeosang Yoon;  Initial implementation

subroutine RAPID_routing_readrst

  use ESMF
  use LIS_fileIOMod
  use LIS_coreMod
  use LIS_logMod
  use RAPID_routingMod, only : RAPID_routing_struc
  use netcdf

  implicit none

  integer       :: n 
  integer       :: ftn
  integer       :: status
  character*100 :: filename
  logical       :: read_restart
  integer       :: varid_Qout

  do n=1, LIS_rc%nnest 

     read_restart = .false.
     if(RAPID_routing_struc(n)%startmode.eq."restart") then
        read_restart = .true.
     endif

     if(read_restart) then
        if(LIS_masterproc) then   
           write(LIS_logunit,*) '[INFO] RAPID restart file used: ', &
                                 trim(RAPID_routing_struc(n)%rstfile)

           status = nf90_open(path=RAPID_routing_struc(n)%rstfile,mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status, "Error opening file "//trim(RAPID_routing_struc(n)%rstfile))

           ! Get the varid of the data variable, based on its name.
           call LIS_verify(nf90_inq_varid(ftn,"Qout",varid_Qout),&
                'Error in nf90_inq_varid in RAPID_routing_readrst')

           ! Read the data.
           call LIS_verify(nf90_get_var(ftn,varid_Qout,&
                RAPID_routing_struc(n)%Qout,(/1,1/),(/RAPID_routing_struc(n)%n_riv_bas,1/)),&
                'Error in nf90_get_var in RAPID_routing_readrst')

           ! Close the file
           call LIS_verify(nf90_close(ftn),'Error in nf90_close in RAPID_routing_readrst') 
        endif !if(LIS_masterproc) then
     endif !f(read_restart) then
  enddo 

end subroutine RAPID_routing_readrst
