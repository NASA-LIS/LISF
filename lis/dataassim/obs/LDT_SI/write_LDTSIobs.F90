!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"

! Write the transformed LDTSI observations to a file.
! TODO: Wrap this in a module
subroutine write_LDTSIobs(n, k, OBS_State)

   ! Imports
   use ESMF
   use LIS_coreMod, only: LIS_masterproc
   use LIS_DAobservationsMod, only: LIS_writevar_gridded_obs
   use LIS_fileIOMod, only: LIS_create_output_directory
   use LIS_logMod, only: LIS_verify, LIS_getNextUnitNumber, &
        LIS_releaseUnitNumber

   ! Defaults
   implicit none

   ! Arguments
   integer, intent(in) :: n
   integer, intent(in) :: k
   type(ESMF_State), intent(inout) :: OBS_State 

   ! NOTE:  intent(inout) above for OBS_State is to maintain compatibility
   ! with ESMF 5.2.0rp3

   ! Local variables
   logical                  :: data_update
   integer                  :: ftn
   character(100)           :: obsname
   type(ESMF_Field)         :: snowField
   real, pointer            :: snowobs(:)
   integer                  :: status

   call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
        data_update, rc=status)
   call LIS_verify(status)

   if (data_update) then
      
      call ESMF_StateGet(OBS_State, "Observation01", snowField, &
           rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(snowField, localDE=0, farrayPtr=snowobs, rc=status)
      call LIS_verify(status)
      
      if (LIS_masterproc) then 
         ftn = LIS_getNextUnitNumber()
         call LDTSI_obsname(n, k, obsname)        
         
         call LIS_create_output_directory('DAOBS')
         open(ftn, file=trim(obsname), form='unformatted')
      endif

      call LIS_writevar_gridded_obs(ftn, n, k, snowobs)

      if (LIS_masterproc) then 
         call LIS_releaseUnitNumber(ftn)
      end if
      
   end if
end subroutine write_LDTSIobs

! Construct output file name
subroutine LDTSI_obsname(n, k, obsname)

   ! Imports
   use LIS_coreMod, only: LIS_rc

   ! Defaults
   implicit none

   ! Arguments
   integer, intent(in) :: n
   integer, intent(in) :: k
   character(len=*), intent(out) :: obsname

   ! Local variables
   character(len=10) :: cda
   character(len=12) :: cdate
   character(len=12) :: cdate1

   write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LIS_rc%yr, LIS_rc%mo, &
        LIS_rc%da, LIS_rc%hr,LIS_rc%mn
   
   write(unit=cda, fmt='(a2,i2.2)') '.a',k
   write(unit=cdate, fmt='(a2,i2.2)') '.d',n

   obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//&
        '/LISDAOBS_'//cdate1// &
        trim(cda)//trim(cdate)//'.1gs4r'

end subroutine LDTSI_obsname
