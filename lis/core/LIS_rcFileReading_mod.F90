#include "LIS_misc.h"
!------------------------------------------------------------------------------
! NASA GSFC 
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: LIS_rcFileReading_mod
!
! !INTERFACE:
!
module LIS_rcFileReading_mod
!
! !USES:
   use ESMF
   use LIS_logMod
!
   implicit none
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public  :: LIS_read_rc_variables

   interface LIS_read_rc_variables
      module procedure read_rc_var_int
      module procedure read_rc_var_list_int
      module procedure read_rc_var_real
      module procedure read_rc_var_list_real
      module procedure read_rc_var_bool
      module procedure read_rc_var_list_bool
      module procedure read_rc_var_char
      module procedure read_rc_var_list_char
   end interface
!
! !DESCRIPTION:
! Basic utility routines for reading resource file using ESMF calls.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!EOP
!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_int 
!
! !INTERFACE:
!
   subroutine read_rc_var_int(config, var_name, label, default, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*),  intent(in) :: label
      integer, optional, intent(in) :: default
!
! !OUTPUT PARAMETERS:
      integer,           intent(out) :: var_name
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in an integer variable from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_int"
!EOP
!------------------------------------------------------------------------------
!BOC
      if (present(default)) then
         call ESMF_ConfigGetAttribute(config, var_name, &
            label=label, default=default, rc=STATUS)
      else
         call ESMF_ConfigGetAttribute(config, var_name, label=label, rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      endif

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_int
!EOC
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_list_int 
!
! !INTERFACE:
!
   subroutine read_rc_var_list_int(config, var_name, label, num, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: label
      integer,          intent(in) :: num
!
! !OUTPUT PARAMETERS:
      integer,           intent(out) :: var_name(num)
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a list of integer variables from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, i
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_list_int"
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Check if variable is set in the file
      call ESMF_ConfigFindLabel(config, label, rc=STATUS)

      do i = 1, num
         call ESMF_ConfigGetAttribute(config, var_name(i), rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      enddo

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_list_int
!EOC
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_real 
!
! !INTERFACE:
!
   subroutine read_rc_var_real(config, var_name, label, default, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: label
      real,   optional, intent(in) :: default
!
! !OUTPUT PARAMETERS:
      real,              intent(out) :: var_name
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a real variable from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_real"
!EOP
!------------------------------------------------------------------------------
!BOC
      if (present(default)) then
         call ESMF_ConfigGetAttribute(config, var_name, &
            label=label, default=default, rc=STATUS)
      else
         call ESMF_ConfigGetAttribute(config, var_name, label=label, rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      endif

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_real
!EOC
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_list_real 
!
! !INTERFACE:
!
   subroutine read_rc_var_list_real(config, var_name, label, num, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: label
      real,             intent(in) :: num
!
! !OUTPUT PARAMETERS:
      real,              intent(out) :: var_name(num)
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a list of real variables from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, i
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_list_real"
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Check if variable is set in the file
      call ESMF_ConfigFindLabel(config, label, rc=STATUS)

      do i = 1, num
         call ESMF_ConfigGetAttribute(config, var_name(i), rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      enddo

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_list_real
!EOC
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_char 
!
! !INTERFACE:
!
   subroutine read_rc_var_char(config, var_name, label, default, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: label
      character(len=*), optional,  intent(in) :: default
!
! !OUTPUT PARAMETERS:
      character(len=*),  intent(out) :: var_name
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a character-based variable from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_char"
!EOP
!------------------------------------------------------------------------------
!BOC
      if (present(default)) then
         call ESMF_ConfigGetAttribute(config, var_name, &
            label=label, default=default, rc=STATUS)
      else
         call ESMF_ConfigGetAttribute(config, var_name, label=label, rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      endif

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_char
!EOC
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_list_char 
!
! !INTERFACE:
!
   subroutine read_rc_var_list_char(config, var_name, label, num, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: label
      real,             intent(in) :: num
!
! !OUTPUT PARAMETERS:
      character(len=*),  intent(out) :: var_name(num)
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a list of character-based variables from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, i
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_list_char"
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Check if variable is set in the file
      call ESMF_ConfigFindLabel(config, label, rc=STATUS)

      do i = 1, num
         call ESMF_ConfigGetAttribute(config, var_name(i), rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      enddo

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_list_char
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_rc_var_bool 
!
! !INTERFACE:
!
   subroutine read_rc_var_bool(config, var_name, label, default, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*),  intent(in) :: label
      logical, optional, intent(in) :: default
!
! !OUTPUT PARAMETERS:
      logical,           intent(out) :: var_name
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a boolean variable from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_bool"
!EOP
!------------------------------------------------------------------------------
!BOC
      if (present(default)) then
         call ESMF_ConfigGetAttribute(config, var_name, &
            label=label, default=default, rc=STATUS)
      else
         call ESMF_ConfigGetAttribute(config, var_name, label=label, rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      endif

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_bool
!EOC
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: read_rc_var_list_bool 
!
! !INTERFACE:
!
   subroutine read_rc_var_list_bool(config, var_name, label, num, rc)
!
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: label
      integer,          intent(in) :: num
!
! !OUTPUT PARAMETERS:
      logical,           intent(out) :: var_name(num)
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a list of boolean variables from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, i
      character(len=ESMF_MAXSTR), parameter :: IAm = "read_rc_var_list_bool"
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Check if variable is set in the file
      call ESMF_ConfigFindLabel(config, label, rc=STATUS)

      do i = 1, num
         call ESMF_ConfigGetAttribute(config, var_name(i), rc=STATUS)
         call LIS_verify(STATUS, TRIM(label)//" not defined")
      enddo

      if (present(rc)) rc = STATUS

      return

   end subroutine read_rc_var_list_bool
!EOC
!------------------------------------------------------------------------------
end module LIS_rcFileReading_mod
