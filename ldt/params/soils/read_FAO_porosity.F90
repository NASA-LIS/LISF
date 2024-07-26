!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_FAO_porosity
! \label{read_FAO_porosity}
!
! !REVISION HISTORY:
!  27 Sept 2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_FAO_porosity(n, array)
! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod, only      : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n            !nest index
  real,    intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),3)

! !DESCRIPTION:
!  This subroutine retrieves FAO porosity fraction data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved porosity fraction
!  \end{description}
!EOP

  integer       :: ftn
  character(len=LDT_CONST_PATH_LEN) :: filename
  logical       :: file_exists
  real          :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
! ________________________________

  array = LDT_rc%udef

  filename = trim(LDT_rc%pofile(n))//trim('.L1.1gd4r')
  write(LDT_logunit,*) "[INFO] Reading FAO porosity file: ",trim(filename)

  inquire(file=trim(filename), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Porosity map ",trim(filename)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  ftn = LDT_getNextUnitNumber()
  open(ftn,file=filename, form='unformatted',status='old',&
       access='direct',recl=4)

  temp = LDT_rc%udef
  call readLISdata(n, ftn, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                   LDT_rc%soil_gridDesc(n,:), 1, temp )  
  array(:,:,1) = temp(:,:,1)

  call LDT_releaseUnitNumber(ftn)

  filename = trim(LDT_rc%pofile(n))//trim('.L2.1gd4r')
  inquire(file=trim(filename), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'Porosity map ',trim(filename),' not found'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif

  write(LDT_logunit,*) '[INFO] Reading FAO porosity file',trim(filename)
  ftn = LDT_getNextUnitNumber()
  open(ftn,file=filename, form='unformatted',status='old',&
       access='direct',recl=4)

  temp = LDT_rc%udef
  call readLISdata(n, ftn, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                  LDT_rc%soil_gridDesc(n,:), 1, temp )
  array(:,:,2) = temp(:,:,1)

  call LDT_releaseUnitNumber(ftn)

  filename = trim(LDT_rc%pofile(n))//trim('.L3.1gd4r')
  inquire(file=trim(filename), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'Porosity map ',trim(filename),' not found'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif

  write(LDT_logunit,*) "[INFO] Reading FAO porosity file: ",trim(filename)
  ftn = LDT_getNextUnitNumber()
  open(ftn,file=filename,form='unformatted',status='old',&
       access='direct',recl=4)

  temp = LDT_rc%udef
  call readLISdata(n, ftn, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                  LDT_rc%soil_gridDesc(n,:), 1, temp)
  array(:,:,3) = temp(:,:,1)

  call LDT_releaseUnitNumber(ftn)

  write(LDT_logunit,*) "[INFO] Read porosity files"

end subroutine read_FAO_porosity
