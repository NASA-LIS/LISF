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
! !ROUTINE: read_HYMAP_fld_height
! \label{read_HYMAP_fld_height}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_HYMAP_fld_height(n, array)
! !USES:
  use ESMF
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata 
  use LDT_paramDataMod
  use HYMAP_parmsMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer,   intent(in) :: n
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                HYMAP_struc(n)%hymap_fld_height%vlevels)

! !DESCRIPTION:
!  This subroutine retrieves the bottom temperature climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r,k
  logical :: file_exists
  integer      :: nc_dom, nr_dom
  integer      :: line1, line2, line
  real         :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real         :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))

  ftn = LDT_getNextUnitNumber()

  inquire(file=trim(HYMAP_struc(n)%fldheightfile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) '[ERR] HYMAP floodplain height map, ',trim(HYMAP_struc(n)%fldheightfile),&
                          ', not found.'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif

  open(ftn, file=trim(HYMAP_struc(n)%fldheightfile), access='direct',&
       status='old', form="unformatted", convert="big_endian", recl=4)
  
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
             rlat(c,r),rlon(c,r))
     enddo
  enddo

  ! Yeosang Yoon: fix typo
  nc_dom = nint((HYMAP_struc(n)%hymapparms_gridDesc(8)-&
       HYMAP_struc(n)%hymapparms_gridDesc(5))/(HYMAP_struc(n)%hymapparms_gridDesc(9)))+1
  nr_dom = nint((HYMAP_struc(n)%hymapparms_gridDesc(7)-&
       HYMAP_struc(n)%hymapparms_gridDesc(4))/(HYMAP_struc(n)%hymapparms_gridDesc(10)))+1

  do k=1,HYMAP_struc(n)%hymap_fld_height%vlevels
!     call readLISdata(n, ftn, LDT_rc%hymap_proj, &
!          LDT_rc%hymap_gridtransform, &
!          HYMAP_struc(n)%hymapparms_gridDesc(:), 1, array(:,:,k))
     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           line1 = nint((rlat(c,r)-HYMAP_struc(n)%hymapparms_gridDesc(4))/&
                HYMAP_struc(n)%hymapparms_gridDesc(10))+1
           line2 = nint((rlon(c,r)-HYMAP_struc(n)%hymapparms_gridDesc(5))/&
                HYMAP_struc(n)%hymapparms_gridDesc(9))+1
           line = (k-1)*nc_dom*nr_dom + (line1-1)*nc_dom + line2
           read(ftn,rec=line) array(c,r,k)
        enddo
     enddo

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(array(c,r,k).lt.0) then
              array(c,r,k) = LDT_rc%udef
           endif
        enddo
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

end subroutine read_HYMAP_fld_height
