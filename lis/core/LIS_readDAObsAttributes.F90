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
! !ROUTINE: LIS_readDAObsAttributes
! \label{LIS_readDAObsAttributes}
!
! !REVISION HISTORY:
!  14 Oct 2006: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine LIS_readDAObsAttributes(k, vname,varmin,varmax)
! !USES: 
    use LIS_coreMod,   only : LIS_rc
    use LIS_logMod,    only : LIS_logunit, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber

    implicit none
! !ARGUMENTS: 
    integer                    :: k
    character(len=*), allocatable  :: vname(:)
    real        , allocatable      :: varmin(:)
    real        , allocatable      :: varmax(:)
    integer                    :: ftn 
! 
! !DESCRIPTION: 
!  This routine reads the observation attributes for each
!  observation type. The format of the attributes file is: 
!
!  Variable name 
!  Min Value  Max Value
!
!  The arguments are: 
!  \begin{description}
!   \item[n]           index of nest
!   \item[vname]       name of the observation variables
!   \item[varmin]      minimum value of the variable
!   \item[varmax]      maximum value of the variable
!  \end{description}
!EOP
    integer            :: i 


    write(LIS_logunit,*) 'Opening attributes for observations ',&
         trim(LIS_rc%obsattribfile(k))
    ftn = LIS_getNextUnitNumber()
    open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
    read(ftn,*)
    read(ftn,*) LIS_rc%nobtypes(k)
    read(ftn,*)
    
    allocate(vname(LIS_rc%nobtypes(k)))
    allocate(varmax(LIS_rc%nobtypes(k)))
    allocate(varmin(LIS_rc%nobtypes(k)))
        
    do i=1,LIS_rc%nobtypes(k)
       read(ftn,fmt='(a40)') vname(i)
       read(ftn,*) varmin(i),varmax(i)
       write(LIS_logunit,*) vname(i),varmin(i),varmax(i)
    enddo
    call LIS_releaseUnitNumber(ftn)
  end subroutine LIS_readDAObsAttributes
