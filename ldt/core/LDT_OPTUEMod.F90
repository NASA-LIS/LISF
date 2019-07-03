!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_OPTUEMod
!BOP
!
! !MODULE: LDT_OPTUEMod
!
! !DESCRIPTION:
!  The code in this file implements routines to process the output from 
!  the LIS OPT/UE system and generate spatially distributed parameters.
! 
!
! !REVISION HISTORY:
!
!  25 Jun 2019: Sujay Kumar; Initial implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_logMod
  use LDT_paramDataMod

  implicit none

!- LSM-specific parameters:
  type, public :: optueparam_type_dec
     integer                             :: selectOpt
     integer                             :: nparam
     type(LDT_paramEntry), allocatable   :: param(:)
     
  end type optueparam_type_dec
  
  type(optueparam_type_dec), allocatable :: LDT_OPTUEparam_struc(:)

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_optue_init  ! initializes data structures and memory
  public :: LDT_optue_writeHeader
  public :: LDT_optue_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

!EOP

contains

!BOP
! 
! !ROUTINE: LDT_optue_init
! \label{LDT_optue_init}
! 
! !INTERFACE:
  subroutine LDT_optue_init()
! !USES:

! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! optue datasets
!
!EOP
    implicit none
    integer :: n, i,t,c,r
    integer :: rc
    integer :: ftn
    character*100 :: filename
    real, allocatable :: fitness(:,:),avgfitness(:,:)
    logical           :: file_exists
! ______________________________________________________________
    
    allocate(LDT_OPTUEparam_struc(LDT_rc%nnest))
    LDT_OPTUEparam_struc(:)%selectOpt = 1
    write(LDT_logunit,*)" - - - - - - - - - OPT/UE Parameters - - - - - - - - - - - -"

    call ESMF_ConfigFindLabel(LDT_config,"LIS OPT/UE output file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,filename,rc=rc)
       call LDT_verify(rc,"LIS OPT/UE output file: not defined")
    enddo

    do n=1,LDT_rc%nnest
       inquire(file=trim(filename),exist=file_exists) 

       if(file_exists) then 
          ftn = LDT_getNextUnitNumber()
          
          open (ftn,file=filename,form='unformatted')
          read(ftn) LDT_OPTUEparam_struc(n)%nparam
          
          allocate(LDT_OPTUEparam_struc(n)%param(&
               LDT_OPTUEparam_struc(n)%nparam))

          do t=1,LDT_OPTUEparam_struc(n)%nparam
             allocate(LDT_OPTUEparam_struc(n)%param(t)%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),1))
          enddo
          
          allocate(fitness(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(avgfitness(LDT_rc%lnc(n),LDT_rc%lnr(n)))

          do t=1,LDT_OPTUEparam_struc(n)%nparam+2
             if(t==LDT_OPTUEparam_struc(n)%nparam+1) then 
                read(ftn) fitness
             elseif(t==LDT_OPTUEparam_struc(n)%nparam+2) then 
                read(ftn) avgfitness
             else
                read(ftn) LDT_OPTUEparam_struc(n)%param(t)%short_name
                write(LDT_logunit,*) '[INFO] Reading OPT/UE parameter ',&
                     trim(LDT_OPTUEparam_struc(n)%param(t)%short_name)
                read(ftn) LDT_OPTUEparam_struc(n)%param(t)%value(:,:,1)
                LDT_OPTUEparam_struc(n)%param(t)%selectOpt = 1
                LDT_OPTUEparam_struc(n)%param(t)%vlevels = 1
                LDT_OPTUEparam_struc(n)%param(t)%standard_name = ''
                LDT_OPTUEparam_struc(n)%param(t)%units = ''
                LDT_OPTUEparam_struc(n)%param(t)%valid_min = 0
                LDT_OPTUEparam_struc(n)%param(t)%valid_max = 0
                LDT_OPTUEparam_struc(n)%param(t)%num_bins = 1
             endif
          enddo
          call LDT_releaseUnitNumber(ftn)
       else
          write(LDT_logunit,*) '[ERR] LIS OPT/UE file ',trim(filename), &
               'not found'
          call LDT_endrun()
       endif

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(fitness(c,r).eq.-1E20) then 
                do t=1,LDT_OPTUEparam_struc(n)%nparam
                   LDT_OPTUEparam_struc(n)%param(t)%value(c,r,1) = LDT_rc%udef
                enddo
             endif
          enddo
       enddo
       
    enddo

  end subroutine LDT_optue_init



  subroutine LDT_optue_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer   :: n 
    integer   :: ftn
    integer   :: dimID(3)
    integer   :: tdimID(3)

    integer   :: t

    if(LDT_OPTUEparam_struc(n)%selectOpt==1) then 
       tdimID(1) = dimID(1)
       tdimID(2) = dimID(2)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

       do t=1,LDT_OPTUEparam_struc(n)%nparam
          call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
               LDT_OPTUEparam_struc(n)%param(t))
       enddo
#endif
    endif
  end subroutine LDT_optue_writeHeader

  subroutine LDT_optue_writeData(n,ftn)

    use LDT_coreMod, only : LDT_rc

    integer  :: n 
    integer  :: ftn
   
    integer   :: t

    if(LDT_OPTUEparam_struc(n)%selectOpt==1) then 
       do t=1,LDT_OPTUEparam_struc(n)%nparam
          call LDT_writeNETCDFdata(n,ftn,LDT_OPTUEparam_struc(n)%param(t))
       enddo
    endif
  end subroutine LDT_optue_writeData

end module LDT_OPTUEMod
