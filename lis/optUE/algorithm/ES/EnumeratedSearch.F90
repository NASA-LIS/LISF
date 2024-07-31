!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module EnumeratedSearch
!BOP
!
! !MODULE: EnumeratedSearch
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!  19 Aug 2009; Sujay Kumar; Initial Specification
!
  use ESMF
  use ES_varctl

  implicit none
  
  PRIVATE
  
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ESOpt_init
  public :: ESOpt_setup
  public :: ESOpt_run
  public :: ESOpt_checkConvergence
  public :: ES_getdecSpaceValues
  public :: ES_setdecSpaceValues
  public :: ES_getNparam
!EOP

  type esstruc
     real, allocatable :: parent(:)
     real, allocatable :: best(:)
     real          :: best_obj
  end type esstruc

  type(esctl)              :: es_ctl
  type(esstruc),   allocatable :: es_struc(:)

contains 

!BOP
! !ROUTINE: ESOpt_init
! \label{ESOpt_init}
! 
! !INTERFACE: 
    subroutine ESOpt_init()
! !USES:   
      use LIS_coreMod,         only : LIS_rc, LIS_config, LIS_vecTile
      use LIS_optUEMod, only : LIS_decisionSpace
      use LIS_logMod,          only : LIS_logunit, LIS_getNextUnitNumber, &
           LIS_releaseUnitNumber, LIS_endrun, LIS_verify

! !DESCRIPTION: 
!   This routine performs the initialization steps for the Enumerated Search.  
!
!EOP      
      implicit none 

      character*100,  allocatable     :: vname(:)
      integer                     :: i,t
      integer                     :: ftn
      integer                     :: n 

      integer                     :: status
      type(ESMF_Field)            :: varField
      type(ESMF_Field)            :: fitField
      type(ESMF_ArraySpec)        :: arrspec1

      n = 1
      call ESMF_ConfigGetAttribute(LIS_config,es_ctl%decspaceAttribsFile,&
           label="ES decision space attributes file:",rc=status)
      call LIS_verify(status,'ES decision space attributes file: not defined')
      

      ftn = LIS_getNextUnitNumber()
      write(LIS_logunit,*) 'Reading decision space attributes ...', &
           es_ctl%decspaceAttribsFile
      open(ftn,file=trim(es_ctl%decspaceAttribsFile),status='old')
      read(ftn,*)
      read(ftn,*) es_ctl%nparam

      allocate(es_ctl%parmax(es_ctl%nparam))
      allocate(es_ctl%parmin(es_ctl%nparam))
      allocate(es_ctl%pardel(es_ctl%nparam))
      allocate(es_ctl%npossbl(es_ctl%nparam))
      allocate(es_ctl%parm_count(es_ctl%nparam))
     
      es_ctl%parm_count = 1 

      allocate(vname(es_ctl%nparam))

      read(ftn,*) 
      do i=1,es_ctl%nparam
         read(ftn,fmt='(a100)') vname(i)
         write(LIS_logunit,*) 'vname ',vname(i)
      enddo
      read(ftn,*) 
      read(ftn,*) (es_ctl%parmax(i),i=1,es_ctl%nparam)
      read(ftn,*) 
      read(ftn,*) (es_ctl%parmin(i),i=1,es_ctl%nparam)
      read(ftn,*) 
      read(ftn,*) (es_ctl%pardel(i),i=1,es_ctl%nparam)
      
      es_ctl%ntotal = 1
      do i=1,es_ctl%nparam
         es_ctl%npossbl(i) = nint((es_ctl%parmax(i)-es_ctl%parmin(i))/&
              es_ctl%pardel(i)) + 1
         es_ctl%ntotal = es_ctl%ntotal*es_ctl%npossbl(i)
      enddo

      write(LIS_logunit,*) 'Total number of enumerations: ',es_ctl%ntotal

      call LIS_releaseUnitNumber(ftn)
      write(LIS_logunit,*) 'Finished reading decision space attributes ..'

      es_ctl%genNo = 0 
      
      allocate(es_struc(LIS_rc%ngrid(n)))

      call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
           rc=status)
      call LIS_verify(status)
      do i=1,es_ctl%nparam
         varField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
              name=trim(vname(i)),rc=status)
         call LIS_verify(status)
         call ESMF_StateAdd(LIS_decisionSpace, (/varField/),rc=status)
         call LIS_verify(status)
      enddo

      deallocate(vname)

      do t=1,LIS_rc%ngrid(n)
! population (for each tile point)
         allocate(es_struc(t)%parent(es_ctl%nparam))
         allocate(es_struc(t)%parent(es_ctl%nparam))
         allocate(es_struc(t)%best(es_ctl%nparam))
      
         es_struc(t)%parent(:) = es_ctl%parmin(:)
         es_struc(t)%best_obj = 1E8
      enddo


    end subroutine ESOpt_init


!BOP
! !ROUTINE: ESOpt_setup
! \label{ESOpt_setup}
!
! !INTERFACE: ESOpt_setup
    subroutine ESOpt_setup()
! !USES: 
!EOP
    end subroutine ESOpt_setup

!BOP
! 
! !ROUTINE: ESOpt_checkConvergence
! \label{ESOpt_checkConvergence}
! 
! !INTERFACE: 
    subroutine ESOpt_checkConvergence(check)
! !ARGUMENTS: 
      logical, intent(INOUT) :: check
!EOP
      if(es_ctl%genNo.ge.es_ctl%ntotal) then 
         check = .true.
      else
         check = .false.
      endif

    end subroutine ESOpt_checkConvergence

!BOP
! !ROUTINE: ESOpt_run
! \label{ESOpt_run}
! 
! !INTERFACE: 
    subroutine ESOpt_run()
! !USES: 
      use LIS_coreMod, only : LIS_rc
      use LIS_optUEMod, only : LIS_ObjectiveFunc, LIS_feasibleSpace
      use LIS_logMod,  only : LIS_verify

      implicit none
!EOP      
      integer            :: i, t, n
      type(ESMF_Field)   :: objField
      type(ESMF_Field)   :: feasField
      real, pointer      :: objValue(:)
      integer, pointer   :: mod_flag(:)
      integer            :: status

      n = 1
     
      es_ctl%genNo = es_ctl%genNo+1

      call evaluateobjfunction(LIS_rc%optuetype)

      call ESMF_StateGet(LIS_ObjectiveFunc,"Min Criteria Value",objField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(objField,localDE=0, farrayPtr=objValue,rc=status)
      call LIS_verify(status)
      
      call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
      call LIS_verify(status)

      n = 1
      do t=1,LIS_rc%ngrid(n)
!penalize infeasible solutions
          if(mod_flag(t).eq.1) then 
             objValue(t) = -1.0
          endif
          write(122,*) (es_struc(t)%parent(i), i=1,es_ctl%nparam), objValue(t)
       end do
    
!adjust decision space 
!for now this is structured for a two parameter problem. 

      if(mod(es_ctl%parm_count(1), es_ctl%npossbl(1)).eq.0) then 
         es_ctl%parm_count(1) =1
         es_ctl%parm_count(2) = es_ctl%parm_count(2) + 1
      else
         es_ctl%parm_count(1) = es_ctl%parm_count(1) + 1
      endif

      do i=1,es_ctl%nparam
         do t=1,LIS_rc%ngrid(n)
            es_struc(t)%parent(i) = es_ctl%parmin(i)+&
                 (es_ctl%parm_count(i)-1) * es_ctl%pardel(i)
         enddo
      enddo
    
!      call setoptuetypedecspace(LIS_rc%optuetype)

!update the best soln
      do t=1,LIS_rc%ngrid(n)
         if(objValue(t).lt.es_struc(t)%best_obj) then 
            es_struc(t)%best = es_struc(t)%parent
            es_struc(t)%best_obj = objValue(t)
         endif
         write(123,*) (es_struc(t)%best(i), i=1,es_ctl%nparam), es_struc(t)%best_obj
      enddo

    end subroutine ESOpt_run
!BOP
! !ROUTINE: ES_getdecSpaceValues
! \label{ES_getdecSpaceValues}
! 
! !INTERFACE: 
    subroutine ES_getdecSpaceValues(n, decvals)
! !USES: 
      use LIS_coreMod,  only : LIS_rc
!EOP
      implicit none

      integer            :: n
      real               :: decvals(es_ctl%nparam, LIS_rc%ntiles(n))
      integer            :: i,t

      do i=1, es_ctl%nparam
         do t=1, LIS_rc%ngrid(n)
            decvals(i,t) = es_struc(t)%parent(i)
         enddo
      enddo

    end subroutine ES_getdecSpaceValues
!BOP
! 
! !ROUTINE: ES_setdecSpaceValues
! 
! !INTERFACE: 
    subroutine ES_setdecSpaceValues(n, decvals)
! !USES: 
      use LIS_coreMod,  only : LIS_rc
!EOP
      implicit none

      integer            :: n
      real               :: decvals(es_ctl%nparam, LIS_rc%ntiles(n))

      integer            :: i,t
      do i=1, es_ctl%nparam
         do t=1, LIS_rc%ngrid(n)
            es_struc(t)%parent(i) = decvals(i,t)
         enddo
      enddo

    end subroutine ES_setdecSpaceValues

   subroutine ES_getNparam(nparam)
      
      integer   :: nparam
      
      nparam = es_ctl%nparam

    end subroutine ES_getNparam

end module EnumeratedSearch
