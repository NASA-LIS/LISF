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
! !ROUTINE: HYMAP2_setPertStates
! \label{HYMAP2_setPertStates}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP2_setPertStates(n,Nstate,pert_State,progpert)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use HYMAP2_routingMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: Nstate
  type(ESMF_State)       :: pert_State
  real                   :: progpert(Nstate,&
       LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!   This routine sets the perturbation object with values
!   from the corresponding ESMF data structure. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]           index of the nest \newline
!  \item[Nstate]      number of state variables
!  \item[pert\_State] ESMF state perturbation object
!  \item[progpert]    prognostic state perturbation state object used
!                     in the perturbation algorithm
!  \end{description}
!EOP
  
  integer                       :: i,t,m,kk,col,row
  character*100, pointer        :: pertobjs(:)
  type(ESMF_Field), allocatable :: pertField(:)
  real,             pointer     :: pertdata1d(:)
  integer                       :: status

  allocate(pertobjs(Nstate))
  allocate(pertField(Nstate))

  call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
  call LIS_verify(status)

  do i=1,Nstate
     call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet failed in gmaoperturb')
     call ESMF_FieldGet(pertField(i),localDE=0,&
          farrayPtr=pertdata1d,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed in gmaoperturb')

     do t=1,HYMAP2_routing_struc(n)%nseqall
        do m=1,LIS_rc%nensem(n)
           kk = (t-1)*LIS_rc%nensem(n)+m
           col = HYMAP2_routing_struc(n)%seqx(t)
           row = HYMAP2_routing_struc(n)%seqy(t)
           pertdata1d(kk) = progpert(i,col,row,m)
        enddo
     enddo
  enddo
end subroutine HYMAP2_setPertStates
