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
! !ROUTINE: NoahMP401_getsmpred
! \label{NoahMP401_getsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari; Modified for NoahMP401 
!
! !INTERFACE:
subroutine NoahMP401_getsmpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use NoahMP401_lsmMod
  use NoahMP401_dasoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the Soil moisture obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6), sm_out
  character*50           :: units_tp(6)
  real                   :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index))


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       obs_pred)
  
#if 0 
  if(LIS_rc%useANNinDA(n).ne.1) then      
!rescaling to relative wetness space. 
     if(LIS_rc%dascaloption(k).eq."Linear scaling") then 

        if(noahmp401_dasm_struc(n)%ntimes.gt.1) then 
           kk = LIS_rc%mo
        else
           kk = 1
        endif

        do i=1,LIS_rc%obs_ngrid(k)
           do m=1,LIS_rc%nensem(n)
              obs_tmp = (obs_pred(i,m) - noahmp401_dasm_struc(n)%model_xrange(i,kk,1))/&
                   (noahmp401_dasm_struc(n)%model_xrange(i,kk,noahmp401_dasm_struc(n)%nbins) - & 
                   noahmp401_dasm_struc(n)%model_xrange(i,kk,1))
              obs_pred(i,m) = obs_tmp
           enddo
        enddo
     endif

  else
     obs_pred = 0.0
     count1 = 0 
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)
           t = i+m-1
           gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
           
           inputs_tp(1) = noahmp401_struc(n)%noahmp401(t)%smc(1)
           inputs_tp(2) = noahmp401_struc(n)%noahmp401(t)%prcp*(1-noahmp401_struc(n)%noahmp401(t)%fpice) !Noah33 rainf
           inputs_tp(3) = noahmp401_struc(n)%noahmp401(t)%prcp*noahmp401_struc(n)%noahmp401(t)%fpice     !Noah33 snowf    
           !MN: NOTE:  noahmp401 --> prcp (total precip), fpice (snow fraction in precipitation [-])
           inputs_tp(4) = noahmp401_struc(n)%noahmp401(t)%fveg !Noah33 shdfac
           inputs_tp(5) = noahmp401_struc(n)%noahmp401(t)%sstc(noahmp401_struc(n)%nsnow+1)  !Noah33 stc(1)
           inputs_tp(6) = noahmp401_struc(n)%noahmp401(t)%sneqv*LIS_CONST_RHOFW

           units_tp(1) = "m^3 m-3"
           units_tp(2) = "kg m-2 s-1"
           units_tp(3) = "kg m-2 s-1"
           units_tp(4) = "1"
           units_tp(5) = "K"
           units_tp(6) = "kg m-2"

           call LIS_forwardEstimate_with_ANN(n, gid, inputs_tp, &
                units_tp, sm_out)
           obs_pred(gid,m) = obs_pred(gid,m) + sm_out
           count1(gid,m) = count1(gid,m) + 1
        enddo
     enddo
     
     do i=1,LIS_rc%ngrid(n)
        do m=1,LIS_rc%nensem(n)
           obs_pred(i,m) = obs_pred(i,m)/(count1(i,m))
        enddo
     enddo
  endif
#endif
end subroutine NoahMP401_getsmpred

