!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_ceop
! \label{read_ceop}
! 
! !REVISION HISTORY:
!
!  08 Dec 2002: Sujay Kumar; Initial Specification 
!
! !INTERFACE:      
subroutine read_ceop(n,findex,order) 
! !USES:
  use LIS_logMod, only        : LIS_logunit
  use LIS_coreMod, only     : LIS_rc,LIS_domain
  use LIS_metforcingMod,  only : LIS_forc
  use ceop_forcingMod, only : ceop_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: n
  integer, intent(in)           :: findex
  integer, intent(in)           :: order
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  the correct CEOP station data (ASCII), transforms into LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    instance, order=2, read the next instance)
!  \item[n]
!    index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[readbon](\ref{readbon}) \newline
!    reads the Bondville, IL data
!  \item[readfpk](\ref{readfpk}) \newline
!    reads the Fort Peck, MT data 
!  \item[readsgp](\ref{readsgp}) \newline
!    reads the SGP data 
!  \item[normalize\_stnwts](\ref{normalize_stnwts}) \newline
!    renormalizes the station weights accounting for
!    missing data
!  \item[interp\_stndata](\ref{interp_stndata}) \newline
!    spatially interpolates the station data onto the 
!    LIS grid.
!  \end{description}
!EOP  

  integer :: c,r,count1,f
  real :: tair(ceop_struc(n)%nstns)
  real :: qair(ceop_struc(n)%nstns)
  real :: u(ceop_struc(n)%nstns),v(ceop_struc(n)%nstns)
  real :: pcp(ceop_struc(n)%nstns)
  real :: psurf(ceop_struc(n)%nstns)
  real :: swdown(ceop_struc(n)%nstns)
  real :: lwdown(ceop_struc(n)%nstns)
  real :: varfield(8,LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer :: errorcode 

  varfield = 0.0
  varfield1 = 0.0
  if(ceop_struc(n)%location.eq.1) then 
!     call readbon(n,tair,qair,swdown,lwdown,u,v,psurf,pcp)
!  elseif(ceop_struc(n)%location.eq.2) then 
!     call readfpk(n,tair,qair,swdown,lwdown,u,v,psurf,pcp)
  elseif(ceop_struc(n)%location.eq.3) then 
     call readsgp(n,tair,qair,swdown,lwdown,u,v,psurf,pcp,errorcode)
  endif
  if(errorcode.eq.0) then 
     do c=1,ceop_struc(n)%nstns
        if(psurf(c).ne.-999.99) then 
           psurf(c) = psurf(c)*100
        endif
        if(tair(c).ne.-999.99) then 
           tair(c) = tair(c)+273.16
        endif
        if(qair(c).ne.-999.99) then 
           qair(c) = qair(c)/1000.0
        endif
     enddo
     call normalize_stnwts(psurf,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)
     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,psurf,varfield(7,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     
     call normalize_stnwts(tair,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)
     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,tair,varfield(1,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     
     call normalize_stnwts(qair,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)
     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,qair,varfield(2,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     
     call normalize_stnwts(u,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)
     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,u,varfield(5,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     
     call normalize_stnwts(v,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)
     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,v,varfield(6,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)

     call normalize_stnwts(swdown,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)

     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,swdown,varfield(3,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     
     call normalize_stnwts(lwdown,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)

     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,lwdown,varfield(4,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     
     call normalize_stnwts(pcp,ceop_struc(n)%nstns,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%undef,ceop_struc(n)%stnwt)
     call interp_stndata(ceop_struc(n)%stnwt,ceop_struc(n)%undef,pcp,varfield(8,:),&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns)
     do f=1,8
        count1 = 0 
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              varfield1(c,r) = varfield(f,c+count1)
           enddo
           count1 = count1 + LIS_rc%lnc(n)
        enddo
        
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then
                    ceop_struc(n)%metdata1(f,LIS_domain(n)%gindex(c,r)) =&
                         varfield1(c,r)
                 elseif(order.eq.2) then 
                    ceop_struc(n)%metdata2(f,LIS_domain(n)%gindex(c,r)) =&
                         varfield1(c,r)
                 endif
              endif
           enddo
        enddo
     enddo
  else
     if(order.eq.1) then 
        ceop_struc(n)%metdata1 = -999.99
     elseif(order.eq.2) then 
        ceop_struc(n)%metdata2 = -999.99
     endif

  endif
333 format(I4,4(1X,I2),1X,I4,4(1X,I2),43X,F11.5,F12.5,12F8.2,&
          11F9.2)    

end subroutine read_ceop
