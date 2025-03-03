!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      subroutine agrlwdn( sfctmp, e, iclamt, rldown )

!RENAMED IN LDAS
!      subroutine longwv( sfctmp, e, iclamt, rldown )  

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  
!     name:  longwave radiation routine    
!     ====                      
!  
!     purpose:
!     ======= 
!     to compute the net downward longwave radiation at the
!     earth's surface.
!  
!     method:  
!     ====== 
!     - calculate the emissivity of the clear sky.
!     - calculate downwelling longwave radiation of the clear sky.
!     - add contribution of low, middle and high clouds to the
!       clear sky portion.
!
!     process narrative:  flux3 - located in the flux3 sdf in dnxm
!     =================
!      
!     references: 
!     ==========
!     dr idso's paper in the j. of geophys. research,   
!     no 74, pp 5397-5403.  
!  
!     dr r.f.wachtmann's paper in the digest of preprints,  
!     topical meeting on remote sensing of the atmosphere,  
!     anaheim,ca, optical society of america, entitled, 
!     "expansion of atmospheric temperature-moisture
!     profiles in empirical orthogonal functions for remote 
!     sensing applications", 1975   
!  
!     called from:   
!     ===========
!     getagrlw   
!
!     interface:
!     =========  
!
!       input variables  
!       ---------------  
!       sfctmp, e, iclamt
!    
!       output variables 
!       ----------------
!       rldown
!
!     files accessed:
!     ============== 
!     filename/unit#           r/w               description 
!     ------------------------ --- -------------------------------------
!     none
!
!     remarks: none
!     =======
!
!     variables: 
!     =========  
!     label     .......................description......................
!  
!     cldfrt    fraction of low, mid, and high cloud amounts
!     clrsky    irradiance contribution from the clear sky  
!     e         sfc vapor pressure (pascals) 
!     emb       sfc vapor pressure (millibars)   
!     emissa    idso clr sky emissivity (all wavelengths)
!     emissb    sasc clr sky emissivity (adjusted emissa value)  
!     hcterm    high cloud emissivity term   
!     iclamt    low, mid, and high cloud amounts in percent
!     lcterm    low cloud emissivity term   
!     mcterm    mid cloud emissivity term  
!     rldown    net downward longwave irradiance.
!     sfctmp    sfc temperature ( deg k )
!     sigma     stefan boltzman constant
!     zh        high cloud height coefficient   
!     zl        low cloud height coefficient
!     zm        mid cloud height coefficient
!
!     updates
!     =======
!     15 may 1988  initial version........................capt rice/sddc  
!     07 sep 1999  ported to ibm sp-2.  added intent attributes to
!                  arguments.  updated prolog..............mr gayno/dnxm
!     25 oct 2001  implement in LDAS.....................jesse meng/ncep
!-----------------------------------------------------------------------
!----------------------------------------------------------------------- 

      implicit none

      real,     intent(in)         :: iclamt  ( 3 ) 

      real                         :: cldfrt  ( 3 )
      real                         :: clrsky   
      real,     intent(in)         :: e
      real                         :: emb  
      real                         :: emissa   
      real                         :: emissb   
      real                         :: hcterm   
      real                         :: lcterm   
      real                         :: mcterm   
      real,     intent(out)        :: rldown   
      real,     intent(in)         :: sfctmp   
      real,     parameter          :: sigma = 5.67e-08
      real,     parameter          :: zh    = 8.0
      real,     parameter          :: zl    = 1.3
      real,     parameter          :: zm    = 3.1  

!     ------------------------------------------------------------------
!     executable code starts here...compute the cloud amount 
!     in fraction of overcast (.1 to 1.0). 
!     ------------------------------------------------------------------

      cldfrt(1) =  iclamt(1)  / 100.0
      cldfrt(2) =  iclamt(2)  / 100.0
      cldfrt(3) =  iclamt(3)  / 100.0

!     ------------------------------------------------------------------
!     convert vapor pressure units from pascals to millibars for use
!     in determining emissivities.  
!     ------------------------------------------------------------------

      emb = e * 0.01

!     ------------------------------------------------------------------
!     compute the effective clr sky emissivity for all wavelengths  
!     (emissa) using idso's equation.   
!     ------------------------------------------------------------------

      emissa = 0.700 + (5.95e-5 * emb * exp(1500 / sfctmp))

!     ------------------------------------------------------------------
!     use emissa in wachtmann's model for sky irradiance to calc a  
!     resultant longwave downward radiation value.  first calc a sasc   
!     emmisivity (emissb), which is an adjusted idso emmisivity.
!     then use emissb to calculate the blackbody irradiance of the  
!     clear sky (the 1st term of wachtmann's equation).
!     ------------------------------------------------------------------

      emissb = -0.792 + (3.161 * emissa) - (1.573 * emissa * emissa) 
      clrsky =  emissb * sigma * ( sfctmp * sfctmp * sfctmp * sfctmp )

!     ------------------------------------------------------------------
!     now compute the irradiance contribution from the low, middle, 
!     and hi cloud layers (the 2nd thru 4th terms in wachtmann' eqn).
!     ------------------------------------------------------------------

      lcterm = (80.0 - (5.0 * zl)) * cldfrt(1) 
      mcterm = (80.0 - (5.0 * zm)) * (1.0 - cldfrt(1)) * cldfrt(2)  
      hcterm = (80.0 - (5.0 * zh)) * (1.0 - cldfrt(1)) * &
               (1.0 - cldfrt(2)) * cldfrt(3)   

!     ------------------------------------------------------------------
!     put it all together to get a resultant downwrd longwave irrad.
!     ------------------------------------------------------------------

      rldown = clrsky + hcterm + mcterm + lcterm

      return

      end   
