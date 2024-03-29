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
module LDT_spatialDownscalingMod
!BOP
!
! !MODULE: LDT_spatialDownscalingMod
! 
! !DESCRIPTION:
! 
! !REVISION HISTORY: 
!  6 Nov 2012: Sujay Kumar; initial specification

  use ESMF
  use LDT_logMod
  use LDT_constantsMod
  use LDT_coreMod
  use LDT_FORC_AttributesMod
  use LDT_fileIOMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LDT_lapseRateCorrection
  PUBLIC :: LDT_slopeAspectCorrection
!  PUBLIC :: LDT_init_pcpclimo
!  PUBLIC :: LDT_init_pcpclimo_native
!  PUBLIC :: LDT_generatePcpClimoRatioField
!  PUBLIC :: LDT_pcpClimoDownscaling
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!----------------------------------------------------------------------------- 
#if 0
  PUBLIC :: LDT_pcpclimo
 
  type, public :: pcpclimo_dec_type
     real, allocatable :: pcpdata1_native(:)
     real, allocatable :: pcpdata2_native(:)
     real, allocatable :: pcpdata1(:)
     real, allocatable :: pcpdata2(:)
     integer       :: month1
     integer       :: month2
  end type pcpclimo_dec_type

  type(pcpclimo_dec_type), allocatable :: LDT_pcpclimo(:,:)
#endif

!EOP
contains


!BOP
! !ROUTINE: LDT_lapseRateCorrection
! \label{LDT_lapseRateCorrection}
!
! !REVISION HISTORY:
!  11 Apr 2000: Brian Cosgrove; Initial Code
!  12 May 2000: Brian Cosgrove; Corrected for zero humidities
!  25 Jan 2001: Matt Rodell; Compute number of input and output
!		grid points, use to allocate local arrays
!  15 Mar 2001: Jon Gottschalck; if-then to handle negative vapor
!		pressures in long wave correction
!  14 Nov 2003: Sujay Kumar; Adopted in LIS
!  14 Oct 2014: KR Arsenault; Adopted in LDT
!
! !INTERFACE:
  subroutine LDT_lapseRateCorrection(nest, modelelev, LDT_FORC_Base_State)
! !USES:
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: nest
    real                :: modelelev(LDT_rc%ngrid(nest))
    type(ESMF_State)    :: LDT_FORC_Base_State

! !DESCRIPTION:
!  Corrects the lowest model level Temperature, Pressure, 
!  Humidity and Longwave Radiation forcing values for 
!  differences in elevation between the LIS running 
!  grid and the native forcing grid. The corrections are based on the 
!  lapse-rate and hypsometric adjustments to these variables described
!  in Cosgrove et. al (2003). 
!  
!  Cosgrove, B.A. et.al, Real-time and retrospective forcing in the 
!  North American Land Data Assimilation (NLDAS) project, Journal of 
!  Geophysical Research, 108(D22), 8842, DOI: 10.1029/2002JD003118, 2003.  
!
!  The arguments are: 
!  \begin{description}
!   \item [nest]
!     index of the domain or nest.
!   \item [findex]
!     index of the forcing dataset
!  \end{description}
!
!EOP
    integer:: t
    real :: force_tmp,force_hum,&
         force_lwd,force_prs
    real :: elevdiff
    integer, parameter :: bb=2016
    integer :: index
    real :: mee,mfe,ee,fe
    real :: lapse, rdry, ratio
    real :: esat,qsat,rh,fesat,fqsat,femiss,emiss
    real :: tcforce,pcforce,hcforce,lcforce,tbar
    integer            :: status
    type(ESMF_Field)   :: tmpField,q2Field,lwdField,psurfField
    real, pointer      :: tmp(:),q2(:),lwd(:),psurf(:)
    
    rdry = 287.
    lapse = -0.0065
    
    call ESMF_StateGet(LDT_FORC_Base_State,&
         trim(LDT_FORC_Tair%varname(1)),tmpField,&
         rc=status)
    call LDT_verify(status)
    
    call ESMF_FieldGet(tmpField,localDE=0, farrayPtr=tmp,rc=status)
    call LDT_verify(status)
    
    
    call ESMF_StateGet(LDT_FORC_Base_State,&
         trim(LDT_FORC_Qair%varname(1)),q2Field,&
         rc=status)
    call LDT_verify(status)
    
    call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
    call LDT_verify(status)

    
    call ESMF_StateGet(LDT_FORC_Base_State,&
         trim(LDT_FORC_LWdown%varname(1)),lwdField,&
         rc=status)
    call LDT_verify(status)
    
    call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
    call LDT_verify(status)
    
    
    call ESMF_StateGet(LDT_FORC_Base_State,&
         trim(LDT_FORC_Psurf%varname(1)),psurfField,&
         rc=status)
    call LDT_verify(status)
    
    call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
    call LDT_verify(status)
    

    do t=1,LDT_rc%ntiles(nest)
       if(tmp(t).gt.0) then 
          force_tmp = tmp(t)
          force_hum = q2(t)
          force_lwd = lwd(t)
          force_prs = psurf(t)
          index = LDT_domain(nest)%tile(t)%index
          elevdiff = LDT_domain(nest)%tile(t)%elev-&
               modelelev(index)
          tcforce=force_tmp+(lapse*elevdiff)
          tbar=(force_tmp+tcforce)/2.
          pcforce=force_prs/(exp((LDT_CONST_G*elevdiff)/(rdry*tbar)))
          if (force_hum .eq. 0) force_hum=1e-08
          ee=(force_hum*force_prs)/0.622               
          esat=611.2*(exp((17.67*(force_tmp-LDT_CONST_TKFRZ))/&
               ((force_tmp-LDT_CONST_TKFRZ)+243.5)))
          qsat=(0.622*esat)/(force_prs-(0.378*esat))
          rh=(force_hum/qsat)*100.
          fesat=611.2*(exp((17.67*(tcforce-LDT_CONST_TKFRZ))/ &
               ((tcforce-LDT_CONST_TKFRZ)+243.5)))
          fqsat=(0.622*fesat)/(pcforce-(0.378*fesat))
          hcforce=(rh*fqsat)/100.
          fe=(hcforce*pcforce)/0.622
          mee=ee/100.
          mfe=fe/100.
       !----------------------------------------------------------------------
       ! correct for negative vapor pressure at very low temperatures at
       ! high latitudes
       !----------------------------------------------------------------------
          if (mee .le. 0) mee = 1e-08
          if (mfe .le. 0) mfe = 1e-08
          emiss  =1.08*(1-exp(-mee**(force_tmp/bb)))
          femiss =1.08*(1-exp(-mfe**(tcforce/bb)))
          ratio=(femiss*(tcforce**4))/(emiss*(force_tmp**4))
          lcforce=force_lwd*ratio
          
          tmp(t)   = tcforce
          q2(t)    = hcforce
          lwd(t)   = lcforce
          psurf(t) = pcforce
       endif
    end do
  end subroutine LDT_lapseRateCorrection


!BOP
! !ROUTINE: LDT_slopeAspectCorrection
! \label{LDT_slopeAspectCorrection}
!
! !REVISION HISTORY:
!  1 Feb 2011: Sujay Kumar; Initial implementation 
!
! !INTERFACE:
  subroutine LDT_slopeAspectCorrection(nest,LDT_FORC_Base_State)
! !USES:
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: nest
    type(ESMF_State)    :: LDT_FORC_Base_State

! !DESCRIPTION:
!
!  This subroutine adjusts the downward shortwave radiation to correct 
!  for the influence of topographic slope and aspect. The method separates
!  the downward shortwave values into direct and diffuse components first. 
!  The topographic corrections are then applied to the direct component. 
!  
!  Adapted from Dingman, S.L., Physical Hydrology - 2nd ed.- Appendix E, 
!  Prentice Hall, 2002. The separation of the downward shortwave forcing
!  into direct and diffuse components is adopted from HySSIB LSM. The
!  slope and aspect based correction of the direct shortwave term is 
!  adopted from FASST LSM. 
!  
!  The arguments are: 
!  \begin{description}
!   \item [nest]
!     index of the domain or nest.
!   \item [findex]
!     index of the forcing dataset
!  \end{description}
!
!EOP
    integer:: t
    integer            :: status
    type(ESMF_Field)   :: swdField
    real, pointer      :: swd(:)
    real               :: sunang
    real               :: cloud
    real               :: difrat, vnrat
    real               :: swddirect, swddiffuse
    real               :: lhour, cosz
    integer            :: zone
    real               :: decl, omega
    real               :: lat
    integer            :: index
    real               :: solzen, delaz, sdircorr
    real               :: saz, costheta, szen, ha, phi,cosphi
    real               :: tst, time_offset, mody, eqtime, thour
    real*8             :: fyear
    real               :: hangle, aslope, deg2rad
    real               :: slope, aspect

    call ESMF_StateGet(LDT_FORC_Base_State,&
         trim(LDT_FORC_SWdown%varname(1)),swdField,&
         rc=status)
    call LDT_verify(status)

    call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
    call LDT_verify(status)

    ! separation into direct and diffuse components adapted from HySSIB
    ! (Xue et al. Journal of Climate 1991, Sellers et al. JAS, 1986)
    !
    do t=1,LDT_rc%ntiles(nest)
       if(swd(t).gt.0) then 
          index = LDT_domain(nest)%tile(t)%index

          call LDT_localtime(LDT_rc%gmt, LDT_domain(nest)%grid(index)%lon, &
               lhour, zone)
          call coszenith( LDT_domain(nest)%grid(index)%lon, &
               LDT_domain(nest)%grid(index)%lat, &
               lhour, zone, LDT_rc%doy, cosz, decl, omega)

          sunang = max(cosz,0.01764)
          cloud  = (1160.0*sunang - swd(t))/(963.0*sunang)
          cloud  = max(cloud, 0.0)
          cloud  = min(cloud, 1.0)

          if(abs(sunang-0.0223).gt.0.0001) then 
             difrat = 0.0604 / (sunang - 0.0223) + 0.0683       
             difrat = max(difrat, 0.0)
             difrat = min(difrat, 1.0)
          else
             difrat = 1.0
          endif

          difrat = difrat + (1.0 - difrat)*cloud

          vnrat = (580.0 - cloud*464.0)/((580.0 - cloud*499.0) + & 
               (580.0 - cloud*464.0))
          swddirect  = swd(t) * &
               ((1.0 - difrat) *vnrat + (1.0 - difrat)*(1.0 - vnrat))
          swddiffuse = swd(t) * (difrat*vnrat + difrat*(1.0 - vnrat))


          deg2rad = LDT_CONST_PI/180.0
          solzen = acos(cosz)
          slope  = LDT_domain(nest)%tile(t)%slope/deg2rad
          aspect = LDT_domain(nest)%tile(t)%aspect/deg2rad
          lat = LDT_domain(nest)%grid(index)%lat*LDT_CONST_PI/180.0

          !        hangle = 15*(12-lhour)*deg2rad
          !        azm_x = sin(hangle)*cos(decl)
          !        azm_y = (-(cos(hangle))*cos(decl)*sin(lat))+&
          !             (cos(lat)*sin(decl))
          !        saz = atan(azm_x/azm_y)
          ! Computation of saz
          thour = LDT_rc%hr
          if((thour-24) <=0) thour = 0
          mody = LDT_rc%yr - dint(LDT_rc%yr*2.5d-1)*4d0
          if(abs(mody).gt.0) then 
             fyear = (2.0*LDT_CONST_PI/365.0)*&
                  (LDT_rc%doy-1.0+(thour-12.0)/24.0)
          else
             fyear = (2.0*LDT_CONST_PI/366.0)*&
                  (LDT_rc%doy-1.0+(thour-12.0)/24.0)
          endif
          !calculate the equation of time in minutes
          eqtime = 229.18d0*(7.5d-5 + 1.868d-3*dcos(fyear)            & 
               - 3.2077d-2*dsin(fyear) - 1.4615d-2*dcos(2d0*fyear) &
               - 4.0849d-2*dsin(2d0*fyear))

          !calculate the true solar time 
          !CHECK if this the timeoffset = lhour - LDT_rc%hr is correct
          time_offset = eqtime + 4.0*LDT_domain(nest)%grid(index)%lon + &
               60.0*abs(lhour-LDT_rc%hr)
          tst = thour*60.0 + LDT_rc%mn + time_offset

          !solar hour angle
          ha = (tst*0.25 - 180.0)*deg2rad

          cosphi = sin(lat)*sin(decl)+cos(lat)*cos(decl)*cos(ha)
          phi = acos(cosphi)
          szen = phi/deg2rad

          costheta = (sin(lat)*cosphi-sin(decl))/(cos(lat)*sin(phi))

          !avoid floating point errors
          if(abs(costheta).gt.1) then 
             if(costheta > 0 ) then 
                costheta = 1
             else
                costheta = -1
             endif
          endif

          if(lat.ge.0.0) then 
             if(ha.lt.0.0) then  
                saz = 180.0-acos(costheta)/deg2rad
             else
                saz = 180.0+acos(costheta)/deg2rad
             endif
          else
             if(ha.lt.0.0) then 
                saz = acos(costheta)/deg2rad
             else
                saz = 360.0-acos(costheta)/deg2rad
             endif
          endif

          if(solzen .ge. 90*deg2rad) then 
             swd(t) = 0.0
          elseif(solzen.lt.90*deg2rad .and. swd(t).gt.0.0) then 
             !CHECK: slope here is in deg? -- how was it computed? 
             aslope = max(0.0, min(1.57,slope*deg2rad))

             if(aslope > 0.0 .and. aslope.lt. 90*deg2rad) then
                !aspect should be in degrees... 
                delaz = abs(aspect - saz) * deg2rad 
                if(solzen <= 85*deg2rad) then 
                   sdircorr = (sin(solzen)*sin(aslope)*cos(delaz))/&
                        cos(solzen)
                else
                   sdircorr = (sin(solzen)*sin(aslope)*cos(delaz))/&
                        cos(85.0)
                endif
                swddirect = max(0.0, swddirect*(cos(aslope) + sdircorr))
             endif
          endif
          !CHECK: contribution of the albedo term? 
          swd(t) = swddirect + swddiffuse
       endif
    enddo

  end subroutine LDT_slopeAspectCorrection

#if 0

  subroutine LDT_init_pcpclimo()

    integer :: n,m

    allocate(LDT_pcpclimo(LDT_rc%nnest,LDT_rc%nmetforc))

    do n=1,LDT_rc%nnest
       do m=1,LDT_rc%nmetforc
          if(LDT_rc%pcp_downscale(m).ne.0) then
             
             allocate(LDT_pcpclimo(n,m)%pcpdata1(&
                  LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(LDT_pcpclimo(n,m)%pcpdata2(&
                  LDT_rc%lnc(n)*LDT_rc%lnr(n)))

             LDT_pcpclimo(n,m)%pcpdata1 = LDT_rc%udef
             LDT_pcpclimo(n,m)%pcpdata2 = LDT_rc%udef

             LDT_pcpclimo(n,m)%month1 = -1
             LDT_pcpclimo(n,m)%month2 = -1
          endif
       enddo
    enddo

  end subroutine LDT_init_pcpclimo

  subroutine LDT_init_pcpclimo_native(n,m,nc,nr)
    
    integer      :: n 
    integer      :: m
    integer      :: nc
    integer      :: nr
    
    allocate(LDT_pcpclimo(n,m)%pcpdata1_native(nc*nr))         
    allocate(LDT_pcpclimo(n,m)%pcpdata2_native(nc*nr))
    
    LDT_pcpclimo(n,m)%pcpdata1_native = LDT_rc%udef
    LDT_pcpclimo(n,m)%pcpdata2_native = LDT_rc%udef
    
  end subroutine LDT_init_pcpclimo_native
! 
!BOP
! !ROUTINE: LDT_generatePcpClimoRatioField
! \label{LDT_generatePcpClimoRatioField}
!
! !REVISION HISTORY:

!
! !INTERFACE:
  subroutine LDT_generatePcpClimoRatioField(n, m, forcing_name, month, &
       input_size, input_data, input_bitmap)
! !USES:
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

    implicit none
! !ARGUMENTS: 
    integer                 :: n 
    integer                 :: m
    character(len=*)        :: forcing_name
    integer                 :: month
    integer                 :: input_size
    real                    :: input_data(input_size)
    logical*1               :: input_bitmap(Input_size)

!
! !DESCRIPTION:
!EOP
    logical                :: read_flag
    logical                :: file_exists
    integer                :: nid,ncId,nrId,pcpId,pcphiresId
    integer                :: nc,nr,c,r,t
    real,     allocatable      :: pcpclimo_native(:,:)
    real,     allocatable      :: pcpclimo_hires(:,:)
    real,     allocatable      :: pcpclimo_hires1(:,:)
    real,     allocatable      :: pcpclimo(:,:)


!Step 1: read climo data in the forcing grid. 
    
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )

    read_flag = .false. 
    
    if(LDT_pcpclimo(n,m)%month1.ne.month) then !time to read both bookends
       LDT_pcpclimo(n,m)%month1 = month
       if(month+1.gt.12) then           
          LDT_pcpclimo(n,m)%month2 = 1
       else
          LDT_pcpclimo(n,m)%month2 = month+1
       endif
       read_flag = .true. 
    endif

    if(read_flag) then 
       inquire(file=LDT_rc%paramfile(n), exist=file_exists)
       if(file_exists) then 
          write(LDT_logunit,*) 'Reading precip climatology data for months ',&
               LDT_pcpclimo(n,m)%month1,' and ',LDT_pcpclimo(n,m)%month2
          
          call LDT_verify(nf90_open(path=LDT_rc%paramfile(n),&
               mode=NF90_NOWRITE,ncid=nid),&
               'Error in nf90_open in LDT_generatePcpClimoRatioField')
          
          call LDT_verify(nf90_inq_dimid(nid,&
               "east_west_"//trim(forcing_name),&
               ncId),&
               'Error in nf90_inq_dimid in LDT_generatePcpClimoRatioField')
          
          call LDT_verify(nf90_inq_dimid(nid,&
               "north_south_"//trim(forcing_name),&
               nrId), 'Error in nf90_inq_dimid in LDT_generatePcpClimoRatioField')
          
          call LDT_verify(nf90_inquire_dimension(nid,ncId, len=nc),&
               'Error in nf90_inquire_dimension in LDT_generatePcpClimoRatioField')
          
          call LDT_verify(nf90_inquire_dimension(nid,nrId, len=nr),&
               'Error in nf90_inquire_dimension in LDT_generatePcpClimoRatioFieldy')
          
          if(nc*nr.ne.input_size) then 
             write(LDT_logunit,*) 'The input dimensions of the '//trim(forcing_name)
             write(LDT_logunit,*) '(',input_size,')'
             write(LDT_logunit,*) 'does not match the dimensions in the LIS parameter file'
             write(LDT_logunit,*) '(',nc*nr,')'
             call LDT_endrun()
          endif
          
          allocate(pcpclimo_native(nc,nr))
          
          call LDT_verify(nf90_inq_varid(nid,"PPTCLIM_"//trim(forcing_name),pcpId),&
               'PPTCLIM_'//trim(forcing_name)//' field not found in LIS param file')
          
          call LDT_verify(nf90_get_var(nid,pcpId,pcpclimo_native, &
               start=(/1,1,LDT_pcpclimo(n,m)%month1/), count=(/nc,nr,1/)),&
               'nf90_get_var failed for PPTCLIM_'//trim(forcing_name))
          
          !convert to 1d field
          do r=1,nr
             do c=1,nc
                LDT_pcpclimo(n,m)%pcpdata1_native(c+(r-1)*nc) = pcpclimo_native(c,r)
             enddo
          enddo

          call LDT_verify(nf90_get_var(nid,pcpId,pcpclimo_native, &
               start=(/1,1,LDT_pcpclimo(n,m)%month2/), count=(/nc,nr,1/)),&
               'nf90_get_var failed for PPTCLIM_'//trim(forcing_name))
          
          !convert to 1d field
          do r=1,nr
             do c=1,nc
                LDT_pcpclimo(n,m)%pcpdata2_native(c+(r-1)*nc) = pcpclimo_native(c,r)
             enddo
          enddo
          deallocate(pcpclimo_native)

          allocate(pcpclimo_hires(LDT_rc%pnc(n),LDT_rc%pnr(n)))
          allocate(pcpclimo_hires1(LDT_rc%gnc(n),LDT_rc%gnr(n)))
          allocate(pcpclimo(LDT_rc%lnc(n),LDT_rc%lnr(n)))

          call LDT_verify(nf90_inq_varid(nid, "PPTCLIM",pcpHiresId),&
               'PPTCLIM field not found in LIS param file')
          call LDT_verify(nf90_get_var(nid,pcphiresId,pcpclimo_hires,&
               start=(/1,1,LDT_pcpclimo(n,m)%month1/), &
               count=(/LDT_rc%gnc(n),LDT_rc%gnr(n),1/)),&
               'n90_get_var failed for PPTCLIM in LIS param file')
          
          call LDT_convertParamDataToLocalDomain(n,&
               pcpclimo_hires, pcpclimo_hires1)

          pcpclimo(:,:) = pcpclimo_hires1(&
               LDT_ews_halo_ind(n,LDT_localPet+1):&         
               LDT_ewe_halo_ind(n,LDT_localPet+1), &
               LDT_nss_halo_ind(n,LDT_localPet+1): &
               LDT_nse_halo_ind(n,LDT_localPet+1))
          
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                LDT_pcpclimo(n,m)%pcpdata1(c+(r-1)*LDT_rc%lnc(n)) = &
                     pcpclimo(c,r)
             enddo
          enddo

          call LDT_verify(nf90_get_var(nid,pcphiresId,pcpclimo_hires,&
               start=(/1,1,LDT_pcpclimo(n,m)%month2/), &
               count=(/LDT_rc%gnc(n),LDT_rc%gnr(n),1/)),&
               'n90_get_var failed for PPTCLIM in LIS param file')
          
          call LDT_convertParamDataToLocalDomain(n,&
               pcpclimo_hires, pcpclimo_hires1)

          pcpclimo(:,:) = pcpclimo_hires1(&
               LDT_ews_halo_ind(n,LDT_localPet+1):&         
               LDT_ewe_halo_ind(n,LDT_localPet+1), &
               LDT_nss_halo_ind(n,LDT_localPet+1): &
               LDT_nse_halo_ind(n,LDT_localPet+1))
          
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                LDT_pcpclimo(n,m)%pcpdata2(c+(r-1)*LDT_rc%lnc(n)) = &
                     pcpclimo(c,r)
             enddo
          enddo

          deallocate(pcpclimo_hires)
          deallocate(pcpclimo_hires1)
          deallocate(pcpclimo)

          call LDT_verify(nf90_close(nid),&
               'nf90_close failed in LDT_generatePcpClimoRatioField')
       
          write(LDT_logunit,*) 'Done reading precip climatology data '
       endif
    endif

    do t=1,input_size
       if(.not.input_bitmap(t)) then         
          input_data(t) = -9999.0
       endif
    enddo


    !step 2: compute ratio field
    do t=1,input_size
       if(input_bitmap(t)) then 
          if(month.eq.LDT_pcpclimo(n,m)%month1) then 
             if(LDT_pcpclimo(n,m)%pcpdata1_native(t).ne.0) then 
                input_data(t) = input_data(t)/&
                     LDT_pcpclimo(n,m)%pcpdata1_native(t)
             else
                input_data(t) = 0.0
             endif
          elseif(month.eq.LDT_pcpclimo(n,m)%month2) then 
             if(LDT_pcpclimo(n,m)%pcpdata2_native(t).ne.0) then 
                input_data(t) = input_data(t)/&
                     LDT_pcpclimo(n,m)%pcpdata2_native(t)
             else
                input_data(t) = 0.0
             endif
          endif
       else
          input_data(t) = LDT_rc%udef
       endif
    enddo
    
#endif
    
  end subroutine LDT_generatePcpClimoRatioField
  
  subroutine LDT_pcpClimoDownscaling(n,m,month, &
       output_size, output_data, output_bitmap)
    
    integer                   :: n 
    integer                   :: m
    integer                   :: month
    integer                   :: output_size
    real                      :: output_data(output_size)
    logical*1                 :: output_bitmap(output_size)

    integer               :: t
!multiply by the high-res climo field. 

  do t=1,output_size
     if(output_bitmap(t)) then 
        if(month.eq.LDT_pcpclimo(n,m)%month1) then 
           if(LDT_pcpclimo(n,m)%pcpdata1(t).ne.-9999.0) then 
              output_data(t) = output_data(t)*&
                   LDT_pcpclimo(n,m)%pcpdata1(t)
           else
              output_data(t) = LDT_rc%udef
           endif

        elseif(month.eq.LDT_pcpclimo(n,m)%month2) then 
           output_data(t) = output_data(t)*&
                LDT_pcpclimo(n,m)%pcpdata2(t)
        endif
     endif
  enddo
 
 end subroutine LDT_pcpClimoDownscaling

#endif

end module LDT_spatialDownscalingMod

