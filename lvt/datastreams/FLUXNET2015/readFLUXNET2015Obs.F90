!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readtemplateObs
! \label{readtemplateObs}
!
! !INTERFACE: 
  subroutine readFLUXNET2015Obs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use LVT_logMod
  use FLUXNET2015_obsMod
  use map_utils
!
! !DESCRIPTION: 
! 
!   This routine provides the data processing code for the FLUXNET2015
!   station dataset. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Feb 2017: Sujay Kumar, Initial specification
! 
!EOP
!--------------------------------------------------------------------------
  implicit none

  integer         :: source 

  integer :: i 
  type(ESMF_Time) :: fluxnettime1, fluxnettime2
  integer         :: t, st, et
  integer         :: yr, mo, da, hr, mn, ss
  real*8          :: lis_prevtime
  integer         :: c,r,stn_col, stn_row
  real            :: col,row
  real            :: gmt
  integer         :: doy
  real            :: qle(LVT_rc%lnc, LVT_rc%lnr)
  real            :: qh(LVT_rc%lnc, LVT_rc%lnr)
  integer         :: status
  real            :: timenow

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)

  qle = LVT_rc%udef
  qh = LVT_rc%udef

  ! Start of new data
  if(FLUXNET2015obs(source)%startflag.or.&
       (FLUXNET2015obs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(FLUXNET2015obs(source)%starttime, &
          yy=LVT_rc%dyr(source), &
          mm=1, dd=1, h=0, m=0,s = 0, calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'fluxnet starttime set failed')

     FLUXNET2015obs(source)%startflag = .false. 
     FLUXNET2015obs(source)%yr = LVT_rc%dyr(source)

     do i = 1,FLUXNET2015obs(source)%n_stns
        call read_fluxnet_station(source,i)
     end do
  endif

  if(mod(timenow, 1800.0).eq.0.0) then 
     call ESMF_TimeSet(fluxnettime1, yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), &
          m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'fluxnettime1 set failed')
     
     t = nint((fluxnettime1 - FLUXNET2015obs(source)%starttime)/&
          FLUXNET2015obs(source)%timestep)+1

     do i=1,FLUXNET2015obs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, FLUXNET2015obs(source)%stnlat(i), &
             FLUXNET2015obs(source)%stnlon(i), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if((stn_col.ge.1.and.stn_col.le.LVT_rc%lnc).and.&
             (stn_row.ge.1.and.stn_row.le.LVT_rc%lnr)) then 
           if(FLUXNET2015obs(source)%qle(i,t).ne.LVT_rc%udef) then 
              qle(stn_col, stn_row) =  FLUXNET2015obs(source)%qle(i,t)
           end if
           
           if(FLUXNET2015obs(source)%qh(i,t).ne.LVT_rc%udef) then 
              qh(stn_col, stn_row) = FLUXNET2015obs(source)%qh(i,t)
           end if
        
        endif
     end do
  else
     qle = LVT_rc%udef
     qh = LVT_rc%udef

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_qle, source, qle,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qh, source,qh,vlevel=1,units="W/m2")


end subroutine readFLUXNET2015Obs

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
subroutine read_fluxnet_station(source, stn_index)
! 
! !USES:
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use FLUXNET2015_obsMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:

! Opens each station file and reads line by line and extracts the 
! relevant variables into the datastructures for later use. 
! 
!  Note that the code only processes the water and energy measurements currently
!
!  TA_F_MDS,,"Air temperature, gapfilled using MDS method"
!  TA_F_MDS_QC,,Quality flag for TA_F_MDS
!  TA_F_MDS_NIGHT,,Average nighttime TA_F_MDS
!  TA_ERA,,"Air temperature, downscaled from ERA, linearly regressed using measured only site data"
!  TA_ERA_NIGHT,,Average nighttime TA_ERA
!  TA_ERA_NIGHT_SD,,Standard deviation for TA_ERA_NIGHT
!  TA_F,,"Air temperature, consolidated from TA_F_MDS and TA_ERA"
!  TA_F_QC,,Quality flag for TA_F
!  TA_F_NIGHT_QC,,Quality flag for TA_F_NIGHT
!  SW_IN_POT,,"Shortwave radiation, incoming, potential (top of atmosphere)"
!  SW_IN_F_MDS,,"Shortwave radiation, incoming, gapfilled using MDS (negative values set to zero, e.g., negative values from instrumentation noise)"
!  SW_IN_F_MDS_QC,,Quality flag for SW_IN_F_MDS
!  SW_IN_ERA,,"Shortwave radiation, incoming, downscaled from ERA, linearly regressed using measured only site data (negative values set to zero)"
!  SW_IN_F,,"Shortwave radiation, incoming consolidated from SW_IN_F_MDS and SW_IN_ERA (negative values set to zero)"
!  SW_IN_F_QC,,Quality flag for SW_IN_F
!  LW_IN_F_MDS,,"Longwave radiation, incoming, gapfilled using MDS"
!  LW_IN_F_MDS_QC,,Quality flag for LW_IN_F_MDS
!  LW_IN_ERA,,"Longwave radiation, incoming, downscaled from ERA, linearly regressed using measured only site data"
!  LW_IN_F,,"Longwave radiation, incoming, consolidated from LW_IN_F_MDS and LW_IN_ERA"
!  LW_IN_F_QC,,Quality flag for LW_IN_F
!  LW_IN_JSB,,"Longwave radiation, incoming, calculated from TA_F_MDS, SW_IN_F_MDS, VPD_F_MDS and SW_IN_POT using the JSBACH algorithm (Sonke Zaehle)"
!  LW_IN_JSB_QC,,Quality flag for LW_IN_JSB
!  LW_IN_JSB_ERA,,"Longwave radiation, incoming, downscaled from ERA, linearly regressed using site level LW_IN_JSB calculated from measured only drivers"
!  LW_IN_JSB_F,,"Longwave radiation, incoming, consolidated from LW_IN_JSB and LW_IN_JSB_ERA"
!  LW_IN_JSB_F_QC,,Quality flag for LW_IN_JSB_F
!  VPD_F_MDS,,"Vapor Pressure Deficit, gapfilled using MDS"
!  VPD_F_MDS_QC,,Quality flag for VPD_F_MDS
!  VPD_ERA,,"Vapor Pressure Deficit, downscaled from ERA, linearly regressed using measured only site data"
!  VPD_F,,Vapor Pressure Deficit consolidated from VPD_F_MDS and VPD_ERA
!  VPD_F_QC,,Quality flag for VPD_F
!  PA,,Atmospheric pressure
!  PA_ERA,,"Atmospheric pressure, downscaled from ERA, linearly regressed using measured only site data"
!  PA_F,,Atmospheric pressure consolidated from PA and PA_ERA
!  PA_F_QC,,Quality flag for PA_F
!  P,,Precipitation
!  P_ERA,,"Precipitation, downscaled from ERA, linearly regressed using measured only site data"
!  P_F,,Precipitation consolidated from P and P_ERA
!  P_F_QC,,Quality flag for P_F
!  WS_ERA,,"Wind speed, downscaled from ERA, linearly regressed using measured only site data"
!  WS_F,,"Wind speed, consolidated from WS and WS_ERA"
!  WS_F_QC,,Quality flag of WS_F
!  WD,,Wind direction
!  WD,,Relative humidity, range 0-100
!  USTAR,,Friction velocity
!  USTAR_QC,,Quality flag of USTAR
!  NETRAD,,Net radiation
!  NETRAD_QC,,Quality flag of NETRAD
!  PPFD_IN,,"Photosynthetic photon flux density, incoming"
!  PPFD_IN_QC,,Quality flag of PPFD_IN
!  PPFD_DIF,,"Photosynthetic photon flux density, diffuse incoming"
!  PPFD_DIF_QC,,Quality flag of PPFD_DIF
!  PPFD_OUT,,"Photosynthetic photon flux density, outgoing"
!  PPFD_OUT_QC,,Quality flag of PPFD_OUT
!  SW_DIF,,"Shortwave radiation, diffuse incoming"
!  SW_DIF_QC,,Quality flag of SW_DIF
!  SW_OUT,,"Shortwave radiation, outgoing"
!  SW_OUT_QC,,Quality flag of SW_OUT
!  LW_OUT,,"Longwave radiation, outgoing"
!  LW_OUT_QC,,Quality flag of LW_OUT
!  CO2_F_MDS,,"CO2 mole fraction, gapfilled with MDS"
!  CO2_F_MDS_QC,,Quality flag for CO2_F_MDS
!  TS_F_MDS_#,,"Soil temperature, gapfilled with MDS (numeric index ""#"" increases with the depth, 1 is shallowest)"
!  TS_F_MDS_#_QC,,Quality flag for TS_F_MDS_#
!  SWC_F_MDS_#,,"Soil water content, gapfilled with MDS (numeric index ""#"" increases with the depth, 1 is shallowest)"
!  SWC_F_MDS_#_QC,,Quality flag for SWC_F_MDS_#
!  G_F_MDS,,Soil heat flux
!  G_F_MDS_QC,,Quality flag of G_F_MDS
!  LE_F_MDS,,"Latent heat flux, gapfilled using MDS method"
!  LE_F_MDS_QC,,"Quality flag for LE_F_MDS, LE_CORR, LE_CORR25, and LE_CORR75"
!  LE_CORR,,"Latent heat flux, corrected LE_F_MDS by energy balance closure correction factor"
!  LE_CORR_25,,"Latent heat flux, corrected LE_F_MDS by energy balance closure correction factor, 25th percentile"
!  LE_CORR_75,,"Latent heat flux, corrected LE_F_MDS by energy balance closure correction factor, 75th percentile"
!  LE_RANDUNC,,"Random uncertainty of LE, from measured only data"
!  LE_RANDUNC_METHOD,,Method used to estimate the random uncertainty of LE
!  LE_RANDUNC_N,,Number of half-hour data points used to estimate the random uncertainty of LE
!  LE_CORR_JOINTUNC,,Joint uncertainty estimation for LE
!  H_F_MDS,,"Sensible heat flux, gapfilled using MDS method"
!  H_F_MDS_QC,,"Quality flag for H_F_MDS, H_CORR, H_CORR25, and H_CORR75"
!  H_CORR,,"Sensible heat flux, corrected H_F_MDS by energy balance closure correction factor"
!  H_CORR_25,,"Sensible heat flux, corrected H_F_MDS by energy balance closure correction factor, 25th percentile"
!  H_CORR_75,,"Sensible heat flux, corrected H_F_MDS by energy balance closure correction factor, 75th percentile"
!  H_RANDUNC,,"Random uncertainty of H, from measured only data"
!  H_RANDUNC_METHOD,,Method used to estimate the random uncertainty of H
!  H_RANDUNC_N,,Number of half-hour data points used to estimate the random uncertainty of H
!  H_CORR_JOINTUNC,,Joint uncertainty estimation for H
!  EBC_CF_N,,Number of data points used to calculate energy closure balance correction factor. Driver data points within sliding window (ECB_CF Method 1) or number of ECB_CF data points (for ECB_CF Methods 2 and 3)
!  EBC_CF_METHOD,,Method used to calculate the energy balance closure correction factor
! 
! !FILES USED:
!
! !REVISION HISTORY:
! Sujay Kumar; Initial version
! 
!EOP
!----------------------------------------------------------------------
  implicit none

  integer,         intent(in)      :: source
  integer,         intent(in)      :: stn_index

  real :: Month     ! Month of current data
  real :: Day       ! Day of current data
  real :: Hour      ! Hour, from 0 to 23.5
  real :: min
  real :: sec
  real :: DoY       ! Decimal day of the year

  ! Initializes other variables. iCurrent is the location of the comma after
  ! the field currently being read. 
  ! i is the index of the arrays being accessed.
  Integer :: iCurrent       
  Integer :: year, arrayLen, stat
  integer :: ftn 
  integer :: n
  integer :: status

  !This initializes the current line, first line, and filename variables. 
  Character (len = 2000) :: currentLine
  Character (len = 500)  :: filename 
  character*4           :: fyr
  logical               :: file_exists

  type(ESMF_Time)       :: fluxnettime
  real                  :: gmt
  integer               :: zone
  integer               :: ios

  character*20          :: startTime,stopTime
  real                  :: TA_F_MDS,TA_F_MDS_QC,TA_ERA,TA_F,TA_F_QC,SW_IN_POT
  real                  :: SW_IN_F_MDS,SW_IN_F_MDS_QC,SW_IN_ERA,SW_IN_F,SW_IN_F_QC
  real                  :: LW_IN_F_MDS,LW_IN_F_MDS_QC,LW_IN_ERA,LW_IN_F,LW_IN_F_QC
  real                  :: LW_IN_JSB,LW_IN_JSB_QC,LW_IN_JSB_ERA,LW_IN_JSB_F,LW_IN_JSB_F_QC
  real                  :: VPD_F_MDS,VPD_F_MDS_QC,VPD_ERA,VPD_F,VPD_F_QC
  real                  :: PA,PA_ERA,PA_F,PA_F_QC,P,P_ERA,P_F,P_F_QC
  real                  :: WS,WS_ERA,WS_F,WS_F_QC,CO2_F_MDS,CO2_F_MDS_QC
  real                  :: TS_F_MDS_1,TS_F_MDS_1_QC,SWC_F_MDS_1,SWC_F_MDS_1_QC
  real                  :: LE_F_MDS,LE_CORR,LE_CORR_25,LE_CORR_75,LE_F_MDS_QC
  real                  :: LE_RANDUNC,LE_RANDUNC_METHOD,LE_RANDUNC_N,LE_CORR_JOINTUNC
  real                  :: H_F_MDS,H_CORR,H_CORR_25,H_CORR_75,H_F_MDS_QC,H_RANDUNC
  real                  :: H_RANDUNC_METHOD,H_RANDUNC_N,H_CORR_JOINTUNC,EBC_CF_N
  real                  :: EBC_CF_METHOD

  write(unit=fyr,fmt='(i4.4)') LVT_rc%dyr(source)

  filename = trim(FLUXNET2015obs(source)%odir)//'/'//trim(fyr)//'/FLX_'//&
       trim(FLUXNET2015obs(source)%stn_name(stn_index))&
       //'_FLUXNET2015_FULLSET_HH_'//trim(fyr)//'.dat'

  inquire(file=trim(filename),exist=file_exists)

  if(file_exists) then 
      write(LVT_logunit,*) '[INFO] Reading ',trim(filename)
      year = LVT_rc%dyr(source)
  
      !Calculates whether or not the current year is a leap year following the
      !method outlined by Microsoft. (Since the data is all recent, this 
      !method could simply be replaced with Mod(year, 4) == 0 but this is far
      !more accurate)
      !The value 17520 is 356 days * 24 hours * 2 measurements per hour
      !Likewise, 17568 is 366 days * 24 hours * 2 measurements per hour
      if (Mod(year, 4) == 0) Then            !Step1
         if (Mod(year, 100) == 0) Then       !Step2
            if (Mod(year, 400) == 0) Then    !Step3
               arrayLen = 17568              !Step4, leap year
            else                          
               arrayLen = 17520              !Step5, not leap year
            end if
         else
            arrayLen = 17568                 !Step4, leap year
         end if
      else
         arrayLen = 17520                    !Step5, not leap year
      End if
    
      !Opens the file and promptly discards the header
      ftn = LVT_getNextUnitNumber()
      Open (ftn, file = filename)
      READ (ftn, '(a)') currentLine
  
      ! Repeatedly reads a line and processes it
      DO         
         READ (ftn, "(A2000)", iostat = stat) currentLine
         if (stat /= 0) exit       
	
         ! looks for the first available comma and indexes it
         ! reads what is between the start of the line and the
         ! comma as the next value sets the current line to the
         ! current line minus that value and its comma. Values 
         ! that were not used are not read into any array but 
         ! are simply discarded. Basically the entire idea 
         ! revolves around the format of the data:
         ! value1,value2,value3,etc


         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *,iostat=ios) startTime
         if(ios.ne.0) exit
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *,iostat=ios) stopTime
         if(ios.ne.0) exit
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
          
         !air temperature gap filled using MDS method
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) ta_f_mds
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         !quality flag for ta_f_mds
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) TA_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) TA_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
    
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) TA_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) TA_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SW_IN_POT
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SW_IN_F_MDS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SW_IN_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SW_IN_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SW_IN_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SW_IN_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_F_MDS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_JSB
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_JSB_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_JSB_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_JSB_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LW_IN_JSB_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) VPD_F_MDS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) VPD_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) VPD_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) VPD_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) VPD_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) PA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) PA_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) PA_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) PA_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) P
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) P_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) P_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) P_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) WS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) WS_ERA
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) WS_F
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) WS_F_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) CO2_F_MDS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) CO2_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) TS_F_MDS_1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) TS_F_MDS_1_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SWC_F_MDS_1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) SWC_F_MDS_1_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_F_MDS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_CORR
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_CORR_25
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_CORR_75
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_RANDUNC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_RANDUNC_METHOD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_RANDUNC_N
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) LE_CORR_JOINTUNC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_F_MDS
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_CORR
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_CORR_25
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_CORR_75
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_F_MDS_QC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_RANDUNC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_RANDUNC_METHOD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_RANDUNC_N
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) H_CORR_JOINTUNC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) EBC_CF_N
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         iCurrent = Index(currentLine, ",") 
         Read(currentLine(1: iCurrent - 1), *) EBC_CF_METHOD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     

!         read(stopTime,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') year, month, day, hour,min
         read(stopTime(1:4),*) year 
         stopTime = stopTime(5:len(stopTime))

         read(stopTime(1:2),*) month
         stopTime = stopTime(3:len(stopTime))

         read(stopTime(1:2),*) day
         stopTime = stopTime(3:len(stopTime))

         read(stopTime(1:2),*) hour
         stopTime = stopTime(3:len(stopTime))

         read(stopTime(1:2),*) min
         stopTime = stopTime(3:len(stopTime))
         

         call ESMF_TimeSet(fluxnettime, yy=year, &
              mm = nint(month), &
              dd = nint(day), &
              h = floor(hour), &
              m = nint(min), &
              calendar = LVT_calendar, &
              rc=status)
         call LVT_verify(status)
         
         n = nint((fluxnettime - FLUXNET2015obs(source)%starttime)/FLUXNET2015obs(source)%timestep) + 1

         if(n.ge.1.and.n.le.17570) then 

            if (LE_F_MDS.ne.-9999.0000.and.LE_F_MDS_QC.eq.0) Then
               FLUXNET2015obs(source)%Qle(stn_index, n) = LE_F_MDS
            end if

            if (H_F_MDS.ne.-9999.0000.and.H_F_MDS_QC.eq.0) Then
               FLUXNET2015obs(source)%Qh(stn_index, n) = H_F_MDS
            end if

         endif
      END DO

      call LVT_releaseUnitNumber(ftn)

      write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)
   end if

end subroutine read_fluxnet_station
