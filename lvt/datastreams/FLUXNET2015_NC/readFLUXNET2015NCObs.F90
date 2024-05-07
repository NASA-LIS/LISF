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
  subroutine readFLUXNET2015NCObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use LVT_logMod
  use FLUXNET2015NC_obsMod
  use map_utils
!
! !DESCRIPTION: 
! 
!   This routine provides the data processing code for the FLUXNET2015NC
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
  real            :: qtau(LVT_rc%lnc, LVT_rc%lnr)
  real            :: gpp(LVT_rc%lnc, LVT_rc%lnr)
  integer         :: status
  real            :: timenow

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)

  qle = LVT_rc%udef
  qh = LVT_rc%udef
  qtau = LVT_rc%udef
  gpp = LVT_rc%udef

  ! Start of new data
  if(FLUXNET2015NCobs(source)%startflag.or.&
       (FLUXNET2015NCobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 

     call ESMF_TimeSet(FLUXNET2015NCobs(source)%starttime, &
          yy=LVT_rc%dyr(source), &
          mm=1, dd=1, h=0, m=0,s = 0, calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'fluxnet starttime set failed')

     FLUXNET2015NCobs(source)%startflag = .false. 
     FLUXNET2015NCobs(source)%yr = LVT_rc%dyr(source)

     FLUXNET2015NCobs(source)%qle(:,:) = LVT_rc%udef
     FLUXNET2015NCobs(source)%qh(:,:) = LVT_rc%udef
     FLUXNET2015NCobs(source)%qtau(:,:) = LVT_rc%udef
     FLUXNET2015NCobs(source)%gpp(:,:) = LVT_rc%udef

     do i = 1,FLUXNET2015NCobs(source)%n_stns
        call read_fluxnet_station_nc(source,i)
     end do
  endif

  if(mod(timenow, 1800.0).eq.0.0) then 
     call ESMF_TimeSet(fluxnettime1, yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), &
          m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'fluxnettime1 set failed')
     
     t = nint((fluxnettime1 - FLUXNET2015NCobs(source)%starttime)/&
          FLUXNET2015NCobs(source)%timestep)+1

     do i=1,FLUXNET2015NCobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, FLUXNET2015NCobs(source)%stnlat(i), &
             FLUXNET2015NCobs(source)%stnlon(i), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if((stn_col.ge.1.and.stn_col.le.LVT_rc%lnc).and.&
             (stn_row.ge.1.and.stn_row.le.LVT_rc%lnr)) then 
           if(FLUXNET2015NCobs(source)%qle(i,t).ne.LVT_rc%udef) then 
              qle(stn_col, stn_row) =  FLUXNET2015NCobs(source)%qle(i,t)
           end if
           
           if(FLUXNET2015NCobs(source)%qh(i,t).ne.LVT_rc%udef) then 
              qh(stn_col, stn_row) = FLUXNET2015NCobs(source)%qh(i,t)
           end if

           if(FLUXNET2015NCobs(source)%qtau(i,t).ne.LVT_rc%udef) then 
              qtau(stn_col, stn_row) = FLUXNET2015NCobs(source)%qtau(i,t)
           end if

           if(FLUXNET2015NCobs(source)%gpp(i,t).ne.LVT_rc%udef) then 
              gpp(stn_col, stn_row) =  FLUXNET2015NCobs(source)%gpp(i,t)
           end if

        
        endif
     end do
  else
     qle = LVT_rc%udef
     qh = LVT_rc%udef
     qtau = LVT_rc%udef
     gpp = LVT_rc%udef

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_qle, source, qle,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qh, source,qh,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qtau, source,qtau,vlevel=1,units="m/s")
  call LVT_logSingleDataStreamVar(LVT_MOC_gpp, source, gpp,vlevel=1,units="umol/m2s")


end subroutine readFLUXNET2015NCObs

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
subroutine read_fluxnet_station_nc(source, stn_index)
! 
! !USES:
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use FLUXNET2015NC_obsMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

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
  integer               :: nid, tid, dim1Id, QLEid, QLEQCid, QHid, QHQCid, QTAUid, GPPId, tdims, t
  integer               :: yr, mo, da, hr, mn
  integer               :: startIndex, endIndex
  real*8, allocatable     :: timestamp(:)
  integer*8, allocatable  :: timestamp1(:)
  real, allocatable     :: QLE(:),QLEQC(:)
  real, allocatable     :: QH(:),QHQC(:)
  real, allocatable     :: QTAU(:)
  real, allocatable     :: GPP(:)
  character*12          :: date
  character*20          :: startTime,stopTime

  write(unit=fyr,fmt='(i4.4)') LVT_rc%dyr(source)

  filename = trim(FLUXNET2015NCobs(source)%odir)//'/FLX_'//&
       trim(FLUXNET2015NCobs(source)%stn_name(stn_index))&
       //'_FLUXNET2015_SUBSET_HH_LVT.nc'
  inquire(file=trim(filename),exist=file_exists)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

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
    
      ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
      call LVT_verify(ios, 'Error opening file'//trim(filename))  

!dimensions
      ios = nf90_inq_dimid(nid, 'time',dim1Id)
      call LVT_verify(ios, 'Error nf90_inq_dimid: time')

      ios = nf90_inquire_dimension(nid, dim1Id, len=tdims)
      call LVT_verify(ios, 'Error nf90_inquire_dimension:')

      allocate(timestamp(tdims))
      allocate(timestamp1(tdims))

      ios = nf90_inq_varid(nid, 'TIMESTAMP_START',tid)
      call LVT_verify(ios, 'Error nf90_inq_varid: TIMESTAMP_START')

      ios = nf90_get_var(nid,tid, timestamp)
      call LVT_verify(ios, 'Error nf90_get_var: timestamp')

      startIndex = -1
      endIndex = -1
      do t=1,tdims
         timestamp1(t) = timestamp(t)
         write(unit=date, fmt='(i12.12)') timestamp1(t)
         read(date, fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') yr, mo, da, hr, mn

         if(LVT_rc%dyr(source).eq.yr.and.&
              mo.eq.1.and.&
              da.eq.1.and.&
              hr.eq.0.and.&
              mn.eq.0) then 
            startIndex = t            
         endif
         if((LVT_rc%dyr(source)+1).eq.yr.and.&
              mo.eq.1.and.&
              da.eq.1.and.&
              hr.eq.0.and.&
              mn.eq.0) then 
            endIndex = t
            exit;
         endif
      enddo
      if(endIndex.eq.-1) then
         endIndex = tdims
      endif
      if(startIndex.ne.-1.and.&
           endIndex.ge.startIndex) then 
!------------------------------------------------------------------
! latent heat flux 
!------------------------------------------------------------------
         allocate(QLE(endIndex-startIndex+1))
         allocate(QLEQC(endIndex-startIndex+1))

         QLE = LVT_rc%udef
         QLEQC = LVT_rc%udef

         !Get LE
         ios = nf90_inq_varid(nid, 'LE_F_MDS',QLEid)
         call LVT_verify(ios, 'Error nf90_inq_varid: QLE')
         
         ios = nf90_get_var(nid,QLEid, QLE,  start=(/startIndex/), &
              count=(/endIndex-startIndex+1/))
         call LVT_verify(ios, 'Error nf90_get_var: QLE')

         !Get LE QC
         ios = nf90_inq_varid(nid, 'LE_F_MDS_QC',QLEQCid)
         call LVT_verify(ios, 'Error nf90_inq_varid: QLEQC')
         
         ios = nf90_get_var(nid,QLEQCid, QLEQC,  start=(/startIndex/), &
              count=(/endIndex-startIndex+1/))
         call LVT_verify(ios, 'Error nf90_get_var: QLEQC')

         FLUXNET2015NCobs(source)%Qle(stn_index,:) = LVT_rc%udef
         do t=1,endIndex-startIndex+1
            if (QLE(t).ne.-9999.0000.and.QLEQC(t).eq.0) then
                FLUXNET2015NCobs(source)%Qle(stn_index,t) = &
                      QLE(t)
            endif
         enddo
!------------------------------------------------------------------
! sensible heat flux 
!------------------------------------------------------------------
         allocate(QH(endIndex-startIndex+1))
         allocate(QHQC(endIndex-startIndex+1))

         QH = LVT_rc%udef
         QHQC = LVT_rc%udef

         !Get H 
         ios = nf90_inq_varid(nid, 'H_F_MDS',QHid)
         call LVT_verify(ios, 'Error nf90_inq_varid: QH')
         
         ios = nf90_get_var(nid,QHid, QH,  start=(/startIndex/), &
              count=(/endIndex-startIndex+1/))
         call LVT_verify(ios, 'Error nf90_get_var: QH')

         !Get HQC 
         ios = nf90_inq_varid(nid, 'H_F_MDS_QC',QHQCid)
         call LVT_verify(ios, 'Error nf90_inq_varid: QHQC')
         
         ios = nf90_get_var(nid,QHQCid, QHQC,  start=(/startIndex/), &
              count=(/endIndex-startIndex+1/))
         call LVT_verify(ios, 'Error nf90_get_var: QHQC')

         FLUXNET2015NCobs(source)%Qh(stn_index,:) = LVT_rc%udef
         do t=1,endIndex-startIndex+1
            if (QH(t).ne.-9999.0000.and.QHQC(t).eq.0) then
                FLUXNET2015NCobs(source)%Qh(stn_index,t) = &
                    QH(t)
            endif
         enddo
!------------------------------------------------------------------
! Momentum flux (friction velocity) 
!------------------------------------------------------------------
         allocate(QTAU(endIndex-startIndex+1))

         QTAU = LVT_rc%udef

         !Get Ustar
         ios = nf90_inq_varid(nid, 'USTAR',QTAUid)
         call LVT_verify(ios, 'Error nf90_inq_varid: QTAU')
         
         ios = nf90_get_var(nid,QTAUid, QTAU,  start=(/startIndex/), &
              count=(/endIndex-startIndex+1/))
         call LVT_verify(ios, 'Error nf90_get_var: QTAU')

         FLUXNET2015NCobs(source)%Qtau(stn_index,:) = LVT_rc%udef
         do t=1,endIndex-startIndex+1
            if (QTAU(t).ne.-9999.0000) then
                FLUXNET2015NCobs(source)%Qtau(stn_index,t) = &
                      QTAU(t)
            endif
         enddo

!------------------------------------------------------------------
! GPP 
!------------------------------------------------------------------
         allocate(GPP(endIndex-startIndex+1))

         GPP = LVT_rc%udef

         !Get GPP
         ios = nf90_inq_varid(nid, 'GPP_NT_VUT_REF',GPPid)
         call LVT_verify(ios, 'Error nf90_inq_varid: GPP')
         
         ios = nf90_get_var(nid,GPPid, GPP,  start=(/startIndex/), &
              count=(/endIndex-startIndex+1/))
         call LVT_verify(ios, 'Error nf90_get_var: GPP')
       
         FLUXNET2015NCobs(source)%GPP(stn_index,:) = LVT_rc%udef
         do t=1,endIndex-startIndex+1
            
            if (GPP(t).ne.-9999.0000) then
                FLUXNET2015NCobs(source)%GPP(stn_index,t) = &
                      GPP(t)
            endif
         enddo
!-

!------------------------------------------------------------------

         deallocate(QLE)
         deallocate(QLEQC)
         deallocate(QH)
         deallocate(QHQC)
         deallocate(QTAU)
         deallocate(GPP)

      endif
      ios = nf90_close(nid)
      deallocate(timestamp)
      deallocate(timestamp1)

      write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)
   end if
#endif

 end subroutine read_fluxnet_station_nc
