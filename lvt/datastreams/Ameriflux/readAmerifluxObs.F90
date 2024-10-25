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
! 
! !ROUTINE: readtemplateObs
! \label{readtemplateObs}
!
! !INTERFACE: 
  subroutine readAmerifluxObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use LVT_logMod
  use Ameriflux_obsMod
  use map_utils
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This initializes variables the entire program will use such as yr, DoY, and 
! status.It calls the read_Ameriflux_Station routine if one or both of the 
! following criteria are met: a flag is set to true, or the LVT year and data
! input year do not match. The routine then time averates the data and converts
! lat, long coordinates to i, j and logs the data. Specifically, the data used
! was Level 3 (L3). This was chosen out of the four available levels (L1, L2, 
! L3, and L4) because L1 was raw data that could not be used for our purposes.
! L2 has the same general values but without the quality flags that could help
! us. Meanwhile L4 filters out some of the fields which are most important to 
! us such as soil heat flux. The data itself belongs to Ameriflux and was used
! with their permission. To see the original data, go to:
! ftp://cdiac.ornl.gov/pub/ameriflux/data/Level3/
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Jun 2010: Teodor Georgiev, Ameriflux reader
! 
!EOP
!--------------------------------------------------------------------------
  implicit none

  integer         :: source 

  integer :: i 
  type(ESMF_Time) :: amerifluxtime1, amerifluxtime2
  integer         :: t, st, et
  integer         :: yr, mo, da, hr, mn, ss
  real*8          :: lis_prevtime
  integer         :: c,r,stn_col, stn_row
  real            :: col,row
  real            :: gmt
  integer         :: doy
  real            :: qle(LVT_rc%lnc, LVT_rc%lnr)
  real            :: qh(LVT_rc%lnc, LVT_rc%lnr)
  real            :: qg(LVT_rc%lnc, LVT_rc%lnr)
  real            :: ta(LVT_rc%lnc, LVT_rc%lnr)
  real            :: sm1(LVT_rc%lnc, LVT_rc%lnr)
  real            :: st1(LVT_rc%lnc, LVT_rc%lnr)
  real            :: pcp(LVT_rc%lnc, LVT_rc%lnr)
  integer         :: status
  real                  :: timenow

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)

  qle = LVT_rc%udef
  qh = LVT_rc%udef
  qg = LVT_rc%udef
  ta = LVT_rc%udef
  sm1 =LVT_rc%udef
  st1 =LVT_rc%udef
  pcp =LVT_rc%udef

  ! Start of new data
  if(Amerifluxobs(source)%startflag.or.&
       (Amerifluxobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(Amerifluxobs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm=1, dd=1, h=0, m=0,s = 0, calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'ameriflux starttime set failed')

     Amerifluxobs(source)%startflag = .false. 
     Amerifluxobs(source)%yr = LVT_rc%dyr(source)

     do i = 1,Amerifluxobs(source)%n_stns
        call read_ameriflux_station(source,i)
     end do
  endif

  if(mod(timenow, 1800.0).eq.0.0) then 
     call ESMF_TimeSet(amerifluxtime1, yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), &
          m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'amerifluxtime1 set failed')
     
     t = nint((amerifluxtime1 - amerifluxobs(source)%starttime)/&
          amerifluxobs(source)%timestep)+1

     do i=1,amerifluxobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, amerifluxobs(source)%stnlat(i), &
             amerifluxobs(source)%stnlon(i), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if((stn_col.ge.1.and.stn_col.le.LVT_rc%lnc).and.&
             (stn_row.ge.1.and.stn_row.le.LVT_rc%lnr)) then 
           if(amerifluxobs(source)%qle(i,t).ne.LVT_rc%udef) then 
              qle(stn_col, stn_row) =  amerifluxobs(source)%qle(i,t)
           end if
           
           if(amerifluxobs(source)%qh(i,t).ne.LVT_rc%udef) then 
              qh(stn_col, stn_row) = amerifluxobs(source)%qh(i,t)
           end if
        
           if(amerifluxobs(source)%qg(i, t).ne.LVT_rc%udef) then
              qg(stn_col, stn_row) = amerifluxobs(source)%qg(i, t)
           end if
           
           if(amerifluxobs(source)%ta(i, t).ne.LVT_rc%udef) then
              ta(stn_col, stn_row) = amerifluxobs(source)%ta(i, t)
           end if
           
           if(amerifluxobs(source)%sfsm(i, t).ne.LVT_rc%udef) then
              sm1(stn_col, stn_row) = amerifluxobs(source)%sfsm(i, t)
           end if
           
           if(amerifluxobs(source)%sfst(i, t).ne.LVT_rc%udef) then
              st1(stn_col, stn_row) = amerifluxobs(source)%sfst(i, t)
           end if
           if(amerifluxobs(source)%precip(i, t).ne.LVT_rc%udef) then
              pcp(stn_col, stn_row) = amerifluxobs(source)%precip(i, t)
           end if
        endif
     end do
  else
     qle = LVT_rc%udef
     qh = LVT_rc%udef
     qg = LVT_rc%udef
     ta = LVT_rc%udef
     sm1 = LVT_rc%udef
     st1 = LVT_rc%udef
     pcp = LVT_rc%udef
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_qle, source, qle,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qh, source,qh,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qg, source,qg,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source,sm1, vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp, source,st1, vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_tairforc, source,ta,vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source,pcp,vlevel=1,units="kg/m2")

end subroutine readAmerifluxObs

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
subroutine read_ameriflux_station(source, stn_index)
! 
! !USES:
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use Ameriflux_obsMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! Opens the file and reads line by line, inputting each variable into its array
! Further information on reading is provided at the code
! The array will then be processed in the above subroutine
! 
! !FILES USED:
!
! !REVISION HISTORY:
! Teodor Georgiev, initial version
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
  real :: Le        ! temporary latent heat flux 
  real :: H         ! temporary sensible heat flux
  real :: Precip    ! Precipitation in mm
  real :: Rg        ! global radiation W/m2
  real :: qf_Rg     ! quality flag for Rg:
                                                ! 0: ok, 1: possible spike
  real :: qf_Rad    ! quality flag for Rg:
                                                ! 0: ok, 1: possible bad values
  real :: Rn        ! net radiation W/m2
  real :: Rd        ! diffuse radiation W/m2
  real :: Ta        ! air temperature C
  real :: Ts1       ! soil temperature depth 1 C
  real :: Ts2       ! soil temperature depth 2 C
  real :: SWC1      ! soil water content depth 1 %vol
  real :: SWC2      ! soil water content depth 2 %vol
  real :: G1	! Soil heat flux
  real :: G2	! Soil heat flux 2
  real :: Rh        ! relative humidity as a %
  real :: WS        ! Wind horizontal speed m/s
  ! Initializes other variables. iCurrent is the location of the comma after
  ! the field currently being read. 
  ! i is the index of the arrays being accessed.
  Integer :: iCurrent, i          
  Integer :: year, arrayLen, stat
  integer :: ftn 
  integer :: n
  integer :: status

  !This initializes the current line, first line, and filename variables. 
  Character (len = 300) :: currentLine
  Character (len = 200) :: filename 
  character*4           :: fyr
  logical               :: file_exists

  type(ESMF_Time)       :: amerifluxtime
  real                  :: gmt
  integer               :: zone
  integer               :: ios

  write(unit=fyr,fmt='(i4.4)') LVT_rc%dyr(source)
  if(Amerifluxobs(source)%version.eq."Level3") then 
     filename = trim(Amerifluxobs(source)%odir)//'/'//&
          trim(Amerifluxobs(source)%site_name(stn_index))&
          //'/'//trim(adjustl(Amerifluxobs(source)%stn_name(stn_index)))//&
          trim(fyr)//'_L3.txt'
  else
     write(LVT_logunit,*) '[ERR] This version of the Ameriflux data is not currently supported'
     call LVT_endrun()
     
  endif
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
      READ (ftn, *)
  
      i = 1

      ! Repeatedly reads a line and processes it
      DO         
         READ (ftn, "(A300)", iostat = stat) currentLine
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
         READ(currentLine(1: iCurrent - 1), *,iostat=ios) Month
         if(ios.ne.0) exit
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *,iostat=ios) Day
         if(ios.ne.0) exit
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *,iostat=ios) Hour
         if(ios.ne.0) exit
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *,iostat=ios) DoY
         if(ios.ne.0) exit
         currentLine = currentLine(iCurrent + 1: Len(currentLine)) 
     
         iCurrent = Index(currentLine, ",") !CO2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !H2O
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !ZL
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !FC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_Fc
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !NEE_st
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_NEE_st
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !NEE_or
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_NEE_or
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) H !Amerifluxobs(source)%Qh(stn_index,i)
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Le !Amerifluxobs(source)%Qle(stn_index,i)
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !ustar
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_ust
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Precip
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rg
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !PPFD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !R_pot
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) qf_Rg
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) qf_Rad
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !Rr
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rn
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rd
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !APAR
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Ta
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Ts1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Ts2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) SWC1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) SWC2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) G1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) G2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rh
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !WD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         ! by the last line the current line consists of only a value,
         ! no comma indexing is required
         Read(currentLine, *) WS
         
         min = (hour - floor(hour))*60.0
         sec = (min- floor(min))*60.0
         
         call LVT_localtime2gmt(gmt, Amerifluxobs(source)%stnlon(stn_index), hour, zone)         
         call ESMF_TimeSet(amerifluxtime, yy=year, &
              mm = nint(month), &
              dd = nint(day), &
              h = floor(gmt), &
              m = nint(min), &
              s = nint(sec), & 
              calendar = LVT_calendar, &
              rc=status)
         call LVT_verify(status)
         
         n = nint((amerifluxtime - Amerifluxobs(source)%starttime)/Amerifluxobs(source)%timestep) + 1
!         print*, year, month, day, floor(gmt), min, sec, n
         
         if(n.ge.1.and.n.le.17570) then 
            if (H.ne.-9999.0000) Then
               Amerifluxobs(source)%Qh(stn_index, n) = H
            end if
            
            if (Le.ne.-9999.0000) Then
               Amerifluxobs(source)%Qle(stn_index, n) = Le
            end if
        
         !The current method only logs G if both values are present
         !This is sometimes ineffective since one of the two is
         !usually inaccurate, so I have included the code to accept
         !either one or the other.
            if ((G1.ne.-9999.0000).and.(G2.ne.-9999.0000)) Then
               Amerifluxobs(source)%Qg(stn_index, n) = (G1 + G2)/2
        ! else if (G1.ne.-9999.0000) Then
        !    Amerifluxobs(source)%Qg(stn_index, n) = G1   
        ! else if (G2.ne.-9999.0000) Then
        !    Amerifluxobs(source)%Qg(stn_index, n) = G2            
            end if
         
            if(Ta.ne.-9999.0000) then 
               Amerifluxobs(source)%Ta(stn_index,n) = Ta + 273.15
            endif

            if(Amerifluxobs(source)%nstlayers.eq.1) then             
               if (ts1.ne.-9999.0000) then 
                  Amerifluxobs(source)%sfst(stn_index,n) = &
                       (Amerifluxobs(source)%sfst_wt(stn_index,1)*(ts1+273.15))
               endif
            elseif(Amerifluxobs(source)%nstlayers.eq.2) then                
               if ((ts1.ne.-9999.0000.and.ts2.ne.-9999.0)) then 
                  Amerifluxobs(source)%sfst(stn_index,n) = &
                       (Amerifluxobs(source)%sfst_wt(stn_index,1)*(ts1+273.15)+&
                       Amerifluxobs(source)%sfst_wt(stn_index,2)*(ts2+273.15) )
               endif
            endif

            if(Amerifluxobs(source)%nsmlayers.eq.1) then             
               if (swc1.ne.-9999.0000.and.swc1/100.0.le.0.5) then 
                  Amerifluxobs(source)%sfsm(stn_index,n) = &
                       (Amerifluxobs(source)%sfsm_wt(stn_index,1)*swc1/100.0)
                  if(stn_index.eq.40) print*, 'sm1',Amerifluxobs(source)%sfsm(stn_index,n)
                  if(Amerifluxobs(source)%sfsm(stn_index,n).gt.0.5) then 
                     print*, 'Note that soil moisture > 0.5 '
                     print*, ' -- likely due to the reporting being '
                     print*, 'in percentages. We need to multiply it '
                     print*, 'with the porosity for this location'
                     print*, stn_index,Amerifluxobs(source)%stn_name(stn_index), &
                       Amerifluxobs(source)%sfsm(stn_index,n), swc1                  
                     stop
                  endif

               endif
            elseif(Amerifluxobs(source)%nsmlayers.eq.2) then                
               if ((swc1.ne.-9999.0000.and.swc2.ne.-9999.0).and.&
                    (swc1/100.0.le.0.5.and.swc2/100.0.le.0.5))then 
                  Amerifluxobs(source)%sfsm(stn_index,n) = &
                       (Amerifluxobs(source)%sfsm_wt(stn_index,1)*swc1/100.0+&
                       Amerifluxobs(source)%sfsm_wt(stn_index,2)*swc2/100.0 )
                  if(stn_index.eq.40) print*, 'sm2',Amerifluxobs(source)%sfsm(stn_index,n)
                  if(Amerifluxobs(source)%sfsm(stn_index,n).gt.0.5) then 
                     print*, 'Note that soil moisture > 0.5 '
                     print*, ' -- likely due to the reporting being '
                     print*, 'in percentages. We need to multiply it '
                     print*, 'with the porosity for this location'
                     print*, stn_index,n,Amerifluxobs(source)%sfsm(stn_index,n), &
                          swc1, swc2, Amerifluxobs(source)%sfsm_wt(stn_index,1), &
                          Amerifluxobs(source)%sfsm_wt(stn_index,2)
                  endif
               endif
            endif
         
            if (precip.ne.-9999.0000) Then
               Amerifluxobs(source)%precip(stn_index, n) = precip/1800.0 !mm/s
            end if
            
            i= i + 1
         endif
      END DO

      call LVT_releaseUnitNumber(ftn)

      write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)
   end if

end subroutine read_ameriflux_station
